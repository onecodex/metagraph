#ifndef __VECTOR_ALGORITHM_HPP__
#define __VECTOR_ALGORITHM_HPP__

#include <functional>
#include <cassert>
#include <cstdlib>

#include <sdsl/int_vector.hpp>
#include <sdsl/select_support_scan.hpp>

class ThreadPool;
class bit_vector;


sdsl::bit_vector to_sdsl(const std::vector<bool> &vector);
sdsl::bit_vector to_sdsl(const std::vector<uint8_t> &vector);

template <class Vector>
sdsl::int_vector<> pack_vector(const Vector &vector, uint8_t bits_per_number) {
    if constexpr(std::is_same_v<Vector, sdsl::int_vector<>>) {
        if (bits_per_number == vector.width())
            return vector;
    }

    sdsl::int_vector<> packed(vector.size(), 0, bits_per_number);
    for (uint64_t i = 0; i < vector.size(); ++i) {
        packed[i] = vector[i];
    }
    return packed;
}

sdsl::int_vector<> pack_vector(sdsl::int_vector<>&& vector,
                               uint8_t bits_per_number);


/**
 * Atomic bit fetching, setting, and unsetting on packed vectors.
 * fetch_and_* return the old values. The default memorder __ATOMIC_SEQ_CST
 * enforces the ordering of writes and reads across threads. See
 * https://gcc.gnu.org/onlinedocs/gcc/_005f_005fatomic-Builtins.html
 * for more details.
 */

inline bool fetch_and_set_bit(uint64_t *v,
                              uint64_t i,
                              bool atomic = false,
                              int mo = __ATOMIC_SEQ_CST) {
    const uint64_t mask = (1llu << (i & 0x3F));

    if (atomic) {
        return __atomic_fetch_or(&v[i >> 6], mask, mo) & mask;
    } else {
        uint64_t &word = v[i >> 6];
        if (word & mask) {
            return true;
        } else {
            word |= mask;
            return false;
        }
    }
}

inline bool fetch_and_unset_bit(uint64_t *v,
                                uint64_t i,
                                bool atomic = false,
                                int mo = __ATOMIC_SEQ_CST) {
    const uint64_t mask = (1llu << (i & 0x3F));

    if (atomic) {
        return __atomic_fetch_and(&v[i >> 6], ~mask, mo) & mask;
    } else {
        uint64_t &word = v[i >> 6];
        if (word & mask) {
            word &= ~mask;
            return true;
        } else {
            return false;
        }
    }
}

inline bool fetch_bit(const uint64_t *v,
                      uint64_t i,
                      bool atomic = false,
                      int mo = __ATOMIC_SEQ_CST) {
    return atomic
        ? ((__atomic_load_n(&v[i >> 6], mo) >> (i & 0x3F)) & 1)
        : ((v[i >> 6] >> (i & 0x3F)) & 1);
}

inline void set_bit(uint64_t *v,
                    uint64_t i,
                    bool atomic = false,
                    int mo = __ATOMIC_SEQ_CST) {
    if (atomic) {
        __atomic_or_fetch(&v[i >> 6], 1llu << (i & 0x3F), mo);
    } else {
        v[i >> 6] |= (1llu << (i & 0x3F));
    }
}

inline void unset_bit(uint64_t *v,
                      uint64_t i,
                      bool atomic = false,
                      int mo = __ATOMIC_SEQ_CST) {
    if (atomic) {
        __atomic_and_fetch(&v[i >> 6], ~(1llu << (i & 0x3F)), mo);
    } else {
        v[i >> 6] &= ~(1llu << (i & 0x3F));
    }
}

inline uint64_t atomic_fetch_and_add(sdsl::int_vector<> &vector, uint64_t i,
                                     uint64_t count,
                                     std::mutex &backup_mutex) {
    // TODO: support adding negatives
#ifdef MODE_TI
    size_t width = vector.width();
    uint64_t bit_pos = i * width;
    uint8_t shift = bit_pos & 0x7F;

    if (shift + width <= 128) {
        __uint128_t *limb = &reinterpret_cast<__uint128_t*>(vector.data())[bit_pos >> 7];
        __uint128_t mask = ((__uint128_t(1) << width) - 1) << shift;
        __uint128_t inc = __uint128_t(count) << shift;
        __uint128_t exp_val;
        __uint128_t new_val;
        do {
            exp_val = *limb;
            new_val = ((exp_val + inc) & mask) | (exp_val & (~mask));
        } while (!__sync_bool_compare_and_swap(limb, exp_val, new_val));

        return (exp_val & mask) >> shift;
    }

#endif

    std::lock_guard<std::mutex> lock(backup_mutex);
    std::atomic_thread_fence(std::memory_order_acquire);
    uint64_t old_val = vector[i];
    vector[i] += count;
    std::atomic_thread_fence(std::memory_order_release);
    return old_val;
}

inline uint64_t atomic_exchange(sdsl::int_vector<> &vector, uint64_t i, uint64_t val,
                                std::mutex &backup_mutex) {
#ifdef MODE_TI
    size_t width = vector.width();
    uint64_t bit_pos = i * width;
    uint8_t shift = bit_pos & 0x7F;

    if (shift + width <= 128) {
        __uint128_t *limb = &reinterpret_cast<__uint128_t*>(vector.data())[bit_pos >> 7];
        __uint128_t mask = ((__uint128_t(1) << width) - 1) << shift;
        __uint128_t val_shift = __uint128_t(val) << shift;
        __uint128_t exp_val;
        __uint128_t new_val;
        do {
            exp_val = *limb;
            new_val = val_shift | (exp_val & (~mask));
        } while (!__sync_bool_compare_and_swap(limb, exp_val, new_val));

        return (exp_val & mask) >> shift;
    }

#endif

    std::lock_guard<std::mutex> lock(backup_mutex);
    std::atomic_thread_fence(std::memory_order_acquire);
    uint64_t old_val = vector[i];
    vector[i] = val;
    std::atomic_thread_fence(std::memory_order_release);
    return old_val;
}


template <class Bitmap, class Callback>
void call_ones(const Bitmap &vector,
               uint64_t begin, uint64_t end,
               Callback callback) {
    assert(begin <= end);
    assert(end <= vector.size());

    uint64_t i = begin;
    for ( ; i < end && (i & 0x3F); ++i) {
        if (vector[i])
            callback(i);
    }
    uint64_t word;
    for (uint64_t j = i + 64; j <= end; j += 64) {
        word = vector.get_int(i, 64);
        if (!word) {
            i += 64;
            continue;
        }

        i += sdsl::bits::lo(word);
        callback(i++);

        for ( ; i < j; ++i) {
            if (vector[i])
                callback(i);
        }
    }
    for ( ; i < end; ++i) {
        if (vector[i])
            callback(i);
    }
}

template <class Callback>
void call_ones(const sdsl::bit_vector &vector,
               uint64_t begin, uint64_t end,
               Callback callback,
               bool atomic,
               int mo = __ATOMIC_SEQ_CST) {
    if (!atomic) {
        call_ones(vector, begin, end, callback);
        return;
    }

    assert(begin <= end);
    assert(end <= vector.size());

    uint64_t i = begin;
    for ( ; i < end && (i & 0x3F); ++i) {
        if (fetch_bit(vector.data(), i, true))
            callback(i);
    }
    uint64_t word;
    for (uint64_t j = i + 64; j <= end; j += 64) {
        word = __atomic_load_n(&vector.data()[i >> 6], mo);
        if (!word) {
            i += 64;
            continue;
        }

        i += sdsl::bits::lo(word);
        callback(i++);

        for ( ; i < j; ++i) {
            if (fetch_bit(vector.data(), i, true))
                callback(i);
        }
    }
    for ( ; i < end; ++i) {
        if (fetch_bit(vector.data(), i, true))
            callback(i);
    }
}

template <class Bitmap, class Callback>
void call_ones(const Bitmap &vector, Callback callback) {
    call_ones(vector, 0, vector.size(), callback);
}

template <class Callback>
void call_ones(const sdsl::bit_vector &vector,
               Callback callback,
               bool atomic,
               int mo = __ATOMIC_SEQ_CST) {
    call_ones(vector, 0, vector.size(), callback, atomic, mo);
}

template <class Bitmap, class Callback>
void call_zeros(const Bitmap &vector,
                uint64_t begin, uint64_t end,
                Callback callback) {
    assert(begin <= end);
    assert(end <= vector.size());

    uint64_t i = begin;
    for ( ; i < end && (i & 0x3F); ++i) {
        if (!vector[i])
            callback(i);
    }
    uint64_t word;
    for (uint64_t j = i + 64; j <= end; j += 64) {
        word = ~vector.get_int(i, 64);
        if (!word) {
            i += 64;
            continue;
        }

        i += sdsl::bits::lo(word);
        callback(i++);

        for ( ; i < j; ++i) {
            if (!vector[i])
                callback(i);
        }
    }
    for ( ; i < end; ++i) {
        if (!vector[i])
            callback(i);
    }
}

template <class Callback>
void call_zeros(const sdsl::bit_vector &vector,
                uint64_t begin, uint64_t end,
                Callback callback,
                bool atomic,
                int mo = __ATOMIC_SEQ_CST) {
    if (!atomic) {
        call_zeros(vector, begin, end, callback);
        return;
    }

    assert(begin <= end);
    assert(end <= vector.size());

    uint64_t i = begin;
    for ( ; i < end && (i & 0x3F); ++i) {
        if (!fetch_bit(vector.data(), i, true))
            callback(i);
    }
    uint64_t word;
    for (uint64_t j = i + 64; j <= end; j += 64) {
        word = ~__atomic_load_n(&vector.data()[i >> 6], mo);
        if (!word) {
            i += 64;
            continue;
        }

        i += sdsl::bits::lo(word);
        callback(i++);

        for ( ; i < j; ++i) {
            if (!fetch_bit(vector.data(), i, true))
                callback(i);
        }
    }
    for ( ; i < end; ++i) {
        if (!fetch_bit(vector.data(), i, true))
            callback(i);
    }
}

template <class Bitmap, class Callback>
void call_zeros(const Bitmap &vector, Callback callback) {
    call_zeros(vector, 0, vector.size(), callback);
}

template <class Callback>
void call_zeros(const sdsl::bit_vector &vector,
                Callback callback,
                bool atomic,
                int mo = __ATOMIC_SEQ_CST) {
    call_zeros(vector, 0, vector.size(), callback, atomic, mo);
}

uint64_t count_ones(const sdsl::bit_vector &vector, uint64_t begin, uint64_t end);

uint64_t inner_prod(const sdsl::bit_vector &first, const sdsl::bit_vector &second);

void compute_or(const std::vector<const bit_vector *> &columns,
                sdsl::bit_vector *result,
                ThreadPool &thread_pool);

// Call this version only for sparse vectors (with the density about 1% or less).
// The buffer must have capacity to store 3 x (number of set bits in all columns)
// 64-bit integers.
std::unique_ptr<bit_vector> compute_or(const std::vector<const bit_vector *> &columns,
                                       uint64_t *buffer,
                                       ThreadPool &thread_pool);

// Assumes that all bits that are set in |column| are set in |reference| too
sdsl::bit_vector generate_subindex(const bit_vector &column,
                                   const sdsl::bit_vector &reference,
                                   uint64_t reference_num_set_bits,
                                   ThreadPool &thread_pool);
// Assumes that all bits that are set in |column| are set in |reference| too
sdsl::bit_vector generate_subindex(const bit_vector &column,
                                   const bit_vector &reference,
                                   ThreadPool &thread_pool);

// Apply the bitwise AND of vector with right-shifts of itself. Only works for
// values of offset < 64
sdsl::bit_vector autocorrelate(const sdsl::bit_vector &vector, uint8_t offset);


// Call (uint64_t index, uint64_t value) for each non-zero value in [begin, end)
template <class Callback>
void call_nonzeros(const sdsl::int_vector<> &vector,
                   uint64_t begin, uint64_t end,
                   Callback callback) {
    if (begin >= end)
        return;

    uint64_t i = begin;
    switch (vector.width()) {
        case 64:
            std::for_each(reinterpret_cast<const uint64_t*>(vector.data()) + begin,
                          reinterpret_cast<const uint64_t*>(vector.data()) + end,
                          [&](auto value) { if (value) callback(i, value); ++i; });
            break;

        case 32:
            std::for_each(reinterpret_cast<const uint32_t*>(vector.data()) + begin,
                          reinterpret_cast<const uint32_t*>(vector.data()) + end,
                          [&](auto value) { if (value) callback(i, value); ++i; });
            break;

        case 16:
            std::for_each(reinterpret_cast<const uint16_t*>(vector.data()) + begin,
                          reinterpret_cast<const uint16_t*>(vector.data()) + end,
                          [&](auto value) { if (value) callback(i, value); ++i; });
            break;

        case 8:
            std::for_each(reinterpret_cast<const uint8_t*>(vector.data()) + begin,
                          reinterpret_cast<const uint8_t*>(vector.data()) + end,
                          [&](auto value) { if (value) callback(i, value); ++i; });
            break;

        default:
            // vector.width() is not a power of two
            assert(vector.width() < 64);

            size_t begin_word = begin * vector.width() / 64;
            size_t end_word = (end * vector.width() + 63) / 64;

            for (uint64_t w = begin_word, it = begin; w < end_word; ++w) {
                if (!vector.data()[w])
                    continue;

                it = std::max(it, w * 64 / vector.width());

                auto it_end = std::min(end,
                    ((w + 1) * 64 + vector.width() - 1) / vector.width());

                for (; it < it_end; ++it) {
                    if (vector[it])
                        callback(it, vector[it]);
                }
            }
    }
}

// Call (uint64_t index, uint64_t value) for each non-zero value in |vector|.
template <class Callback>
void call_nonzeros(const sdsl::int_vector<> &vector, Callback callback) {
    return call_nonzeros(vector, 0, vector.size(), callback);
}

template <typename BitVector>
inline uint64_t next1(const BitVector &v,
                      uint64_t pos,
                      size_t num_steps) {
    assert(pos < v.size());

    for (size_t t = 1; t < num_steps; ++t, ++pos) {
        if (pos == v.size() || v[pos])
            return pos;
    }
    if (pos == v.size())
        return pos;

    uint64_t rk;

    if (num_steps >= 1) {
        auto pair = v.inverse_select(pos);
        if (pair.first)
            return pos;

        rk = pair.second + 1;

    } else {
        rk = pos ? v.rank1(pos - 1) + 1 : 1;
    }

    return rk <= v.num_set_bits() ? v.select1(rk) : v.size();
}

template <typename BitVector>
inline uint64_t prev1(const BitVector &v,
                      uint64_t pos,
                      size_t num_steps) {
    assert(pos < v.size());

    for (size_t t = 1; t < num_steps; ++t, --pos) {
        if (v[pos])
            return pos;

        if (pos == 0)
            return v.size();
    }

    uint64_t rk;

    if (num_steps >= 1) {
        auto pair = v.inverse_select(pos);
        if (pair.first)
            return pos;

        rk = pair.second;

    } else {
        rk = v.rank1(pos);
    }

    return rk ? v.select1(rk) : v.size();
}


// taken from https://github.com/xxsds/sdsl-lite/blob/master/include/sdsl/util.hpp
// this function has been modified.
template <class t_int_vec>
typename t_int_vec::size_type
next_bit(const t_int_vec &v,
         uint64_t idx,
         uint64_t max_steps = std::numeric_limits<uint64_t>::max()) {
    assert(idx < v.bit_size());

    uint64_t pos  = idx >> 6;
    uint64_t node = v.data()[pos];
    node >>= (idx & 0x3F);
    if (node)
        return std::min(idx + sdsl::bits::lo(node), v.bit_size());

    uint64_t end = idx + std::min(max_steps, v.bit_size() - idx);
    for (++pos; (pos << 6) < end; ++pos) {
        if (v.data()[pos])
            return std::min((pos << 6) | sdsl::bits::lo(v.data()[pos]), v.bit_size());
    }
    return v.bit_size();
}

// taken from https://github.com/xxsds/sdsl-lite/blob/master/include/sdsl/util.hpp
// this function has been modified.
template <class t_int_vec>
typename t_int_vec::size_type
prev_bit(const t_int_vec &v,
         uint64_t idx,
         uint64_t max_steps = std::numeric_limits<uint64_t>::max()) {
    assert(idx < v.bit_size());

    int64_t pos  = idx >> 6;
    uint64_t node = v.data()[pos];
    node <<= 63 - (idx & 0x3F);
    if (node)
        return idx - (63 - sdsl::bits::hi(node));

    // last position to visit: 0 or (idx + 1 - max_steps)
    int64_t r_end_word = ((idx + 1 - std::min(idx + 1, max_steps)) >> 6) - 1;
    assert(r_end_word >= -1);
    for (--pos; pos > r_end_word; --pos) {
        if (v.data()[pos])
            return (pos << 6) | sdsl::bits::hi(v.data()[pos]);
    }
    return v.bit_size();
}


// Return an int_vector whose limbs are aligned to alignment bytes.
// This is useful for algorithms requiring access to larger words in the underlying
// vector (atomic increment, SIMD, etc.)
template <int width = 0>
inline sdsl::int_vector<width> aligned_int_vector(size_t size = 0, uint64_t val = 0,
                                                  uint8_t var_width = 0,
                                                  size_t alignment = 8) {
    assert(__builtin_popcountll(alignment) == 1);

    // This is a dirty hack to allow for reallocating an int_vector<width>'s
    // underlying storage
    struct int_vector_access {
        typename sdsl::int_vector<width>::size_type m_size;
        uint64_t *m_data;
        typename sdsl::int_vector<width>::int_width_type m_width;
    };
    static_assert(sizeof(sdsl::int_vector<width>) == sizeof(int_vector_access));
    assert(!width || width == var_width);

    sdsl::int_vector<width> v;
    auto &v_cast = reinterpret_cast<int_vector_access&>(v);
    v_cast.m_size = size * var_width;
    v_cast.m_width = var_width;
    free(v_cast.m_data);

    size_t capacity_bytes = (((((v_cast.m_size + 63) >> 6) << 3) + alignment - 1) / alignment) * alignment;
    if (posix_memalign((void**)&v_cast.m_data, alignment, capacity_bytes) || !v_cast.m_data)
        throw std::bad_alloc();

    memset(v_cast.m_data, 0, capacity_bytes);

    if (val)
        sdsl::util::set_to_value(v, val);

    return v;
}

inline sdsl::bit_vector aligned_bit_vector(size_t size = 0, bool val = false,
                                           size_t alignment = 8) {
    return aligned_int_vector<1>(size, val, 1, alignment);
}


namespace sdsl {

// based on sdsl::select_support_scan
template <uint8_t t_b = 1, uint8_t t_pat_len = 1>
class select_support_scan_offset : public select_support_scan<t_b, t_pat_len> {
  public:
    using typename select_support_scan<t_b, t_pat_len>::size_type;

    explicit select_support_scan_offset(const bit_vector *v = nullptr)
          : select_support_scan<t_b, t_pat_len>::select_support_scan(v) {}

    select_support_scan_offset(const select_support_scan<t_b,t_pat_len> &ss)
          : select_support_scan<t_b, t_pat_len>::select_support_scan(ss.m_v) {}

    inline size_type select_offset(size_type i, size_type offset = 0) const {
        using trait = select_support_trait<t_b, t_pat_len>;
        const uint64_t *data = this->m_v->data() + (offset >> 6);
        size_type word_pos = offset >> 6;
        size_type word_off = offset & 0x3F;
        uint64_t carry = trait::init_carry(data, word_pos);
        size_type args = trait::args_in_the_first_word(*data, word_off, carry);
        if (args >= i) {
            return (word_pos << 6)
                + trait::ith_arg_pos_in_the_first_word(*data, i, word_off, carry);
        }
        word_pos++;
        size_type sum_args = args;
        carry = trait::get_carry(*data);
        uint64_t old_carry = carry;
        args = trait::args_in_the_word(*(++data), carry);
        while (sum_args + args < i) {
            sum_args += args;
            assert(data + 1 < this->m_v->data() + (this->m_v->capacity() >> 6));
            old_carry = carry;
            args = trait::args_in_the_word(*(++data), carry);
            word_pos++;
        }
        return (word_pos << 6)
            + trait::ith_arg_pos_in_the_word(*data, i - sum_args, old_carry);
    }
};

} // namespace sdsl

// Predict the memory footprint in bits for sdsl::sd_vector<>
uint64_t footprint_sd_vector(uint64_t size, uint64_t num_set_bits);

// Predict the memory footprint in bits for sdsl::select_support_mcl<>
uint64_t footprint_select_support_mcl(uint64_t size, uint64_t num_set_bits);

// Predict the memory footprint in bits for sdsl::rank_support_v5<>
uint64_t footprint_rank_support_v5(uint64_t size);

#endif // __VECTOR_ALGORITHM_HPP__
