#ifndef __ANNOBLOOM__
#define __ANNOBLOOM__

#include <iostream>
#include <cassert>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <unordered_set>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include <MurmurHash3.h>
#include <ntHashIterator.hpp>


namespace annotate {

//MERGING

std::vector<uint64_t> merge_or(const std::vector<uint64_t> &a,
                                      const std::vector<uint64_t> &b);
std::vector<uint64_t> merge_and(const std::vector<uint64_t> &a,
                                       const std::vector<uint64_t> &b);

uint64_t popcount(const std::vector<uint64_t> &a);

bool equal(const std::vector<uint64_t> &a, const std::vector<uint64_t> &b);

bool test_bit(const std::vector<uint64_t> &a, const size_t col);

void set_bit(std::vector<uint64_t> &a, const size_t col);

//multihash container
class MultiHash {
    public:
      MultiHash(size_t num_hash = 0) : hashes_(num_hash) { }

      template <typename T>
      MultiHash(const T *data, const size_t num_hash = 0) : hashes_(num_hash) {
          //std::copy(data, data + num_hash, hashes_.begin());
          for (size_t i = 0; i < num_hash; ++i) {
              hashes_[i] = data[i];
          }
      }

      size_t size() const { return hashes_.size(); }

      size_t& operator[](const size_t ind) { return hashes_[ind]; }

      const size_t& operator[](const size_t ind) const { return hashes_[ind]; }

      std::vector<size_t>::iterator begin() { return hashes_.begin(); }

      std::vector<size_t>::iterator end() { return hashes_.end(); }

      std::vector<size_t>::const_iterator begin() const { return hashes_.begin(); }

      std::vector<size_t>::const_iterator end() const { return hashes_.end(); }

    //private:
    public:
      //vector to store hashes (one per seed/hash function)
      std::vector<size_t> hashes_;
};

//HASH WRAPPER STRUCTS
struct {
    template <typename T, typename S>
    void operator()(const T *data, const size_t len, const uint32_t seed, S *hash) {
        //WARNING: make sure that the size of space allocated to hash is at least 8
        uint64_t bigint[2];
        const char *begin = reinterpret_cast<const char*>(data);
        MurmurHash3_x64_128(
                begin,
                len * sizeof(T),
                seed,
                &bigint[0]
        );
        *hash = bigint[0];
    }
} Murmur3Hasher;

//same interface as ntHash, for turning any hash function into an iterative multihasher
class HashIterator {
    public:
      HashIterator& operator++() {
          for (size_t i = 0; i < hashes_.size(); ++i) {
              //use index as seed
              Murmur3Hasher(seq_cur, k_, i, &hashes_[i]);
          }
          seq_cur++;
          return *this;
      }

      HashIterator(const std::string &sequence, const size_t num_hash, const size_t k)
          : seq_begin(sequence.c_str()),
            seq_cur(seq_begin),
            seq_end(sequence.c_str() + sequence.length()),
            hashes_(num_hash),
            k_(k) {
          assert(sequence.length() >= k_);
          operator++();
      }

      bool operator==(const char *that) const {
          return seq_cur == that;
      }

      bool operator!=(const char *that) const {
          return seq_cur != that;
      }

      const char* end() const {
          return seq_end - k_ + 2;
      }

      size_t size() const {
          return hashes_.size();
      }

      const uint64_t* operator*() const {
          return hashes_.data();
      }

      size_t pos() const {
          return seq_cur - seq_begin - 1;
      }

    private:
      const char *seq_begin, *seq_cur, *seq_end;
      std::vector<uint64_t> hashes_;
      size_t k_;
};


template <class HashIt>
std::vector<MultiHash> hash(HashIt&& hash_it, const size_t num_hash);

std::vector<MultiHash> hash_murmur(
        const std::string &sequence,
        const size_t num_hash,
        const size_t k);

std::vector<size_t> annotate(const MultiHash &multihash, const size_t max_size);

class ExactFilter {
  public:
    ExactFilter(size_t num_hash_functions = 0, size_t size = 0) {
        num_hash_functions++;
        size++;
    }

    template <typename Object>
    bool find(Object *a, Object *b) const {
        return set_.find(hash_(a, b)) != set_.end();
    }

    bool find(const MultiHash &hash) const {
        assert(hash.size());
        return set_.find(hash[0]) != set_.end();
    }

    template <typename Object>
    bool insert(Object *a, Object *b) {
        return !set_.insert(hash_(a, b)).second;
    }

    bool insert(const MultiHash &hash) {
        assert(hash.size());
        return !set_.insert(hash[0]).second;
    }

    bool operator==(const ExactFilter &that) const {
        if (set_.size() != that.set_.size())
            return false;
        auto it = set_.begin();
        auto jt = that.set_.begin();
        for (; it != set_.end(); ++it, ++jt) {
            if (*it != *jt)
                return false;
        }
        return true;
    }

    bool operator!=(const ExactFilter &that) const {
        if (set_.size() != that.set_.size())
            return true;
        auto it = set_.begin();
        auto jt = that.set_.begin();
        for (; it != set_.end(); ++it, ++jt) {
            if (*it == *jt)
                return false;
        }
        return true;
    }

  private:
    std::unordered_set<uint64_t> set_;

    template <typename Object>
    static uint64_t hash_(Object *a, Object *b) {
        __uint128_t hash = 0;
        Murmur3Hasher(a, b - a, 0, &hash);
        uint64_t first_hash;
        std::copy(reinterpret_cast<const uint64_t*>(&hash),
                  reinterpret_cast<const uint64_t*>(&hash) + 1, &first_hash);
        return first_hash;
    }
};


class BloomFilter {
  public:
    BloomFilter(size_t num_hash_functions, size_t _n_bits = 0) {
        seeds.resize(num_hash_functions);
        std::iota(seeds.begin(), seeds.end(), 0);

        n_bits = _n_bits;
        if (n_bits > 0)
            bits.resize((n_bits >> 6) + 1);
    }

    size_t size() {
        return n_bits;
    }

    void resize(size_t new_size) {
        n_bits = new_size;
        bits.resize((n_bits >> 6) + 1);
    }

  public:
    template <typename T>
    bool find(T *a, T *b) const {
        return hash_helper(a, b);
    }

    template <typename T>
    bool find(const std::vector<T> &a) const {
        return find(&(*(a.begin())), &(*(a.end())));
    }
    template <typename T>
    bool find(const T &a) const {
        return find(&a, &a + sizeof(a));
    }

    bool find(const MultiHash &multihash) const {
        for (auto it = multihash.begin(); it != multihash.end(); ++it) {
            if (!test_bit(bits, *it % n_bits)) {
                return false;
            }
        }
        return true;
        //return equal(merge_or(bits, annotate(hash, n_bits)), bits);
    }

    bool insert(const MultiHash &multihash) {
        bool might_contain = true;
        for (auto it = multihash.begin(); it != multihash.end(); ++it) {
            if (might_contain && !test_bit(bits, *it % n_bits)) {
                might_contain = false;
            }
            set_bit(bits, *it % n_bits);
        }
        return might_contain;
        /*
        auto merged = merge_or(bits, annotate(hash, n_bits));
        bool might_contain = equal(merged, bits);
        std::copy(merged.begin(), merged.end(), bits.begin());
        return might_contain;
        */
    }

    template <typename T>
    bool insert(T *a, T *b) {
        return hash_helper_insert(a, b);
    }

    void serialize(std::ostream &out) const {
        out << n_bits << "\n";
        boost::archive::binary_oarchive oar(out);
        oar & seeds;
        oar & bits;
    }

    void deserialize(std::istream &in) {
        in >> n_bits;
        boost::archive::binary_iarchive iar(in);
        iar & seeds;
        iar & bits;
    }

    bool operator==(const BloomFilter &a) const {
        if (n_bits != a.n_bits) {
            std::cerr << "Different number of bits\n";
            std::cerr << n_bits << " " << a.n_bits << "\n";
            return false;
        }
        if (seeds.size() != a.seeds.size()) {
            std::cerr << "Different number of seeds\n";
            std::cerr << seeds.size() << " " << a.seeds.size() << "\n";
            return false;
        }
        auto jt = a.seeds.begin();
        for (auto it = seeds.begin(); it != seeds.end(); ++it, ++jt) {
            if (*it != *jt) {
                std::cerr << "Seeds not equal\n";
                return false;
            }
        }
        if (bits.size() != a.bits.size()) {
            std::cerr << "Different number of blocks\n";
            std::cerr << bits.size() << " " << a.bits.size() << "\n";
            return false;
        }
        jt = a.bits.begin();
        for (auto it = bits.begin(); it != bits.end(); ++it, ++jt) {
            if (*it != *jt) {
                std::cerr << "Different bits\n";
                std::cerr << *it << " " << *jt << "\n";
                return false;
            }
        }
        return true;
    }

    bool operator!=(const BloomFilter &a) const {
        return !(*this == a);
    }

    double occupancy() const {
        size_t count = 0;
        for (auto it = bits.begin(); it != bits.end(); ++it) {
            count += __builtin_popcountl(*it);
        }
        return static_cast<double>(count) / (bits.size() * 64);
    }

  private:
    template <typename T>
    inline bool hash_helper(T *a, T *b) const {
        __uint128_t hash = 0;
        uint64_t select;

        for (auto it = seeds.begin(); it != seeds.end(); ++it) {
            Murmur3Hasher(a, b - a, *it, &hash);
            select = hash % n_bits;
            /*
            auto jt = bits.begin() + (select >> 6);
            if (jt >= bits.end()) {
                std::cerr << "Out of bounds\n";
                exit(1);
            }
            if ((*jt | (1llu << (select % 64))) != *jt) {
                return false;
            }
            */
            if (!test_bit(bits, select)) {
                return false;
            }
        }
        return true;
    }

    template <typename T>
    inline bool hash_helper_insert(T *a, T *b) {
        assert(b >= a);

        __uint128_t hash = 0;
        uint64_t select;
        bool might_contain = true;
        //uint64_t last;
        for (auto it = seeds.begin(); it != seeds.end(); ++it) {
            Murmur3Hasher(a, b - a, *it, &hash);
            select = hash % n_bits;

            /*
            auto jt = bits.begin() + (select >> 6);
            assert(jt < bits.end());

            last = *jt;
            *jt |= 1llu << (select % 64);
            if (might_contain && *jt != last) {
                might_contain = false;
            }
            */
            if (might_contain && !test_bit(bits, select)) {
                might_contain = false;
            }
            set_bit(bits, select);
        }
        return might_contain;
    }

    std::vector<uint64_t> bits;
    uint64_t n_bits = 0;
    std::vector<uint64_t> seeds;
};


template <class Filter = ExactFilter>
class HashAnnotation {
  public:
    HashAnnotation(size_t num_hash_functions = 0)
          : num_hash_functions_(num_hash_functions) {}

    void resize(size_t size) {
        assert(size > color_bits.size());
        color_bits.resize(size, Filter(num_hash_functions_));
    }

    Filter& operator[](size_t i) {
        return color_bits[i];
    }

    const Filter& operator[](size_t i) const {
        return color_bits[i];
    }

    template <typename T, typename S>
    std::vector<size_t> insert(T *a, T *b, S begin, S end) {
        std::vector<size_t> annot((color_bits.size() >> 6) + 1);
        for (auto it = begin; it != end; ++it) {
            if (*it >= color_bits.size()) {
                std::cerr << "ERROR: Index " << *it << " >= " << color_bits.size() << "\n";
                exit(1);
            }
            if (color_bits[*it].insert(a, b)) {
                //might contain
                set_bit(annot, *it);
                //annot[*it >> 6] |= (1llu << (*it % 64));
            }
        }
        return annot;
    }

    template <typename T>
    std::vector<size_t> insert(T *a, T *b, size_t ind) {
        return insert(a, b, &ind, &ind + 1);
    }

    template <typename S>
    std::vector<size_t> insert(const MultiHash &hash, S begin, S end) {
        std::vector<size_t> annot((color_bits.size() >> 6) + 1);
        for (auto it = begin; it != end; ++it) {
            if (*it >= color_bits.size()) {
                std::cerr << "ERROR: Index " << *it << " >= " << color_bits.size() << "\n";
                exit(1);
            }
            if (color_bits[*it].insert(hash)) {
                set_bit(annot, *it);
                //annot[*it >> 6] |= (1llu << (*it % 64));
            }
        }
        return annot;
    }

    std::vector<size_t> insert(const MultiHash &hash, size_t ind) {
        return insert(hash, &ind, &ind + 1);
    }

    template <typename T, typename S>
    std::vector<size_t> find(T *a, T *b, S begin, S end) const {
        std::vector<size_t> annot((color_bits.size() >> 6) + 1, 0);
        for (auto it = begin; it != end; ++it) {
            if (*it >= color_bits.size()) {
                std::cerr << "ERROR: Index " << *it << " >= " << color_bits.size() << "\n";
                exit(1);
            }
            if (color_bits[*it].find(a, b)) {
                //might contain
                set_bit(annot, *it);
                //annot[*it >> 6] |= (1llu << (*it % 64));
            }
        }
        return annot;
    }

    template <typename T>
    std::vector<size_t> find(T *a, T *b) const {
        std::vector<size_t> annot((color_bits.size() >> 6) + 1);
        for (size_t i = 0; i < color_bits.size(); ++i) {
            if (color_bits[i].find(a, b)) {
                set_bit(annot, i);
                //annot[i >> 6] |= (1llu << (i % 64));
            }
        }
        return annot;
    }

    template <typename T>
    std::vector<size_t> find(T *a, T *b, size_t ind) const {
        return find(a, b, &ind, &ind + 1);
    }

    template <typename S>
    std::vector<size_t> find(const MultiHash &hash, S begin, S end) const {
        std::vector<size_t> annot((color_bits.size() >> 6) + 1);
        for (auto it = begin; it != end; ++it) {
            if (*it >= color_bits.size()) {
                std::cerr << "ERROR: Index " << *it << " >= " << color_bits.size() << "\n";
                exit(1);
            }
            if (color_bits[*it].find(hash)) {
                set_bit(annot, *it);
                //annot[*it >> 6] |= (1llu << (*it % 64));
            }
        }
        return annot;
    }

    std::vector<size_t> find(const MultiHash &hash, size_t ind) const {
        return find(hash, &ind, &ind + 1);
    }

    std::vector<size_t> find(const MultiHash &hash) const {
        std::vector<size_t> annot((color_bits.size() >> 6) + 1);
        for (size_t i = 0; i < color_bits.size(); ++i) {
            if (color_bits[i].find(hash)) {
                set_bit(annot, i);
            }
        }
        return annot;
    }


    size_t size() const {
        return color_bits.size();
    }

    size_t num_hash_functions() const {
        return num_hash_functions_;
    }

    void serialize(std::ostream &out) const {
        size_t size = color_bits.size();
        out.write(reinterpret_cast<const char*>(&size), sizeof(size));
        //out << color_bits.size() << "\n";
        for (auto it = color_bits.begin(); it != color_bits.end(); ++it) {
            it->serialize(out);
        }
    }

    void deserialize(std::istream &in) {
        size_t size;
        //in >> size;
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        color_bits.resize(size);
        for (auto it = color_bits.begin(); it != color_bits.end(); ++it) {
            it->deserialize(in);
        }
    }

    bool operator==(const HashAnnotation<Filter> &a) const {
        if (color_bits.size() != a.color_bits.size()) {
            std::cerr << "Different number of filters\n";
            return false;
        }
        auto jt = a.color_bits.begin();
        for (auto it = color_bits.begin(); it != color_bits.end(); ++it, ++jt) {
            if (*it != *jt) {
                std::cerr << "Color filter " << (it - color_bits.begin()) << " different\n";
                return false;
            }
        }
        return true;
    }

    bool operator!=(const HashAnnotation<Filter> &a) const {
        return !(*this == a);
    }

    void append_bit(size_t num_bits = 0) {
        color_bits.push_back(Filter(num_hash_functions_, num_bits));
    }

    std::vector<double> occupancy() const {
        std::vector<double> occ(color_bits.size());
        for (size_t i = 0; i < color_bits.size(); ++i) {
            occ[i] = color_bits[i].occupancy();
        }
        return occ;
    }



  private:
    std::vector<Filter> color_bits;
    size_t num_hash_functions_;
};

} // namespace annotate

#endif // __ANNOBLOOM__
