#include "hashers.hpp"

#include "cyclichash.h"


namespace annotate {

std::vector<uint64_t> merge_or(const std::vector<uint64_t> &a,
                               const std::vector<uint64_t> &b) {
    assert(a.size() == b.size() && "ORing different sizes");

    std::vector<uint64_t> merged(a.size());
    for (size_t i = 0; i < merged.size(); ++i) {
        merged[i] = a[i] | b[i];
    }
    return merged;
}

std::vector<uint64_t> merge_and(const std::vector<uint64_t> &a,
                                const std::vector<uint64_t> &b) {

    assert(a.size() == b.size() && "ANDing different sizes");

    std::vector<uint64_t> merged(a.size());
    for (size_t i = 0; i < merged.size(); ++i) {
        merged[i] = a[i] & b[i];
    }
    return merged;
}

uint64_t popcount(const std::vector<uint64_t> &a) {
    uint64_t popcount = 0;
    for (auto value : a) {
        popcount += __builtin_popcountl(value);
    }
    return popcount;
}

bool equal(const std::vector<uint64_t> &a,
           const std::vector<uint64_t> &b) {
    assert(a.size() == b.size() && "Checking different sizes");
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] != b[i])
            return false;
    }
    return true;
}

bool test_bit(const std::vector<uint64_t> &a, size_t i) {
    return a[i >> 6] & (1llu << (i % 64));
}

void set_bit(std::vector<uint64_t> &a, size_t i) {
    a[i >> 6] |= 1llu << (i % 64);
}

void print(const std::vector<uint64_t> &a) {
    for (auto it = a.begin(); it != a.end(); ++it) {
        std::cout << *it << " ";
    }
    std::cout << "\n";
}


//CyclicHash
CyclicMultiHash::CyclicMultiHash(const char *data, size_t k, size_t num_hash)
      : hashes_(num_hash),
        k_(k),
        cache_(data, k),
        begin_(0),
        chashers_(num_hash, NULL) {
    assert(k_);

    for (uint32_t j = 0; j < hashes_.size(); ++j) {
        auto *cyclic_hash = new CyclicHash<uint64_t>(k_, j, j + 1, 64lu);
        for (size_t i = 0; i < k_; ++i) {
            cyclic_hash->eat(data[i]);
        }
        chashers_[j] = cyclic_hash;
        hashes_[j] = cyclic_hash->hashvalue;
    }
}

CyclicMultiHash::~CyclicMultiHash() {
    for (size_t i = 0; i < chashers_.size(); ++i) {
        delete reinterpret_cast<CyclicHash<uint64_t>*>(chashers_[i]);
    }
}

bool CyclicMultiHash::reinitialize(const char *data, size_t k, size_t num_hash) {
    if (k != k_ || num_hash != hashes_.size())
        return false;

    for (size_t i = 0; i < k_; ++i) {
        update(data[i]);
    }
    return true;
}

void CyclicMultiHash::update(char next) {
    for (size_t i = 0; i < chashers_.size(); ++i) {
        auto *cyclic_hash = reinterpret_cast<CyclicHash<uint64_t>*>(chashers_[i]);
        cyclic_hash->update(cache_[begin_], next);
        hashes_[i] = cyclic_hash->hashvalue;
    }
    cache_[begin_] = next;
    begin_ = (begin_ == cache_.size() - 1 ? 0 : begin_ + 1);
}

void CyclicMultiHash::reverse_update(char prev) {
    begin_ = (begin_ == 0 ? cache_.size() - 1 : begin_ - 1);
    for (size_t i = 0; i < chashers_.size(); ++i) {
        auto *cyclic_hash = reinterpret_cast<CyclicHash<uint64_t>*>(chashers_[i]);
        cyclic_hash->reverse_update(prev, cache_[begin_]);
        hashes_[i] = cyclic_hash->hashvalue;
    }
    cache_[begin_] = prev;
}


//CyclicHashIterator
CyclicHashIterator::CyclicHashIterator(const char *begin, const char *end,
                                       size_t k, size_t num_hash)
      : hasher_(begin, k, num_hash),
        next_(begin + k),
        end_(end) {
    assert(next_ < end_ || (next_ == end_ && !is_end()));
}

CyclicHashIterator::CyclicHashIterator(const std::string &sequence,
                                       size_t k, size_t num_hash)
      : CyclicHashIterator(&sequence.front(), &sequence.back() + 1, k, num_hash) {}

CyclicHashIterator& CyclicHashIterator::operator++() {
    if (next_ < end_) {
        hasher_.update(*next_);
    }
    next_++;
    return *this;
}


} // namespace annotate
