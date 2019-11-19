#ifndef __VECTORS_HPP__
#define __VECTORS_HPP__

#if _USE_FOLLY
#include <folly/FBVector.h>
#include <folly/small_vector.h>
    template <typename... Args>
    using Vector = folly::fbvector<Args...>;

    template <typename T, size_t NumReserved = 2, typename SizeType = uint32_t>
    using SmallVector = folly::small_vector<T, NumReserved, SizeType>;
#else
#include <vector>
#include <cstdint>
    template <typename... Args>
    using Vector = std::vector<Args...>;

    template <typename T, size_t NumReserved = 2, typename SizeType = uint32_t>
    using SmallVector = std::vector<T>;
#endif

#endif // __VECTORS_HPP__
