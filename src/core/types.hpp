#pragma once

#include <cstddef>
#include <cstdint>

namespace hypo {

using Real  = double;
using Index = std::size_t;
using Int   = int;

constexpr Int kDefaultHaloWidth = 1;
constexpr Real kDefaultTolerance = 1e-6;
constexpr Index kDefaultMaxIter = 10000;

} // namespace hypo
