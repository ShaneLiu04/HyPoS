#include "core/memory_pool.hpp"
#include <cstring>

namespace hypo {

MemoryPool::MemoryPool(std::size_t initialBytes)
    : buffer_(initialBytes)
    , total_(initialBytes)
{
}

MemoryPool::~MemoryPool() = default;

void* MemoryPool::allocate(std::size_t bytes) {
    if (bytes == 0) return nullptr;
    // Align to 64 bytes
    std::size_t alignedBytes = (bytes + 63) & ~std::size_t(63);
    if (used_ + alignedBytes > total_) {
        return nullptr; // Pool exhausted
    }
    void* ptr = buffer_.data() + used_;
    used_ += alignedBytes;
    return ptr;
}

void MemoryPool::reset() noexcept {
    used_ = 0;
}

} // namespace hypo
