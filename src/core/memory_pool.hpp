#pragma once

#include "core/aligned_buffer.hpp"
#include <vector>
#include <cstddef>

namespace hypo {

/**
 * @brief Simple bump allocator for temporary grid-sized allocations.
 * All memory is freed on destruction. No individual deallocation.
 */
class MemoryPool {
public:
    explicit MemoryPool(std::size_t initialBytes = 64 * 1024 * 1024);
    ~MemoryPool();

    MemoryPool(const MemoryPool&) = delete;
    MemoryPool& operator=(const MemoryPool&) = delete;
    MemoryPool(MemoryPool&&) = delete;
    MemoryPool& operator=(MemoryPool&&) = delete;

    /**
     * @brief Allocate 'bytes' from the pool. Aligned to 64 bytes.
     * @return Pointer to allocated memory, or nullptr if pool exhausted.
     */
    void* allocate(std::size_t bytes);

    /**
     * @brief Reset the pool to reuse all memory.
     */
    void reset() noexcept;

    std::size_t usedBytes() const noexcept { return used_; }
    std::size_t totalBytes() const noexcept { return total_; }

private:
    AlignedBuffer<std::byte> buffer_;
    std::size_t used_ = 0;
    std::size_t total_ = 0;
};

} // namespace hypo
