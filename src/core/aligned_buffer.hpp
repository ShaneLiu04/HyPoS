#pragma once

#include "core/types.hpp"
#include <memory>
#include <new>
#include <utility>
#include <cstdlib>

namespace hypo {

/**
 * @brief Aligned buffer with RAII. Uses 64-byte alignment for cache-line
 * optimization and to prevent false sharing.
 */
template <typename T>
class AlignedBuffer {
public:
    AlignedBuffer() noexcept = default;

    explicit AlignedBuffer(Index count) { allocate(count); }

    ~AlignedBuffer() { deallocate(); }

    // Non-copyable
    AlignedBuffer(const AlignedBuffer&) = delete;
    AlignedBuffer& operator=(const AlignedBuffer&) = delete;

    // Movable
    AlignedBuffer(AlignedBuffer&& other) noexcept
        : data_(other.data_)
        , size_(other.size_)
        , capacity_(other.capacity_)
    {
        other.data_ = nullptr;
        other.size_ = 0;
        other.capacity_ = 0;
    }

    AlignedBuffer& operator=(AlignedBuffer&& other) noexcept {
        if (this != &other) {
            deallocate();
            data_ = other.data_;
            size_ = other.size_;
            capacity_ = other.capacity_;
            other.data_ = nullptr;
            other.size_ = 0;
            other.capacity_ = 0;
        }
        return *this;
    }

    void allocate(Index count) {
        deallocate();
        if (count > 0) {
            void* ptr = nullptr;
            constexpr std::size_t kAlignment = 64;
            if (::posix_memalign(&ptr, kAlignment, count * sizeof(T)) != 0) {
                throw std::bad_alloc();
            }
            data_ = static_cast<T*>(ptr);
            size_ = count;
            capacity_ = count;
        }
    }

    void deallocate() noexcept {
        if (data_) {
            ::free(data_); // posix_memalign uses free
            data_ = nullptr;
        }
        size_ = 0;
        capacity_ = 0;
    }

    void resize(Index count) {
        if (count <= capacity_) {
            size_ = count;
            return;
        }
        AlignedBuffer tmp(count);
        // Move-construct elements if needed; for POD types just copy
        for (Index i = 0; i < size_ && i < count; ++i) {
            tmp.data_[i] = data_[i];
        }
        *this = std::move(tmp);
    }

    T* data() noexcept { return data_; }
    const T* data() const noexcept { return data_; }

    Index size() const noexcept { return size_; }
    Index capacity() const noexcept { return capacity_; }

    T& operator[](Index i) noexcept { return data_[i]; }
    const T& operator[](Index i) const noexcept { return data_[i]; }

    T* begin() noexcept { return data_; }
    T* end() noexcept { return data_ + size_; }
    const T* begin() const noexcept { return data_; }
    const T* end() const noexcept { return data_ + size_; }

    bool empty() const noexcept { return size_ == 0; }

    void fill(T value) noexcept {
        for (Index i = 0; i < size_; ++i) {
            data_[i] = value;
        }
    }

private:
    T* data_ = nullptr;
    Index size_ = 0;
    Index capacity_ = 0;
};

} // namespace hypo
