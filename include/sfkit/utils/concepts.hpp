#pragma once

#include <type_traits>

template <typename T>
concept TriviallyCopyable = std::is_trivially_copyable_v<T>;

template <typename T>
concept PlainStorageContainer = requires(T t) {
    { t.data() } -> std::convertible_to<void const*>;
    { t.size() } -> std::convertible_to<std::size_t>;
};

template <typename T>
concept IterableInput = requires(T const& t) {
    { t.begin() } -> std::input_iterator;
    { t.end() } -> std::input_iterator;
};
