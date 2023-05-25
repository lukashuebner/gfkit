#pragma once

#include <type_traits>

// TODO Change documentation
// Declaration of the concept "Hashable", which is satisfied by any type 'T'
// such that for values 'a' of type 'T', the expression std::hash<T>{}(a)
// compiles and its result is convertible to std::size_t
template <typename T>
concept TriviallyCopyable = std::is_trivially_copyable_v<T>;

// TODO Add expected return type
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
