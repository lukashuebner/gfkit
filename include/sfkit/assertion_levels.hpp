#pragma once

/// @brief Assertion levels
namespace sfkit::assert {
/// @defgroup assertion-levels Assertion levels
///
/// @{

/// @brief Assertion level for lightweight assertions.
#define SFKIT_ASSERTION_LEVEL_LIGHT 10

/// @brief Assertion level for lightweight assertions.
constexpr int light = SFKIT_ASSERTION_LEVEL_LIGHT;

/// @brief Default assertion level. This level is used if no assertion level is specified.
#define SFKIT_ASSERTION_LEVEL_NORMAL 30

/// @brief Default assertion level. This level is used if no assertion level is specified.
constexpr int normal = SFKIT_ASSERTION_LEVEL_NORMAL;

/// @brief Assertion level for heavyweight assertions.
#define SFKIT_ASSERTION_LEVEL_HEAVY 60

/// @brief Assertion level for heavyweight assertions.
constexpr int heavy = SFKIT_ASSERTION_LEVEL_HEAVY;

/// @}
} // namespace sfkit::assert
