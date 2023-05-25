#pragma once

/// @brief Assertion levels
namespace tdt::assert {
/// @defgroup assertion-levels Assertion levels
///
/// @{

/// @brief Assertion level for lightweight assertions.
#define TDT_ASSERTION_LEVEL_LIGHT 10

/// @brief Assertion level for lightweight assertions.
constexpr int light = TDT_ASSERTION_LEVEL_LIGHT;

/// @brief Default assertion level. This level is used if no assertion level is specified.
#define TDT_ASSERTION_LEVEL_NORMAL 30

/// @brief Default assertion level. This level is used if no assertion level is specified.
constexpr int normal = TDT_ASSERTION_LEVEL_NORMAL;

/// @brief Assertion level for heavyweight assertions.
#define TDT_ASSERTION_LEVEL_HEAVY 60

/// @brief Assertion level for heavyweight assertions.
constexpr int heavy = TDT_ASSERTION_LEVEL_HEAVY;

/// @}
} // namespace tdt::assert
