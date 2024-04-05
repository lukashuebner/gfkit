#pragma once

#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <optional>
#include <tuple>

#ifdef ENABLE_MALLOC_COUNT
    #include <malloc_count.h>
    #include <stack_count.h>
#endif
#include <sys/sysinfo.h>
#include <sys/types.h>

// This assembly magic is from Google Benchmark
// see https://github.com/google/benchmark/blob/v1.7.1/include/benchmark/benchmark.h#L443-L446
template <typename T>
void do_not_optimize(T const& val) {
    asm volatile("" : : "r,m"(val) : "memory");
}

// TODO Why doesn't this work?
// template <typename T>
// inline void do_not_optimize(T& val) {
//     asm volatile("" : "+r,m"(val) : : "memory");
// }

// template <typename T>
// inline void do_not_optimize(T&& val) {
//     asm volatile("" : "+r,m"(val) : : "memory");
// }

// TODO Do we really need special handling for clang and GCC?
// #if defined(__clang__)
//     asm volatile("" : "+r,m"(val) : : "memory");
// #else
//     asm volatile("" : "+m,r"(val) : : "memory");
// #endif

[[nodiscard]] inline bool approx_eq(double a, double b, double eps = 1e-6) {
    return std::abs(a - b) < eps;
}

// TODO Enable the use of timer labels and timers which can be paused and restarted
class Timer {
public:
    using duration = decltype(std::chrono::steady_clock::now() - std::chrono::steady_clock::now());

    Timer() {
        start();
    }

    duration stop() {
        asm("" ::: "memory");
        auto end = std::chrono::steady_clock::now();
        asm("" ::: "memory");

        return end - _start;
    }

    void start() {
        asm("" ::: "memory"); // prevent compiler reordering
        _start = std::chrono::steady_clock::now();
        asm("" ::: "memory");
    }

    template <class Func>
    static duration time_func(Func func) {
        Timer timer;
        func();
        return timer.stop();
    }

private:
    std::chrono::time_point<std::chrono::steady_clock> _start;
};

class MemoryUsage {
public:
    struct Report {
        int64_t               virtmem_delta_byte;
        int64_t               rss_delta_byte;
        int64_t               heap_delta_byte;
        size_t                heap_peak_byte;
        size_t                heap_current_byte;
        size_t                stack_peak_byte;
        std::optional<size_t> data_structrue_byte = std::nullopt;

        Report& data_structure_bytes(size_t size_bytes) {
            data_structrue_byte = size_bytes;
            return *this;
        }
    };

    MemoryUsage() {
        start();
    }

    void start() {
#ifdef ENABLE_MALLOC_COUNT
        asm("" ::: "memory"); // prevent compiler reordering
        _stack_count_base   = stack_count_clear();
        _malloc_count_start = malloc_count_current();
        malloc_count_reset_peak();
#endif
        std::tie(_virtual_base_kib, _rss_base_kib) = virtual_and_rss_kib();
    }

    Report stop(std::optional<size_t> data_structure_bytes = std::nullopt) {
        Report report;

        asm("" ::: "memory");
#ifdef ENABLE_MALLOC_COUNT
        report.stack_peak_byte   = stack_count_usage(_stack_count_base);
        report.heap_peak_byte    = malloc_count_peak();
        report.heap_delta_byte   = malloc_count_current() - _malloc_count_start;
        report.heap_current_byte = malloc_count_current();
#endif
        auto [virtual_kib, rss_kib] = virtual_and_rss_kib();
        report.virtmem_delta_byte   = kib_to_byte(_virtual_base_kib - virtual_kib);
        report.rss_delta_byte       = kib_to_byte(_rss_base_kib - rss_kib);
        asm("" ::: "memory");
        return report;
    }

private:
    std::tuple<int64_t, int64_t> virtual_and_rss_kib() const {
        std::optional<size_t> virtual_kib = std::nullopt;
        std::optional<size_t> rss_kib     = std::nullopt;

        // Virtual memory and Resident Set Size (RSS)
        std::ifstream file("/proc/self/status");
        if (file.is_open()) {
            std::string line;

            while (std::getline(file, line)) {
                std::istringstream iss(line);
                if (!(line.starts_with("VmSize:") || line.starts_with("VmRSS:"))) {
                    continue;
                }
                std::string key, value, unit;
                iss >> key >> value >> unit;

                if (unit != "kB") {
                    std::cerr << "WARNING ! Unit of VmSize: in /proc/self/status is not 'kB'" << std::endl;
                }

                try {
                    if (key == "VmSize:") {
                        virtual_kib = std::stoul(value);
                    } else if (key == "VmRSS:") {
                        rss_kib = std::stoul(value);
                    }
                } catch (std::exception& e) {
                    std::cerr << "WARNING ! Failed to parse VmSize or VmRSS from /proc/self/status (" << key << ": "
                              << value << ")." << std::endl;
                }
            }
        } else {
            std::cerr << "WARNING ! Failed to open /proc/self/status" << std::endl;
        }

        if (!virtual_kib || !rss_kib) {
            std::cerr << "WARNING ! Failed to read VmSize or VmRSS from /proc/self/status" << std::endl;
        }

        return std::make_tuple(*virtual_kib, *rss_kib);
    }

    int64_t round_to_kib(int64_t bytes) {
        return (bytes + 1023) / 1024;
    }

    int64_t kib_to_byte(int64_t kib) {
        return kib * 1024;
    }

#ifdef ENABLE_MALLOC_COUNT
    void*   _stack_count_base   = nullptr;
    int64_t _malloc_count_start = 0;
#endif
    int64_t _virtual_base_kib = 0;
    int64_t _rss_base_kib     = 0;
};
