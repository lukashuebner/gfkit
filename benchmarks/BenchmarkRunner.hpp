#pragma once

#include <catch2/catch_approx.hpp>

#include "timer.hpp"

template <typename LogMem, typename LogTime>
class BenchmarkRunner {
public:
    static constexpr double FLOAT_EQ_EPS = 1e-4;

    BenchmarkRunner(MemoryUsage& memory_usage, Timer& timer, LogMem log_mem, LogTime log_time)
        : _memory_usage(memory_usage),
          _timer(timer),
          _log_mem(log_mem),
          _log_time(log_time) {}

    template <typename Func>
    auto bench(std::string stat_name, Func func, std::string const& variant_name) {
        this->start();
        double const value = func();
        do_not_optimize(value);
        this->stop(stat_name, variant_name);

        return value;
    }

    void start() {
        _memory_usage.start();
        _timer.start();
    }

    void stop(std::string const& section, std::string const& variant) {
        _log_time(section, variant, _timer.stop());
        _log_mem(section, variant, _memory_usage.stop());
    }

private:
    Timer&       _timer;
    MemoryUsage& _memory_usage;
    LogTime      _log_time;
    LogMem       _log_mem;
};
