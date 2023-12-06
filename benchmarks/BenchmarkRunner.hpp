#pragma once

#include <catch2/catch_approx.hpp>

#include "timer.hpp"

template <typename LogMem, typename LogTime, typename CheckFailed>
class BenchmarkRunner {
public:
    static constexpr double FLOAT_EQ_EPS = 1e-4;

    BenchmarkRunner(MemoryUsage& memory_usage, Timer& timer, LogMem log_mem, LogTime log_time, CheckFailed check_failed)
        : _memory_usage(memory_usage),
          _timer(timer),
          _log_mem(log_mem),
          _log_time(log_time),
          _check_failed(check_failed) {}

    template <typename TSKitFunc, typename SFKitDAGFunc, typename SFKitBPFunc>
    void bench_and_check(
        std::string stat_name, TSKitFunc tskit_func, SFKitDAGFunc sfkit_dag_func, SFKitBPFunc sfkit_bp_func
    ) {
        this->start();
        double const tskit_value = tskit_func();
        do_not_optimize(tskit_value);
        this->stop(stat_name, "tskit");

        this->start();
        double const sfkit_dag_value = sfkit_dag_func();
        do_not_optimize(sfkit_dag_value);
        this->stop(stat_name, "sfkit_dag");

        this->start();
        double const sfkit_bp_value = sfkit_bp_func();
        do_not_optimize(sfkit_bp_value);
        this->stop(stat_name, "sfkit_bp");

        // Does our value match the one computed by tskit?
        if (sfkit_dag_value != Catch::Approx(tskit_value).epsilon(FLOAT_EQ_EPS)) {
            std::cerr << "ERROR !! " << stat_name << " mismatch between tskit and sfkit (DAG)" << std::endl;
            std::cerr << "    " << tskit_value << " vs. " << sfkit_dag_value << " (tskit vs. sfkit)" << std::endl;
            _check_failed();
        }

        if (sfkit_bp_value != Catch::Approx(tskit_value).epsilon(FLOAT_EQ_EPS)) {
            std::cerr << "ERROR !! " << stat_name << " mismatch between tskit and sfkit (BP)" << std::endl;
            std::cerr << "    " << tskit_value << " vs. " << sfkit_bp_value << " (tskit vs. sfkit)" << std::endl;
            _check_failed();
        }
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
    CheckFailed  _check_failed;
};
