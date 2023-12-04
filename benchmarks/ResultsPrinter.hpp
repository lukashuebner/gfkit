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

#include "timer.hpp"

// TODO Split this up in .hpp and .cpp file
class ResultsPrinter {
public:
    ResultsPrinter(std::ostream& stream, std::string const& revision, std::string const& machine_id)
        : _stream(stream),
          _revision(revision),
          _machine_id(machine_id) {}

    void print_header() {
        _stream << "section,variant,dataset,revision,machine_id,iteration,variable,value,unit\n";
    }

    template <typename Value>
    void print(
        bool const         warmup,
        std::string const& section,
        std::string const& variant,
        std::string const& dataset,
        std::string const& variable,
        Value const&       value,
        std::string const& unit,
        size_t             iteration
    ) {
        if (!warmup) {
            _stream << section << ",";
            _stream << variant << ",";
            _stream << dataset << ",";
            _stream << _revision << ",";
            _stream << _machine_id << ",";
            _stream << iteration << ",";
            _stream << variable << ",";
            _stream << value << ",";
            _stream << unit << "\n";
        }
    }

    // TODO Implement print_time with typechecking (that the given time is in nanoseconds)
    // alternatively work directly with the chrono::duration type
    void print_timer(
        bool const             warmup,
        std::string const&     section,
        std::string const&     variant,
        std::string const&     dataset,
        Timer::duration const& duration,
        size_t                 iteration
    ) {
        auto duration_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count();
        print(warmup, section, variant, dataset, "walltime", duration_ns, "ns", iteration);
    }

    void print_memory(
        bool const                 warmup,
        std::string const&         section,
        std::string const&         variant,
        std::string const&         dataset,
        MemoryUsage::Report const& report,
        size_t                     iteration
    ) {
        auto _print = [this, warmup, section, variant, dataset, iteration](std::string const& variable, auto value) {
            print(warmup, section, variant, dataset, variable, value, "byte", iteration);
        };

        _print("virtmem_delta", report.virtmem_delta_byte);
        _print("rss_delta", report.rss_delta_byte);
#ifdef ENABLE_MALLOC_COUNT
        _print("stack_peak", report.heap_current_byte);
        _print("heap_peak", report.heap_peak_byte);
        _print("heap_delta", report.heap_delta_byte);
        _print("stack_peak", report.stack_peak_byte);
#endif
        if (report.data_structrue_byte.has_value()) {
            _print("data_structure", report.data_structrue_byte.value());
        }
    }

private:
    std::ostream& _stream;
    std::string   _revision;
    std::string   _machine_id;
};
