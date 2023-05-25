// Description: A simple wrapper around perf to profile a function
// Adapted from: https://muehe.org/posts/profiling-only-parts-of-your-code-with-perf/
#pragma once

#include <cstring>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>

#include <fcntl.h>
#include <kassert/kassert.hpp>
#include <signal.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

inline void profile(std::string name, std::function<void()> body) {
    char const* ctl_fifo_name = ".perf.ctl.fifo";
    char const* ack_fifo_name = ".perf.ack.fifo";
    auto        ctl_fifo      = mkfifo(ctl_fifo_name, 0600);
    auto        ack_fifo      = mkfifo(ack_fifo_name, 0600);
    KASSERT(ctl_fifo >= 0, "Failed to create fifo: " << strerror(errno));
    KASSERT(ack_fifo >= 0, "Failed to create fifo: " << strerror(errno));
    std::stringstream control_ss;
    control_ss << "fifo:" << ctl_fifo_name << "," << ack_fifo_name;
    std::cout << "Control: " << control_ss.str() << std::endl;

    std::string perf_out_file = name.find(".data") == std::string::npos ? (name + ".data") : name;

    // Fork off a new process to run perf on this process
    pid_t             pid;
    std::stringstream parent_pid;
    parent_pid << getpid();
    pid = fork();
    if (pid == 0) { // If we're the new process
        // Redirect output to our stderr and stdout to /dev/null
        // auto fd = open("/dev/null", O_RDWR);
        // dup2(fd, 1);
        // dup2(fd, 2);

        // Run perf on the parent process and exit
        exit(execl(
            "/usr/bin/perf",
            "perf",
            "record",
            "--output",
            perf_out_file.c_str(),
            "--pid",
            parent_pid.str().c_str(),
            "--event=cpu-clock,faults,cache-misses",
            "--call-graph=dwarf",
            "--control",
            control_ss.str().c_str(),
            "--delay=-1",
            nullptr
        ));
    } else { // If we're the parent process
        // Run the code to be profiled
        auto ctl_fifo_fd = open(ctl_fifo_name, O_WRONLY);
        auto ack_fifo_fd = open(ack_fifo_name, O_RDONLY);
        unlink(ctl_fifo_name); // Delete the FIFOs once this program crashes or exits.
        unlink(ack_fifo_name);
        char ack_fifo_buf[5];

        write(ctl_fifo_fd, "enable\n", 7);
        read(ack_fifo_fd, ack_fifo_buf, 5);
        KASSERT(std::string(ack_fifo_buf) == "ack\n", "Failed to enable perf: " << std::string(ack_fifo_buf));

        body();

        write(ctl_fifo_fd, "disable\n", 8);
        read(ack_fifo_fd, ack_fifo_buf, 5);
        KASSERT(std::string(ack_fifo_buf) == "ack\n", "Failed to disable perf: " << std::string(ack_fifo_buf));
        close(ctl_fifo_fd);
        close(ack_fifo_fd);

        // Kill the profiler
        kill(pid, SIGINT);
        waitpid(pid, nullptr, 0);
    }
}

inline void profile(std::function<void()> body) {
    profile("perf.data", body);
}
