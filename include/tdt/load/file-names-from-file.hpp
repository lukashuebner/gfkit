#pragma once

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

inline std::vector<std::string> file_names_from_file(std::string file_name) {
    auto is_whitespaces_only = [](std::string const& s) {
        return s.find_first_not_of(" \t") == std::string::npos;
    };

    auto is_comment = [](std::string const& s) {
        return s.find_first_not_of(" \t") == '#';
    };

    auto is_include = [](std::string const& s) {
        return s.find_first_not_of(" \t") == '@';
    };

    std::vector<std::string> file_names;
    std::ifstream            in_stream(file_name.c_str());

    if (!in_stream) {
        throw std::runtime_error("Could not open file " + file_name);
    }

    std::string str;
    while (std::getline(in_stream, str)) {
        if (is_include(str)) {
            auto const include_file_name       = str.substr(str.find_first_not_of(" \t@"));
            auto const file_names_from_include = file_names_from_file(include_file_name);
            file_names.insert(file_names.end(), file_names_from_include.begin(), file_names_from_include.end());
        } else if (is_whitespaces_only(str) || is_comment(str)) {
            continue;
        } else {
            file_names.push_back(str);
        }
    }

    in_stream.close();
    return file_names;
}
