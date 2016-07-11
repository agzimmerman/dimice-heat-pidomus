#include "test_tools.h"
#include <iostream>
#include <fstream>
bool are_files_equal(const std::string filepath_a, const std::string filepath_b) {
    std::ifstream file_stream;
    file_stream.open(filepath_a, std::ifstream::in);
    std::string string_a( (std::istreambuf_iterator<char>(file_stream) ), (std::istreambuf_iterator<char>() ) );
    file_stream.close();
    file_stream.open(filepath_b, std::ifstream::in);
    std::string string_b( (std::istreambuf_iterator<char>(file_stream) ), (std::istreambuf_iterator<char>() ) );
    file_stream.close();
    bool are_equal = string_a == string_b;
    return are_equal;
}
