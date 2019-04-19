#include "dataio.hpp"

size_t count_line(string file_path) {
    std::ifstream ifs(file_path); 
    return std::count(std::istreambuf_iterator<char>(ifs), 
            std::istreambuf_iterator<char>(), '\n') + 1;
}
