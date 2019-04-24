#ifndef DATAIO_HPP
#define DATAIO_HPP

#include"defs.hpp"

#include<Eigen/Eigen>
#include<tuple>
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include "stamp_lists.hpp"


size_t count_line(string file_path);

InputData
initialize_trace(string trace_fpath, size_t n_topics, size_t num_iter,
        size_t from_, tl::optional<size_t> to, tl::optional<vector<int>> initial_assign, 
        int random_seed=42);

#endif
