#ifndef DATAIO_HPP
#define DATAIO_HPP

#include<Eigen/Eigen>
#include<optional>
#include<tuple>
#include<string>
#include<fstream>

using namespace std;
using namespace Eigen;
//using InitializedData = tuple<>;

size_t count_line(string file_path);

void initialize_trace(string trace_fpath, size_t n_topics, size_t num_iter,
        size_t from_, optional<size_t> to, optional<MatrixXi> initial_assign);

#endif
