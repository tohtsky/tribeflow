#ifndef DATAIO_HPP
#define DATAIO_HPP

#include<Eigen/Eigen>
#include<optional>
#include<tuple>
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>
#include "stamp_lists.hpp"

using namespace std;
using namespace Eigen;

size_t count_line(string file_path);

using InitializeTraceReturnType = std::tuple<
    Eigen::MatrixXd, //Dts_mat
    Eigen::MatrixXi, // Trace_mat
    StampLists, //stamp_lists
    Eigen::MatrixXi, // Count_zh
    Eigen::MatrixXi, // Count_oz
    Eigen::VectorXi, // count_h
    Eigen::VectorXi, // count_z
    Eigen::VectorXd, // prob_topics_aux
    Eigen::MatrixXd, // Theta_zh
    Eigen::MatrixXd, // Psi_oz
    map<string, int>, // hyper2id
    map<string, int> //obj2id
>;

InitializeTraceReturnType
initialize_trace(string trace_fpath, size_t n_topics, size_t num_iter,
        size_t from_, optional<size_t> to, optional<vector<int>> initial_assign, 
        int random_seed=42);

#endif
