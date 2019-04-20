#ifndef TRIBEFLOW_DEFINISIONS_HPP
#define TRIBEFLOW_DEFINISIONS_HPP

#include<Eigen/Eigen>
#include<map>
#include"stamp_lists.hpp"
//using namespace Eigen;
using namespace std;

using IntegerVector = Eigen::Matrix<int32_t, -1, Eigen::Dynamic>;
//using IntegerMatrix = Eigen::Matrix<int32_t, -1, -1, Eigen::Dynamic>;
constexpr int DTS=0, TRACE=1, TRACE_HYPER=2, TRACE_TOPIC=3, 
          STAMP_LIST=4, COUNT_ZH=5, COUNT_sz=6, HYPER2ID=7,
          SITE2ID=8;

using HyperParams = std::tuple<
    size_t, // n_topics
    size_t, // n_iter
    size_t, // burn_in
    bool, // dynamic
    size_t, // n_batches,
    double, // alpha_zh
    double, // beta_zs
    string, // kernel_name
    vector<double> // residency_priors
>;


using InputData = std::tuple<
    Eigen::MatrixXd, //Dts_mat
    Eigen::MatrixXi, // Trace_mat
    Eigen::VectorXi, // trace_hyper_ids
    Eigen::VectorXi, // trace_topics
    StampLists, //stamp_lists
    Eigen::MatrixXi, // Count_zh
    Eigen::MatrixXi, // Count_sz
    map<string, int>, // hyper2id
    map<string, int> //site2id
>;


#endif

