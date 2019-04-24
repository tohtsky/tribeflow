#ifndef TRIBEFLOW_DEFINISIONS_HPP
#define TRIBEFLOW_DEFINISIONS_HPP

#include<Eigen/Eigen>
#include<map>
#include"stamp_lists.hpp"
#include"include/optional.hpp"
//using namespace Eigen;
using namespace std;

using DoubleMatrix = Eigen::Matrix<double, -1, -1, Eigen::RowMajor>;
using IntegerMatrix = Eigen::Matrix<int32_t, -1, -1, Eigen::RowMajor>;
using DoubleVector = Eigen::VectorXd;
using IntegerVector = Eigen::VectorXi;

template <typename T>
using optional = tl::optional<T>; 


//using IntegerVector = Eigen::Matrix<int32_t, -1, Eigen::Dynamic>;
//using IntegerMatrix = Eigen::Matrix<int32_t, -1, -1, Eigen::Dynamic>;
constexpr int DTS=0, TRACE=1, TRACE_HYPER=2, TRACE_TOPIC=3, 
          STAMP_LIST=4, COUNT_ZH=5, COUNT_SZ=6, HYPER2ID=7,
          SITE2ID=8;

struct HyperParams {
    inline HyperParams(size_t n_topics, size_t n_iter, size_t burn_in,
    bool dynamic, size_t n_batches, double alpha_zh,
    double beta_zs, string kernel_name, const vector<double> & residency_priors,
    int random_seed=42):
        n_topics(n_topics), n_iter(n_iter), burn_in(burn_in),
        dynamic(dynamic), n_batches(n_batches), alpha_zh(alpha_zh),
        beta_zs(beta_zs), kernel_name(kernel_name), residency_priors(residency_priors),
        random_seed(random_seed){
        }

    const size_t n_topics;
    const size_t n_iter;
    const size_t burn_in;
    const bool dynamic;
    const size_t n_batches;
    const double alpha_zh;
    const double beta_zs;
    const string kernel_name;
    const vector<double> residency_priors; 
    const int random_seed;
};

using InputData = std::tuple<
    DoubleMatrix, //Dts_mat
    IntegerMatrix, // Trace_mat
    vector<size_t>, // trace_hyper_ids
    IntegerVector, // trace_topics
    StampLists, //stamp_lists
    IntegerMatrix, // Count_zh
    IntegerMatrix, // Count_sz
    map<string, int>, // hyper2id
    map<string, int> //site2id
>;

struct OutPutData {
    inline OutPutData () {} // all default constructible
    size_t n_topics;
    double alpha_zh, beta_zs;
    vector<double> residency_priors; 
    size_t n_iter;
    size_t burn_in;
    DoubleMatrix Dts;
    IntegerMatrix Count_zh, Count_sz;
    IntegerVector count_h, count_z;
    DoubleMatrix Theta_zh, Psi_sz; 
    string kernel_name;
    DoubleMatrix P;
    IntegerVector assign;
    map<string, int> hyper2id, site2id;
};


#endif

