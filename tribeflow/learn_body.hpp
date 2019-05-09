#ifndef LEARN_BODY_HPP
#define LEARN_BODY_HPP
#include <Eigen/Eigen>
#include "defs.hpp"
#include "kernels/base.hpp"
#include <random>

constexpr double REGULATOR = 1e-10;

void fast_populate(
    const IntegerMatrix & Trace,
    const vector<size_t> & trace_hyper_ids,
    const IntegerVector & trace_topics, 
    IntegerMatrix & Count_zh,
    IntegerMatrix & Count_sz,
    IntegerVector & count_h,
    IntegerVector & count_z);

void m_step(
    const DoubleMatrix& Dts,
    const IntegerMatrix& Trace,
    const vector<size_t>& trace_hyper_ids,
    const IntegerVector& trace_topics, 
    StampLists& previous_stamps,
    std::unique_ptr<KernelBase> && kernel,
    std::mt19937 & gen
);

void e_step(
    const DoubleMatrix& Dts,
    const IntegerMatrix& Trace,
    const vector<size_t>& trace_hyper_ids,
    IntegerVector& trace_topics, // change
    const StampLists& previous_stamps,
    IntegerMatrix& Count_zh, // change!
    IntegerMatrix& Count_sz, // change
    IntegerVector& count_h, //change
    IntegerVector& count_z,  //change
    double alpha_zh,
    double beta_zs,
    vector<double> & prob_topics_aux, // change
    std::unique_ptr<KernelBase> && kernel,
    std::mt19937 & gen);

inline double dir_posterior(double joint_count, double global_count, 
        double num_occurences, double smooth) { 
    double numerator = joint_count + smooth ;
    double denominator = global_count + smooth * num_occurences;
    if (denominator == 0) return 0;
    return numerator / denominator;
}

inline Eigen::ArrayXd dir_posterior(
        const IntegerVector & joint_count, 
        const IntegerVector & global_count, 
        double num_occurences, double smooth) { 
    DoubleVector num =  (joint_count.cast<double>().array() + smooth);
    DoubleVector denom = (global_count.cast<double>().array() + smooth * num_occurences + REGULATOR);
    return num.array() / denom.array();
}

void fast_em(const DoubleMatrix & Dts, const IntegerMatrix & Trace,
        const vector<size_t> & trace_hyper_ids,
        IntegerVector& trace_topics, // change 
        StampLists & previous_stamps, 
        IntegerMatrix & Count_zh, IntegerMatrix & Count_sz,
        IntegerVector & count_h, IntegerVector & count_z,
        double alpha_zh, double beta_zs, vector<double> & prob_topics_aux,
        DoubleMatrix & Theta_zh, DoubleMatrix & Psi_sz, size_t n_iter,
        size_t burn_in, std::unique_ptr<KernelBase> && kernel,
        std::mt19937 & gen);

void em(const DoubleMatrix & Dts, const IntegerMatrix & Trace,
        const vector<size_t> & trace_hyper_ids,
        IntegerVector& trace_topics, // change 
        StampLists & previous_stamps, 
        IntegerMatrix & Count_zh, IntegerMatrix & Count_sz,
        IntegerVector & count_h, IntegerVector & count_z,
        double alpha_zh, double beta_zs, vector<double> & prob_topics_aux,
        DoubleMatrix & Theta_zh, DoubleMatrix & Psi_sz, size_t n_iter,
        size_t burn_in, std::unique_ptr<KernelBase> && kernel, std::mt19937 & gen,
        bool average_and_normalize=true);

void col_normalize(DoubleMatrix & target);

size_t sample(
    size_t i,
    double dt,
    const IntegerMatrix& Trace,
    int hyper_id,
    IntegerVector& trace_topics, // change
    const StampLists& previous_stamps,
    IntegerMatrix& Count_zh, // change!
    IntegerMatrix& Count_sz, // change
    IntegerVector& count_h, //change
    IntegerVector& count_z,  //change
    double alpha_zh,
    double beta_zs,
    vector<double> & prob_topics_aux,
    std::unique_ptr<KernelBase> && kernel,
    std::mt19937 & gen);

void aggregate(
        const IntegerMatrix & Count_zh, 
    const IntegerMatrix & Count_sz, 
    const IntegerVector & count_h,
    const IntegerVector & count_z,
    double alpha_zh, double beta_zs,
    DoubleMatrix & Theta_zh, 
    DoubleMatrix & Psi_sz); 

inline int binary_search(const vector<double> & array, double value) {
    int lower = 0, upper = array.size() -1; // closed interval
    int half = 0;
    int idx = -1;
    while (upper >= lower){
        half = lower + (upper - lower) / 2;
        double trial = array[half];
        if (value == trial) {
            idx = half;
            break;
        }
        else if ( value > trial ) {
            lower = half + 1;
        } else {
            upper = half - 1;
        }
    }
    if (idx == - 1) // Element not found, return where it should be
        return lower;
    return idx;
}


#endif
