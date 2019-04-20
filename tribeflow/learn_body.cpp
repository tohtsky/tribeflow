#include"learn_body.hpp"

void fast_populate(
    const Eigen::MatrixXi & Trace,
    const Eigen::VectorXi & trace_hyper_ids,
    const Eigen::VectorXi & trace_topics, 
    Eigen::MatrixXi & Count_zh,
    Eigen::MatrixXi & Count_sz,
    Eigen::VectorXi & count_h,
    Eigen::VectorXi & count_z) {
    size_t n_trace = Trace.rows();
    for (size_t i = 0; i < n_trace; i++) {
        int h = trace_hyper_ids(i);
        int z = trace_topics(i);
        Count_zh(z, h)++;
        count_h(h)++;
        for(size_t j = 0; j < Trace.cols(); j++ ) {
            int s = Trace(i, j);
            Count_sz(s, z)++;
        }
        count_z(z) += Trace.cols(); 
    }
}
