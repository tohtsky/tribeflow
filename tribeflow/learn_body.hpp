#include <Eigen/Eigen>

void fast_populate(
    const Eigen::MatrixXi & Trace,
    const Eigen::VectorXi & trace_hyper_ids,
    const Eigen::VectorXi & trace_topics, 
    Eigen::MatrixXi & Count_zh,
    Eigen::MatrixXi & Count_sz,
    Eigen::VectorXi & count_h,
    Eigen::VectorXi & couht_z);
