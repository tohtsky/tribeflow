#include<iostream>
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

void e_step(
    const Eigen::MatrixXd& Dts,
    const Eigen::MatrixXi& Trace,
    const Eigen::VectorXi& trace_hyper_ids,
    Eigen::VectorXi& trace_topics, // change
    const StampLists& previous_stamps,
    Eigen::MatrixXi& Count_zh, // change!
    Eigen::MatrixXi& Count_sz, // change
    Eigen::VectorXi& count_h, //change
    Eigen::VectorXi& count_z,  //change
    double alpha_zh,
    double beta_zs,
    vector<double>& prob_topics_aux,
    std::unique_ptr<KernelBase> && kernel,
    std::mt19937 & gen) {

    size_t mem_size = Dts.cols();
    size_t n_trace = Trace.rows();

    for (size_t i = 0; i < n_trace; i++ ) {
        double dt = Dts(i,  Dts.cols() -1);
        int hyper_iid = trace_hyper_ids(i);
        int topic_old = trace_topics(i);
        Count_zh(topic_old, hyper_iid)--;
        count_h(hyper_iid) --;
        //Count_sz(Trace(i).transpose().array(), topic_old) -= 1;
        //count_z(Trace(i).array().transpose()) -= 1;
        count_z(topic_old) -= (mem_size + 1); 
        for (size_t j = 0; j < (mem_size + 1); j++) { 
            int s = Trace(i, j);
            Count_sz(s, topic_old)--;
        }
        size_t new_topic = sample(i, dt, Trace, hyper_iid,
                trace_topics, previous_stamps, Count_zh,
            Count_sz, count_h, count_z, alpha_zh, beta_zs,
            prob_topics_aux, std::move(kernel), gen);
        //cout << "new topic is " << new_topic << endl;
        trace_topics(i) = new_topic;
        Count_zh(new_topic, hyper_iid) ++;
        count_h(hyper_iid) ++;

        for ( size_t j = 0; j < (mem_size + 1); j++ ) {
            int site = Trace(i, j);
            Count_sz(site, new_topic) += (mem_size + 1); 
        }
        count_z(new_topic) += (mem_size + 1); 
    }
}

void m_step(
    const Eigen::MatrixXd& Dts,
    const Eigen::MatrixXi& Trace,
    const Eigen::VectorXi& trace_hyper_ids,
    const Eigen::VectorXi& trace_topics, 
    StampLists& previous_stamps,
    std::unique_ptr<KernelBase> && kernel,
    std::mt19937 & gen
) {
    for (size_t i = 0; i < previous_stamps.size(); i++) {
        previous_stamps[i].resize(0);
    }
    for (size_t i = 0; i < Trace.rows(); i++) {
        double dt = Dts(i, Dts.cols() - 1); 
        size_t topic = trace_topics(i);
        previous_stamps.at(topic).push_back(dt);
    }
    kernel -> m_step(previous_stamps);
} 

inline void aggregate(
    const Eigen::MatrixXi & Count_zh, 
    const Eigen::MatrixXi & Count_sz, 
    const Eigen::VectorXi & count_h,
    const Eigen::VectorXi & count_z,
    double alpha_zh, double beta_zs,
    Eigen::MatrixXd & Theta_zh, 
    Eigen::MatrixXd & Psi_sz) {
    int nz = Theta_zh.rows();
    int nh = Theta_zh.cols();
    int ns = Psi_sz.rows();
    for (int z = 0; z < nz; z++) {
        for (int h = 0; h < nh; h++ ){
            Theta_zh(z, h) += dir_posterior(Count_zh(z, h), count_h(h), nz, alpha_zh);
        }

        for (int s = 0; s < ns; s++) {
            Psi_sz(s, z) += dir_posterior(Count_sz(s, z), count_z(z), ns, beta_zs);
        }
    }
}

void fast_em(const Eigen::MatrixXd & Dts, const Eigen::MatrixXi & Trace,
        const Eigen::VectorXi & trace_hyper_ids,
        Eigen::VectorXi& trace_topics, // change 
        StampLists & previous_stamps, 
        Eigen::MatrixXi & Count_zh, Eigen::MatrixXi & Count_sz,
        Eigen::VectorXi & count_h, Eigen::VectorXi & count_z,
        double alpha_zh, double beta_zs, vector<double> & prob_topics_aux,
        Eigen::MatrixXd & Theta_zh, Eigen::MatrixXd Psi_sz, size_t n_iter,
        size_t burn_in, std::unique_ptr<KernelBase> && kernel, std::mt19937 & gen) {
    for (size_t i = 0; i < n_iter; i++) {
        e_step(Dts, Trace, trace_hyper_ids, trace_topics, previous_stamps,
                Count_zh, Count_sz, count_h, count_z,
                alpha_zh, beta_zs, prob_topics_aux, std::move(kernel), gen);
        m_step(Dts, Trace, trace_hyper_ids, trace_topics, previous_stamps, 
                std::move(kernel), gen);

        if (i >= burn_in) {
            aggregate(Count_zh, Count_sz, count_h, count_z,
                alpha_zh, beta_zs, Theta_zh, Psi_sz);
        }

    }
}

void col_normalize(Eigen::MatrixXd & target) {
    Eigen::ArrayXd norm = target.colwise().sum().array().transpose();
    norm = (norm == 0.0).select(1, norm);
    target.array().rowwise() /= norm.transpose();
}

void em(const Eigen::MatrixXd & Dts, const Eigen::MatrixXi & Trace,
        const Eigen::VectorXi & trace_hyper_ids,
        Eigen::VectorXi& trace_topics, // change 
        StampLists & previous_stamps, 
        Eigen::MatrixXi & Count_zh, Eigen::MatrixXi & Count_sz,
        Eigen::VectorXi & count_h, Eigen::VectorXi & count_z,
        double alpha_zh, double beta_zs, vector<double> & prob_topics_aux,
        Eigen::MatrixXd & Theta_zh, Eigen::MatrixXd Psi_sz, size_t n_iter,
        size_t burn_in, std::unique_ptr<KernelBase> && kernel, std::mt19937 & gen,
        bool average_and_normalize) {
    fast_em(Dts, Trace, trace_hyper_ids, trace_topics, 
        previous_stamps, 
        Count_zh, Count_sz, count_h, count_z,
        alpha_zh, beta_zs, prob_topics_aux,
        Theta_zh, Psi_sz, n_iter,
        burn_in, std::move(kernel), gen);

    if (average_and_normalize) {
        if (n_iter > burn_in) {
            Theta_zh /= (static_cast<int>(n_iter) - static_cast<int>(burn_in));
            Psi_sz /= (static_cast<int>(n_iter) - static_cast<int>(burn_in)); 
        }
        col_normalize(Theta_zh);
        col_normalize(Psi_sz);
    }
}



size_t sample(
    size_t i,
    double dt,
    const Eigen::MatrixXi& Trace,
    int hyper_id,
    Eigen::VectorXi& trace_topics, // change
    const StampLists& previous_stamps,
    Eigen::MatrixXi& Count_zh, // change!
    Eigen::MatrixXi& Count_sz, // change
    Eigen::VectorXi& count_h, //change
    Eigen::VectorXi& count_z,  //change
    double alpha_zh,
    double beta_zs,
    vector<double> & prob_topics_aux,
    std::unique_ptr<KernelBase> && kernel,
    std::mt19937 & gen
    ) {

    size_t nz = Count_zh.rows();
    size_t ns = Count_sz.rows();
    int prev;
    int site;
    double prev_prob;
    
    for (size_t z = 0; z < nz; z++ ) {
        site = Trace(i, 0);
        prob_topics_aux[z] = kernel->pdf(dt, z, previous_stamps) * 
            dir_posterior(Count_zh(z, hyper_id), count_h(hyper_id), nz, alpha_zh) *
            dir_posterior(Count_sz(site, z), count_z(z), ns, beta_zs);
        for(size_t j = 1; j < Trace.cols(); j++ ) {
            site = Trace(i, j);
            prev = Trace(i, j - 1);
            prev_prob = dir_posterior(Count_sz(prev, z), count_z(z), ns, beta_zs);
            prob_topics_aux[z] = prob_topics_aux[z] * 
                dir_posterior(Count_sz(site, z), count_z(z), ns, beta_zs ) /
                (1 - prev_prob); 
        } 
        //if (z > 0) {
        //    prob_topics_aux(z) += prob_topics_aux(z-1);
        //}
    }
    return std::discrete_distribution(
        prob_topics_aux.begin(), 
        prob_topics_aux.end()
    )(gen);
}


