#include<iostream>
#include"learn_body.hpp"

void fast_populate(
    const IntegerMatrix & Trace,
    const IntegerVector & trace_hyper_ids,
    const IntegerVector & trace_topics, 
    IntegerMatrix & Count_zh,
    IntegerMatrix & Count_sz,
    IntegerVector & count_h,
    IntegerVector & count_z) {
    size_t n_trace = Trace.rows();
    size_t trace_col = Trace.cols();
    for (size_t i = 0; i < n_trace; i++) {
        int h = trace_hyper_ids(i);
        int z = trace_topics(i);
        Count_zh(z, h)++;
        count_h(h)++;
        for(size_t j = 0; j < trace_col; j++ ) {
            int s = Trace(i, j);
            Count_sz(s, z)++;
        }
        count_z(z) += Trace.cols(); 
    }
}

void e_step(
    const DoubleMatrix& Dts,
    const IntegerMatrix& Trace,
    const IntegerVector& trace_hyper_ids,
    IntegerVector& trace_topics, // change
    const StampLists& previous_stamps,
    IntegerMatrix& Count_zh, // change!
    IntegerMatrix& Count_sz, // change
    IntegerVector& count_h, //change
    IntegerVector& count_z,  //change
    double alpha_zh,
    double beta_zs,
    vector<double>& prob_topics_aux,
    std::unique_ptr<KernelBase> && kernel,
    std::mt19937 & gen) {

    size_t mem_size = Dts.cols();
    size_t n_trace = Trace.rows();

    size_t Dt_last_col = Dts.cols() - 1;
    for (size_t i = 0; i < n_trace; i++ ) {
        double dt = Dts(i,  Dt_last_col);
        int hyper_iid = trace_hyper_ids(i);
        int topic_old = trace_topics(i);
        Count_zh(topic_old, hyper_iid)--;
        count_h(hyper_iid) --;
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

        for (size_t j = 0; j < (mem_size + 1); j++) {
            int site = Trace(i, j);
            Count_sz(site, new_topic) += 1; 
        }
        count_z(new_topic) += (mem_size + 1); 
    }
}

void m_step(
    const DoubleMatrix& Dts,
    const IntegerMatrix& Trace,
    const IntegerVector& trace_hyper_ids,
    const IntegerVector& trace_topics, 
    StampLists& previous_stamps,
    std::unique_ptr<KernelBase> && kernel,
    std::mt19937 & gen
) {
    for (size_t i = 0; i < previous_stamps.size(); i++) {
        previous_stamps[i].resize(0);
    }
    size_t n_trace = Trace.rows();
    size_t Dt_col_last = Dts.cols() - 1;
    for (size_t i = 0; i < n_trace; i++) {
        double dt = Dts(i, Dt_col_last); 
        size_t topic = trace_topics(i);
        previous_stamps.at(topic).push_back(dt);
    }
    kernel -> m_step(previous_stamps);
} 

inline void aggregate(
    const IntegerMatrix & Count_zh, 
    const IntegerMatrix & Count_sz, 
    const IntegerVector & count_h,
    const IntegerVector & count_z,
    double alpha_zh, double beta_zs,
    DoubleMatrix & Theta_zh, 
    DoubleMatrix & Psi_sz) {
    size_t nz = Theta_zh.rows();
    size_t nh = Theta_zh.cols();
    size_t ns = Psi_sz.rows();
    for (size_t z = 0; z < nz; z++) {
        for (size_t h = 0; h < nh; h++ ){
            Theta_zh(z, h) += dir_posterior(Count_zh(z, h), count_h(h), nz, alpha_zh);
        }

        for (size_t s = 0; s < ns; s++) {
            Psi_sz(s, z) += dir_posterior(Count_sz(s, z), count_z(z), ns, beta_zs);
        }
    }
}

void fast_em(const DoubleMatrix & Dts, const IntegerMatrix & Trace,
        const IntegerVector & trace_hyper_ids,
        IntegerVector& trace_topics, // change 
        StampLists & previous_stamps, 
        IntegerMatrix & Count_zh, IntegerMatrix & Count_sz,
        IntegerVector & count_h, IntegerVector & count_z,
        double alpha_zh, double beta_zs, vector<double> & prob_topics_aux,
        DoubleMatrix & Theta_zh, DoubleMatrix Psi_sz, size_t n_iter,
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

void col_normalize(DoubleMatrix & target) {
    Eigen::ArrayXd norm = target.colwise().sum().array().transpose();
    norm = (norm == 0.0).select(1, norm);
    target.array().rowwise() /= norm.transpose();
}

void em(const DoubleMatrix & Dts, const IntegerMatrix & Trace,
        const IntegerVector & trace_hyper_ids,
        IntegerVector& trace_topics, // change 
        StampLists & previous_stamps, 
        IntegerMatrix & Count_zh, IntegerMatrix & Count_sz,
        IntegerVector & count_h, IntegerVector & count_z,
        double alpha_zh, double beta_zs, vector<double> & prob_topics_aux,
        DoubleMatrix & Theta_zh, DoubleMatrix Psi_sz, size_t n_iter,
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
    std::mt19937 & gen
    ) {

    size_t nz = Count_zh.rows();
    size_t ns = Count_sz.rows();
    int prev;
    int site;
    double prev_prob;

    size_t trace_cols = Trace.cols(); 
    for (size_t z = 0; z < nz; z++ ) {
        site = Trace(i, 0);
        prob_topics_aux[z] = kernel->pdf(dt, z, previous_stamps) * 
            dir_posterior(Count_zh(z, hyper_id), count_h(hyper_id), nz, alpha_zh) *
            dir_posterior(Count_sz(site, z), count_z(z), ns, beta_zs);
        for(size_t j = 1; j < trace_cols; j++ ) {
            site = Trace(i, j);
            prev = Trace(i, j - 1);
            prev_prob = dir_posterior(Count_sz(prev, z), count_z(z), ns, beta_zs);
            prob_topics_aux[z] = prob_topics_aux[z] * 
                dir_posterior(Count_sz(site, z), count_z(z), ns, beta_zs ) /
                (1 - prev_prob); 
        } 
    }
    return std::discrete_distribution(
        prob_topics_aux.begin(), 
        prob_topics_aux.end()
    )(gen);
}


