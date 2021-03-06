#include "plearn.hpp"
#include "dataio.hpp"
#include "debug.hpp"
#include "learn_body.hpp"

int main () {
    string trace_fpath = "/Users/tomoki.otsuki/src/tribeflow/tribeflow/test_data/sampled_trace.data";
    size_t n_workers = 1;
    unsigned n_topics = 20;
    string model_output_path = "model.data";
    unsigned n_iter = 2000;
    unsigned burn_in = 0;
    bool dynamic = false;
    unsigned n_batches = 10;
    double alpha_zh = 50.0 / n_topics ;
    double beta_zs = 1e-3;
    string kernel_name = "noop";
    vector<double> residency_priors;

    HyperParams hyper_params ( 
        n_topics, n_iter, burn_in, dynamic, n_batches,
        alpha_zh, beta_zs, kernel_name, residency_priors
    );


    MasterWorker worker(n_workers, hyper_params, initialize_trace(
        trace_fpath, n_topics, n_iter, 0, tl::nullopt, tl::nullopt
    ));
    worker.create_slaves();
    OutPutData result = worker.do_manage();
}
