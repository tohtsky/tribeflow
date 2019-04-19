#include "plearn.hpp"
#include "dataio.hpp"

int main () {
    string trace_fpath = "test_data/sampled_trace.data";
    unsigned n_workers = 4;
    unsigned n_topics = 50;
    string model_output_path = "model.data";
    unsigned n_iter = 20;
    unsigned burn_in = 5;
    unsigned dynamic = false;
    unsigned n_batches = 10;
    double alpha_zh = 50.0 / n_topics ;
    double beta_zs = 1e-3;
    string kernel = "noop";
    vector<double> residency_priors;

    double leaveout = 3e-1;

    unsigned n_traces = count_line(trace_fpath);

    MasterWorker worker(n_workers); 
    //worker.set_problem(
    //    trace_fpath, n_topics, alpha_zh, beta_zs, kernel, residency_priors,
    //    n_iter, from_, to_
    //)
    worker.create_slaves();
    worker.do_manage();
}
