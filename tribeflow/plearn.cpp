#define DEBUG

#include "plearn.hpp"
#include "learn_body.hpp"
#include "kernels/base.hpp"
#include "debug.hpp"

constexpr size_t PRINT_FREQ = 50;

constexpr size_t CACHE_SIZE = 1;
MasterWorker::Slave::Slave(MasterWorker* parent, size_t id) :
    parent_(parent), my_id(id), slave_mutex(new std::mutex), 
    slave_condition(new std::condition_variable){
    }


void MasterWorker::Slave::start_work () {
    this->working_thread = std::thread(
        [this] {
            this->main_job();
        });
}

void MasterWorker::Slave::join() {
    this->working_thread->join();
}

void MasterWorker::Slave::set_message_from_master(SlaveStatus msg) {
    message_from_master = msg;
    this->slave_condition->notify_one();
}

void MasterWorker::Slave::learn () {
    const vector<vector<int>> & workloads = this->parent_->workloads();
    const vector<int> & relevant_trace_ind =  workloads.at(my_id);

    const DoubleMatrix & Dts_ref = this->parent_->Dts();
    const IntegerMatrix & Trace_ref = this->parent_->Trace();
    const vector<size_t> & trace_hyper_ids_ref = this->parent_->trace_hyper_ids();
    const IntegerMatrix & trace_topics_ref = this->parent_->trace_topics(); 

    auto kernel = kernel_factory(this->parent_->hyper_params.kernel_name);

    DoubleMatrix Dts(relevant_trace_ind.size(), Dts_ref.cols()); 
    IntegerMatrix Trace(relevant_trace_ind.size(), Trace_ref.cols());
    vector<size_t> trace_hyper_ids(relevant_trace_ind.size());
    IntegerVector trace_topics(relevant_trace_ind.size()); 

    for(size_t i = 0; i < relevant_trace_ind.size(); i++) {
        size_t orig_ind = relevant_trace_ind[i];
        Dts.row(i) =  Dts_ref.row(orig_ind);
        Trace.row(i) = Trace_ref.row(orig_ind); 
        trace_hyper_ids[i] = trace_hyper_ids_ref[orig_ind];
        trace_topics(i) = trace_topics_ref(orig_ind);
    }

    kernel->build(
            Trace.rows(), this->parent_->Count_zh().rows(),
            this->parent_->hyper_params.residency_priors);

    DoubleMatrix P(this->parent_->P());
    DoubleMatrix P_pair(P);
    kernel->update_state(P);


    std::mt19937 gen(this->parent_->hyper_params.random_seed + my_id);
    const double alpha_zh = this->parent_->hyper_params.alpha_zh;
    const double beta_zs = this->parent_->hyper_params.beta_zs;

    const int nz = this->parent_->Count_zh().rows();
    const int nh = this->parent_->hyper2id().size();
    const int ns = this->parent_->site2id().size();
    int mem_size = this->parent_->Dts().cols();
    if (mem_size <=0) {
        throw std::runtime_error("mem_size <=0 somehow ");
    }
    IntegerMatrix Count_zh = IntegerMatrix::Zero(nz, nh);
    IntegerMatrix Count_sz = IntegerMatrix::Zero(ns, nz);
    IntegerVector count_h = IntegerVector::Zero(nh);
    IntegerVector count_z = IntegerVector::Zero(nz);

    fast_populate(
        Trace, trace_hyper_ids, trace_topics,
        Count_zh, Count_sz, count_h, count_z
    );


    // sample function in plearn.py in original Python version.

    for (size_t worker_id = 0; worker_id < parent_->n_slaves(); worker_id++) {
        previous_encounters_s.insert(
            {worker_id, IntegerMatrix::Zero(ns, nz)}
        );
    }

    StampLists stamps(nz);
    size_t n_trace = Trace.rows();
    for (size_t i = 0; i < n_trace; i++) {
        auto z = trace_topics(i);
        stamps.at(z).push_back(Dts(i, Dts.cols() - 1));
    }

    vector<double> prob_topics_aux(nz, 0.0);

    IntegerMatrix Count_sz_others = IntegerMatrix::Zero(ns, nz);
    IntegerMatrix Count_sz_sum = IntegerMatrix::Zero(ns, nz); 

    DoubleMatrix Theta_zh = DoubleMatrix::Zero(nz, nh);
    DoubleMatrix Psi_sz = DoubleMatrix::Zero(ns, nz);
    bool can_pair = true;
    for (size_t i = 0; i < (parent_->hyper_params.n_iter / CACHE_SIZE) ; i++) {
        if ( i > 0 && (i % PRINT_FREQ == 0) ){ 
            this->print(str_concat("Current iteration: ", i));
        }
        Count_sz_sum = Count_sz + Count_sz_others;
        count_z = Count_sz_sum.colwise().sum().transpose();
        // em
        em(
            Dts, Trace, trace_hyper_ids, trace_topics,
            stamps, Count_zh, Count_sz_sum, count_h, count_z,
            alpha_zh, beta_zs, prob_topics_aux,
            Theta_zh, Psi_sz, CACHE_SIZE, 0,  std::move(kernel), 
            gen, false
        );

        // update local counts
        Count_sz = Count_sz_sum - Count_sz_others;
        count_z = Count_sz.colwise().sum().transpose(); 

        // Update expected belief of other processors 
        if (can_pair) {
            P = kernel-> get_state();
            send_data_count = & Count_sz;
            send_data_p = & P;
            can_pair = paired_update(Count_sz, Count_sz_others, P);
            kernel->update_state(P);
        }
    }
    this->result = ResultData{trace_topics, Count_zh, Count_sz, count_h, count_z, P};
}

bool MasterWorker::Slave::paired_update(
        IntegerMatrix & Count_sz, 
        IntegerMatrix & Count_sz_others, 
        DoubleMatrix & P_local) {
    send_message_to_master(SlaveStatus::PAIRME);
    SlaveStatus am_i_paired = receive_message_from_master();
    if (am_i_paired != SlaveStatus::PAIRED) 
        throw std::runtime_error("Error in control flow...");
    if (!pair_id_) {
        throw std::runtime_error("Pair id not set..."); 
    }

    if (*pair_id_ == my_id) {
        this->pair_id_ = std::nullopt;
        this->received_data = std::nullopt; 
        return false;
    }

    if (!received_data) {
        throw std::runtime_error("Exchanged data not set...");
    }

    const IntegerMatrix & Count_sz_pair = std::get<0>(*received_data);
    const DoubleMatrix & P_pair = std::get<1>(*received_data);

    Count_sz_others += Count_sz_pair;
    Count_sz_others -= this->previous_encounters_s.at(*pair_id_);
    this->previous_encounters_s.at(*pair_id_) = Count_sz_pair;
    P_local += P_pair;
    P_local /= 2;

    this->pair_id_ = std::nullopt;
    this->received_data = std::nullopt;
    return true;
}


void MasterWorker::Slave::main_job () {
    while (true) {
        auto event = receive_message_from_master();
        if (event == SlaveStatus::LEARN ) {
            send_message_to_master(SlaveStatus::STARTED);
            learn();
            send_message_to_master(SlaveStatus::FINISHED);
            break;
        } else {
            throw std::runtime_error("controll should not reach here..");
        }
    }
};

void MasterWorker::Slave::send_message_to_master (SlaveStatus status_code) {
    this->parent_->add_message(MasterMessage{this->my_id, status_code});
};

void MasterWorker::Slave::print(const std::string & message) {
    this->parent_->print_with_lock(this->my_id, message);
}

void MasterWorker::Slave::sample() {
    //bool can_pair = true;
}


SlaveStatus MasterWorker::Slave::receive_message_from_master () {
    std::unique_lock lock(*(this->slave_mutex));
    this->slave_condition->wait(lock, [this]{return message_from_master;});
    auto result = *message_from_master;
    message_from_master = std::nullopt;
    return result;
}

void MasterWorker::Slave::set_pair_id(size_t pair_id) {
    if (this->pair_id_) 
        throw std::runtime_error("pair_id is already set, and tried to overload..");
    this->pair_id_ = pair_id;
}

void MasterWorker::Slave::set_data_other_thread(const InterThreadData & other_data) {
    this->received_data = other_data;
}


MasterWorker::MasterWorker (size_t n_slaves, HyperParams hyper_params,
        InputData && input_data):
    n_slaves_(n_slaves), hyper_params(hyper_params),
    input_data(input_data), gen(hyper_params.random_seed), workloads_(n_slaves),
    kernel_(kernel_factory(hyper_params.kernel_name)) {

        kernel_->build(
            get<TRACE>(input_data).rows(),
            get<COUNT_ZH>(input_data).rows(),
            hyper_params.residency_priors
        );
    }

MasterWorker::~MasterWorker () {
    for(auto & slave: slaves) {
        slave.join();
    }
}

void MasterWorker::create_slaves () {
    size_t nh = hyper2id().size();
    // size_t ns = site2id().size();
    // size_t nz = Count_zh().rows();
    size_t n_paths = Trace().rows();
    vector<int> hyper_ids(nh);
    for (size_t i = 0; i < nh; i++) hyper_ids[i] = i;

    std::shuffle(hyper_ids.begin(), hyper_ids.end(), gen);
    vector<int> hyper_to_slave_id(nh);
    for (size_t i = 0; i < nh; i++) {
        size_t slave_id = i % this->n_slaves_;
        hyper_to_slave_id[hyper_ids[i]] = slave_id;
    }

    const vector<size_t> trace_hyper_ids_ref = trace_hyper_ids();
    for (size_t i = 0; i < n_paths; i++) {
        int h = trace_hyper_ids_ref[i];
        size_t slave_id = hyper_to_slave_id[h];
        workloads_[slave_id].push_back(i);
    }

    for(size_t i = 0; i < n_slaves_; i++) {
        slaves.emplace_back(
            this, i
        );
    }
    //LearningData workload = "Debug message";
    for(size_t i = 0; i < n_slaves_; i++) {
        slaves[i].start_work();
        slaves[i].set_message_from_master(SlaveStatus::LEARN);
    }
}

OutPutData MasterWorker::do_manage () {
    int available_to_pair = -1;
    unsigned n_finished = 0;
    vector<bool> finished(n_slaves_, false);

    while (n_finished != n_slaves_) {
        int worker_id;
        SlaveStatus what;
        std::tie(worker_id, what) = receive_message();
        switch (what) {
            case SlaveStatus::STARTED:
                cout << "Worker[" << worker_id << "] has started it's job." << endl;
                break;
            case SlaveStatus::PAIRME:
                if (n_finished == (n_slaves_ -1) ) {
                    // reserve it!
                    slaves[worker_id].set_pair_id(worker_id); 
                    slaves[worker_id].set_message_from_master(
                        SlaveStatus::PAIRED
                    );
                    available_to_pair = -1;
                } else if (available_to_pair < 0) {
                    available_to_pair = worker_id;
                } else {
                    // paired!
                    if (worker_id == available_to_pair) {
                        std::runtime_error("Some bug must have happened.");
                    }
                    slaves[worker_id].set_pair_id(available_to_pair); 
                    slaves[worker_id].set_exchange_data(
                        slaves[available_to_pair].exchanged_data()
                    );

                    slaves[available_to_pair].set_pair_id(worker_id);
                    slaves[available_to_pair].set_exchange_data(
                        slaves[worker_id].exchanged_data()
                    );

                    slaves[worker_id].set_message_from_master(
                        SlaveStatus::PAIRED
                    );
                    slaves[available_to_pair].set_message_from_master(
                        SlaveStatus::PAIRED
                    );

                    available_to_pair = -1;
                }
                break;
            case SlaveStatus::FINISHED:
                cout << "Worker[" << worker_id <<
                    "] has finished it's iterations." << endl;
                finished[worker_id] = true;
                if (!slaves[worker_id].result) {
                    throw std::runtime_error("A slave says it's done, but result not set..");
                }
                slave_results[worker_id] = *(slaves[worker_id].result);
                n_finished++;
                // The last worker who has requested pairing should be paired with itself.
                if ( (n_finished == n_slaves_ -1) && available_to_pair != -1) {
                    slaves[available_to_pair].set_pair_id(available_to_pair);
                    slaves[available_to_pair].set_message_from_master(
                        SlaveStatus::PAIRED
                    );
                }
                break;
            default:
                throw std::runtime_error(
                        "This message is inappropriate to be received"
                        );
        }
    }
    IntegerMatrix Count_zh = IntegerMatrix::Zero(
        this->Count_zh().rows(), this->Count_zh().cols()
    );
    IntegerMatrix Count_sz = IntegerMatrix::Zero(
        this->Count_sz().rows(), this->Count_sz().cols() 
    );

    IntegerVector count_h = IntegerVector::Zero(
        this->hyper2id().size()
    );
    IntegerVector count_z = IntegerVector::Zero(
        this->Count_zh().rows()
    );

    DoubleMatrix P = kernel_->get_state();
    P.setZero(P.rows(), P.cols());
    IntegerVector trace_topics(Trace().rows());
    for (auto iter = slave_results.begin(); iter != slave_results.end(); iter++) {
        const size_t worker_id = iter->first;
        {
            const IntegerVector & trace_topics_slave = iter->second.trace_topics;
            size_t i = 0;
            for (auto original_i : workloads_[worker_id] ) {
                trace_topics(original_i) = trace_topics_slave(i);
                i++;
            }
        }
        Count_zh += iter->second.Count_zh;
        Count_sz += iter->second.Count_sz;
        count_h += iter->second.count_h;
        count_z += iter->second.count_z;
        P += iter->second.P;
    }
    P /= n_slaves_;
    kernel_ -> update_state(P);

    DoubleMatrix Theta_zh(count_z.rows(), count_h.rows());
    DoubleMatrix Psi_sz(site2id().size(), count_z.rows());
    aggregate(Count_zh, Count_sz, count_h, count_z,
        hyper_params.alpha_zh, hyper_params.beta_zs, Theta_zh, Psi_sz);
    col_normalize(Theta_zh);
    col_normalize(Psi_sz);
    StampLists stamps(Theta_zh.rows());
    size_t nz = Theta_zh.rows();
    for (size_t z = 0; z < nz; z++ ){
        stamps[z].resize(0);
    }

    const DoubleMatrix & Dts_ref = Dts();
    size_t Dts_last_col = Dts_ref.cols() - 1;
    size_t n_trace = trace_topics.rows();
    for(size_t i = 0; i < n_trace; i++) {
        int z = trace_topics(i);
        stamps.at(z).push_back(Dts_ref(i, Dts_last_col));
    }

    OutPutData result;
    result.alpha_zh = hyper_params.alpha_zh;
    result.beta_zs = hyper_params.beta_zs;
    result.n_iter = hyper_params.n_iter;
    result.burn_in = hyper_params.burn_in;

    result.n_topics = nz;

    result.Theta_zh = std::move(Theta_zh);
    result.Psi_sz = std::move(Psi_sz);
    result.hyper2id = std::move(get<HYPER2ID>(input_data));
    result.site2id = std::move(get<SITE2ID>(input_data));
    result.Count_zh = std::move(Count_zh);
    result.Count_sz = std::move(Count_sz);
    result.count_z = std::move(count_z);
    result.count_h = std::move(count_h);
    result.kernel_name = hyper_params.kernel_name;
    return result;
}

void MasterWorker::add_message (MasterMessage message ) {
    std::unique_lock<std::mutex> lock(this->message_add_mutex);
    this->messages_from_slaves.push(message);
    this->condition.notify_one();
}

void MasterWorker::print_with_lock(size_t worker_id, const string & message) {
    std::unique_lock<std::mutex> lock(this->message_add_mutex);
    cout << "Worker[" << worker_id << "] says:" << endl;
    cout << message << endl;
}

MasterMessage MasterWorker::receive_message() {
    std::unique_lock<std::mutex> lock(this->master_mutex); 
    condition.wait(
            lock, [this] { return !this->messages_from_slaves.empty(); }
            );
    MasterMessage message;
    {
        std::unique_lock<std::mutex> lock__(message_add_mutex);
        message = std::move(this->messages_from_slaves.front()); 
        this->messages_from_slaves.pop();
    }
    return message;
}

OutPutData plearn(
    string trace_fpath, size_t n_workers, size_t n_topics,
    size_t n_iter, double alpha_zh, double beta_zs, string kernel_name,
    vector<double> residency_priors) { 

    size_t burn_in = 0;
    bool dynamic = false;
    size_t n_batches = 0;
#ifdef DEBUG
    cout << "Start learning with " << endl << "trace_fpath: " << trace_fpath << endl
        << "n_topics : " << n_topics << endl 
        << "n_iter : " << n_iter << endl 
        << "alpha_zh : " << alpha_zh << endl
        << "beta_zs : " << beta_zs << endl 
        << "residency_priors: [";
    for (auto & v: residency_priors) 
        cout << v << ", ";
    cout << "]" << endl;
#endif

    HyperParams hyper_params( 
            n_topics, n_iter, burn_in, dynamic, n_batches,
            alpha_zh, beta_zs, kernel_name, residency_priors
            ); 

    auto input = initialize_trace(
                trace_fpath, n_topics, n_iter, 0, std::nullopt, std::nullopt
                );
    MasterWorker worker(n_workers, hyper_params,std::move(input));
    worker.create_slaves();
    OutPutData result = worker.do_manage();
    return result;
}


