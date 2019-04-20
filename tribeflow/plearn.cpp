#define DEBUG
#include "plearn.hpp"
#include "learn_body.hpp"

constexpr size_t CACHE_SIZE = 1;

MasterWorker::Slave::Slave(MasterWorker* parent, size_t id) :
    parent_(parent), my_id(id) {
    }


void MasterWorker::Slave::start_work () {
    this->working_thread = std::thread(
        [this] {
            this->main_job();
        });
}

void MasterWorker::Slave::set_pair(int id) {
    if (this->pair_id) 
        throw std::runtime_error("pair_id is already set, and tried to overload..");
    this->pair_id = id;
}

void MasterWorker::Slave::join() {
    this->working_thread->join();
}

void MasterWorker::Slave::set_message_from_master(SlaveStatus msg) {
    message_from_master = msg;
}

void MasterWorker::Slave::learn () {
    const Eigen::MatrixXd & Dts_ref = this->parent_->Dts();
    const Eigen::MatrixXi & Trace_ref = this->parent_->Trace();
    const Eigen::VectorXi & trace_hyper_ids = this->parent_->trace_hyper_ids();
    Eigen::VectorXi trace_topics(this->parent_->trace_topics());
    int nz = this->parent_->Count_zh().rows();
    int nh = this->parent_->hyper2id().size();
    int ns = this->parent_->site2id().size();
    int mem_size = this->parent_->Dts().cols();
    if (mem_size <=0) {
        throw std::runtime_error("mem_size <=0 somehow ");
    }
    Eigen::MatrixXi Count_zh = Eigen::MatrixXi::Zero(nz, nh);
    Eigen::MatrixXi Count_sz = Eigen::MatrixXi::Zero(ns, nz);
    Eigen::VectorXi count_h = Eigen::VectorXi::Zero(nh);
    Eigen::VectorXi count_z = Eigen::VectorXi::Zero(nz);
    fast_populate(
        this->parent_->Trace(), this->parent_->trace_hyper_ids(),
        this->parent_->trace_topics(),
        Count_zh, Count_sz, count_h, count_z
    );

    map<size_t, Eigen::MatrixXi> previous_encounters_s;
    for (size_t worker_id = 0; worker_id < parent_->n_slaves(); worker_id++) {
        previous_encounters_s.insert(
            {worker_id, Eigen::MatrixXi::Zero(ns, nz)}
        );
    }
    StampLists stamps(nz);
    for (size_t i = 0; i < Trace_ref.rows(); i++) {
        auto z = trace_topics(i);
        stamps.at(z).push_back(Dts_ref(i, Dts_ref.cols() - 1));
    }

    Eigen::VectorXd aux = Eigen::VectorXd::Zero(nz);

    Eigen::MatrixXi Count_sz_pair = Eigen::MatrixXi::Zero(ns, nz);
    Eigen::MatrixXi Count_sz_others = Eigen::MatrixXi::Zero(ns, nz);
    Eigen::MatrixXi Count_sz_sum = Eigen::MatrixXi::Zero(ns, nz); 

    Eigen::MatrixXd Theta_zh(nz, nh), Psi_sz(ns, nz);
    bool can_pair = true;
    for (size_t i = 0; i < (parent_->hyper_params.n_iter / CACHE_SIZE) ; i++) {
        Count_sz_sum = Count_sz + Count_sz_others;
        count_z = Count_sz_sum.colwise().sum().transpose();
        // em

        // update local counts
        Count_sz = Count_sz_sum - Count_sz_others;
        count_z = Count_sz.colwise().sum().transpose(); 

        // Update expected belief of other processors

    }
}

void MasterWorker::Slave::main_job () {
    while (true) {
        auto event = receive_message_from_master();
        if (event == SlaveStatus::LEARN ) {
            send_message_to_master(SlaveStatus::STARTED);
            learn();
            //if (my_id == 0) {
            //    cout << parent_->Count_zh().transpose() << endl;
            //}
            send_message_to_master(SlaveStatus::FINISHED);
        } else if (event == SlaveStatus::SENDRESULTS) {
        } else if (event == SlaveStatus::STOP) {
            break;
        }
    }
};

void MasterWorker::Slave::send_message_to_master (SlaveStatus status_code) {
    this->parent_->add_message(MasterMessage{this->my_id, status_code});
};

void MasterWorker::Slave::sample() {
    // map<size_t, bool> previous_encounters_s ;
    bool can_pair = true;
}


SlaveStatus MasterWorker::Slave::receive_message_from_master () {
    for(;;) {
        if (message_from_master) break; 
    }
    auto result = *message_from_master;
    message_from_master = std::nullopt;
    return result;
}

void MasterWorker::Slave::send_data_to_other(size_t other_id, InterThreadData data) {
    parent_->slaves[other_id].set_data(data);
}

void MasterWorker::Slave::set_data(InterThreadData data) {
    this->received_data = data; 
}



MasterWorker::MasterWorker (size_t n_slaves, HyperParams hyper_params,
        InputData && input_data, int random_seed):
    n_slaves_(n_slaves), hyper_params(hyper_params),
    input_data(input_data), gen(random_seed) {}

MasterWorker::~MasterWorker () {
    for(auto & slave: slaves) {
        slave.join();
    }
}

void MasterWorker::create_slaves () {
    size_t nh = hyper2id().size();
    size_t ns = site2id().size();
    size_t nz = Count_zh().rows();
    size_t n_paths = Trace().rows();
    vector<int> hyper_ids(nh);
    for (size_t i = 0; i < nh; i++) hyper_ids[i] = i;

    std::shuffle(hyper_ids.begin(), hyper_ids.end(), gen);
    vector<int> hyper_to_slave_id(nh);
    for (size_t i = 0; i < nh; i++) {
        size_t slave_id = i % this->n_slaves_;
        hyper_to_slave_id[hyper_ids[i]] = slave_id;
    }

    Eigen::MatrixXi workloads = Eigen::MatrixXi::Zero(n_slaves_, n_paths);

    for (size_t i = 0; i < n_paths; i++) {
        int h = (trace_hyper_ids())(i);
        size_t slave_id = hyper_to_slave_id[h];
        workloads(slave_id, i) = 1;
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

void MasterWorker::do_manage () {
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
                if (n_finished == n_slaves_ -1 ) {
                    // reserve it!
                    available_to_pair = worker_id;
                } else {
                    // paired!
                    slaves[worker_id].set_pair(available_to_pair);
                    slaves[available_to_pair].set_pair(worker_id);
                    available_to_pair = -1;
                }
                break;
            case SlaveStatus::FINISHED:
                cout << "Worker[" << worker_id <<
                    "] has finished it's iterations." << endl;
                finished[worker_id] = true;
                n_finished++;
                // The last worker who has requested pairing should be paired with itself.
                if ( (n_finished == n_slaves_ -1) && available_to_pair != -1) {
                    slaves[available_to_pair].set_pair(available_to_pair);
                }
                break;
            default:
                throw std::runtime_error(
                        "This message is inappropriate to be received"
                        );
        }
    }
    for (size_t i = 0; i < n_slaves_; i++)
        slaves[i].set_message_from_master(SlaveStatus::STOP);

}

void MasterWorker::add_message (MasterMessage message ) {
    std::unique_lock<std::mutex> lock(this->message_add_mutex);
    this->messages_from_slaves.push(message);
    this->condition.notify_one();
}

MasterMessage MasterWorker::receive_message() {
    std::unique_lock<std::mutex> lock(this->master_mutex); 
    condition.wait(
            lock, [this] { return !this->messages_from_slaves.empty(); }
            );
    MasterMessage message = std::move(this->messages_from_slaves.front());
    this->messages_from_slaves.pop();
    return message;
}


