#define DEBUG
#include "plearn.hpp"

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

void MasterWorker::Slave::set_workload(LearningData data) {
    workload = data;
}

void MasterWorker::Slave::main_job () {
    while (true) {
        auto event = receive_message_from_master();
        if (event == SlaveStatus::LEARN ) {
            send_message_to_master(SlaveStatus::STARTED);
            for(;;) {
                if (workload) {
                    auto workload_body = std::move(*workload);
                    workload = std::nullopt;
                    //sample(std::move(workload_body));
                    break;
                }
            }

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



MasterWorker::MasterWorker (size_t n_slaves): n_slaves(n_slaves) {}

MasterWorker::~MasterWorker () {
    for(auto & slave: slaves) {
        slave.join();
    }
}

void MasterWorker::create_slaves () {
    for(size_t i = 0; i < n_slaves; i++) {
        slaves.emplace_back(this, i);
    }
    LearningData workload = "Debug message";
    for(size_t i = 0; i < n_slaves; i++) {
        slaves[i].start_work();
        slaves[i].set_workload(workload);
        slaves[i].set_message_from_master(SlaveStatus::LEARN);
    }
}

void MasterWorker::do_manage () {
    int available_to_pair = -1;
    unsigned n_finished = 0;
    vector<bool> finished(n_slaves, false);

    while (n_finished != n_slaves) {
        int worker_id;
        SlaveStatus what;
        std::tie(worker_id, what) = receive_message();
        switch (what) {
            case SlaveStatus::STARTED:
                cout << "Worker[" << worker_id << "] has started it's job." << endl;
                break;
            case SlaveStatus::PAIRME:
                if (n_finished == n_slaves -1 ) {
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
                if ( (n_finished == n_slaves -1) && available_to_pair != -1) {
                    slaves[available_to_pair].set_pair(available_to_pair);
                }
                break;
            default:
                throw std::runtime_error(
                        "This message is inappropriate to be received"
                        );
        }
    }
    for (size_t i = 0; i < n_slaves; i++)
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


