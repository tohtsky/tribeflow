#include "thread_pool.hpp"
#include<iostream>
#include<string>
#include<vector>
#include<map>
#include<atomic>
#include<tuple>
#include<optional>

using namespace std;
enum class SlaveStatus {
    STARTED,
    FINISHED,
    PAIRME,
    PAIRED,
    LEARN,
    SENDRESULTS,
    STOP
};

// thread id, message
using MasterMessage = std::tuple<int, SlaveStatus>;
using InterThreadData = string;

struct MasterWorker {
    struct Slave {
        public:
            Slave(MasterWorker* parent, size_t id): parent_(parent), my_id(id) {
            }

            void start_work () {
                this->working_thread = std::thread([this] {
                    this->main_job();
                });

            }

            void set_pair(int id) {
                if (this->pair_id) 
                    throw std::runtime_error("pair_id is already set, and tried to overload..");
                this->pair_id = id;
            }

            void join() {
                this->working_thread->join();
            }

            void set_message_from_master(SlaveStatus msg) {
                message_from_master = msg;
            }

        private:
            void main_job () {
                while (true) {
                    auto event = receive_message_from_master();
                    if (event == SlaveStatus::LEARN ) {
                        send_message_to_master(SlaveStatus::STARTED);
                        std::this_thread::sleep_for(std::chrono::seconds(1));
                        send_message_to_master(SlaveStatus::FINISHED);
                    } else if (event == SlaveStatus::SENDRESULTS) {
                    } else if (event == SlaveStatus::STOP) {
                        break;
                    }
                }
            };

            void send_message_to_master (SlaveStatus status_code) {
                this->parent_->add_message(MasterMessage{this->my_id, status_code});
            };


            SlaveStatus receive_message_from_master () {
                for(;;) {
                    if (message_from_master) break; 
                }
                auto result = *message_from_master;
                message_from_master = std::nullopt;
                return result;
            }

            void send_data_to_other(size_t other_id, InterThreadData data) {
                parent_->slaves[other_id].set_data(data);
            }

            void set_data(InterThreadData data) {
                this->received_data = data; 
            }


            MasterWorker* parent_;
            const size_t my_id ; 
            optional<std::thread> working_thread;
            optional<InterThreadData> received_data;
            optional<int> pair_id;
            optional<SlaveStatus> message_from_master;
    };

    MasterWorker (size_t n_slaves): n_slaves(n_slaves) {
    }

    ~MasterWorker () {
        for(auto & slave: slaves) {
            slave.join();
        }
    }

    void create_slaves () {
        for(size_t i = 0; i < n_slaves; i++) {
            slaves.emplace_back(this, i);
        }
        for(size_t i = 0; i < n_slaves; i++) {
            slaves[i].start_work();
            slaves[i].set_message_from_master(SlaveStatus::LEARN);
        }
    }

    void do_manage () {
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
    void add_message (MasterMessage message ) {
        std::unique_lock<std::mutex> lock(this->message_add_mutex);
        this->messages_from_slaves.push(message);
        this->condition.notify_one();
    }
    private:

    MasterMessage receive_message() {
        std::unique_lock<std::mutex> lock(this->master_mutex); 
        condition.wait(
            lock, [this] { return !this->messages_from_slaves.empty(); }
        );
        MasterMessage message = std::move(this->messages_from_slaves.front());
        this->messages_from_slaves.pop();
        return message;
    }

    size_t n_slaves;
    vector<Slave> slaves;
    std::mutex master_mutex;
    std::mutex message_add_mutex;
    std::condition_variable condition;
    queue<MasterMessage> messages_from_slaves;

};

int main () {
    MasterWorker worker(2); 
    worker.create_slaves();
    worker.do_manage();
}
