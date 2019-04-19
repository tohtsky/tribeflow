#ifndef PLEAN_HPP
#define PLEAN_HPP

#include<iostream>
#include<string>
#include<vector>
#include<map>
#include<thread>
#include<queue>
#include<atomic>
#include<tuple>
#include<optional>
#include<Eigen/Eigen> 

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

using MasterMessage = std::tuple<int, SlaveStatus>;
using InterThreadData = string;
using LearningData = string;




struct MasterWorker {
    struct Slave {
        public:
            Slave(MasterWorker* parent, size_t id);
            void start_work();
            void set_pair(int id);
            void join();
            void set_message_from_master(SlaveStatus msg);
            void set_workload(LearningData data);
        private:
            void main_job();
            void send_message_to_master(SlaveStatus status_code);
            void sample();
            SlaveStatus receive_message_from_master();
            void send_data_to_other(size_t other_id, InterThreadData data); 
            void set_data(InterThreadData data);

            MasterWorker* parent_;
            const size_t my_id ; 
            optional<std::thread> working_thread;
            optional<InterThreadData> received_data;
            optional<int> pair_id;
            optional<SlaveStatus> message_from_master;
            optional<LearningData> workload; 
    };

    MasterWorker(size_t n_slaves);
    ~MasterWorker();
    void create_slaves();
    void do_manage();
    void add_message(MasterMessage message);
    private: 
    MasterMessage receive_message();
    size_t n_slaves;
    vector<Slave> slaves;
    std::mutex master_mutex;
    std::mutex message_add_mutex;
    std::condition_variable condition;
    queue<MasterMessage> messages_from_slaves;

};

#endif
