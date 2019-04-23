#ifndef PLEAN_HPP
#define PLEAN_HPP

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <thread>
#include <queue>
#include <atomic>
#include <tuple>
#include <optional>
#include <random>
#include <Eigen/Eigen> 
#include "stamp_lists.hpp"
#include "dataio.hpp"
#include "kernels/base.hpp"

using namespace Eigen;
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

struct ResultData{
    inline ResultData(
        const IntegerVector & trace_topics,
        const Eigen::MatrixXi & Count_zh,
        const Eigen::MatrixXi & Count_sz,
        const IntegerVector & count_h,
        const IntegerVector & count_z,
        const Eigen::MatrixXd & P):
        trace_topics(trace_topics), Count_zh(Count_zh), Count_sz(Count_sz),
        count_h(count_h), count_z(count_z), P(P) {} 
    inline ResultData() {}
    ResultData(const ResultData &) = default;
    ResultData & operator = (const ResultData &) = default;
    IntegerVector trace_topics;
    Eigen::MatrixXi Count_zh;
    Eigen::MatrixXi Count_sz;
    IntegerVector count_h;
    IntegerVector count_z;
    Eigen::MatrixXd P;
};

// Count_sz and P
using InterThreadData = std::tuple<Eigen::MatrixXi, Eigen::MatrixXd>; 

struct MasterWorker {
    struct Slave {
        public:
            Slave(MasterWorker* parent, size_t id);
            void start_work();
            void set_pair_id(size_t pair_id); 
            void set_data_other_thread(const InterThreadData & pair_data);

            void join();
            void set_message_from_master(SlaveStatus msg);
            Slave(Slave&&) = default;
            Slave& operator = (Slave&&) = default;

            inline InterThreadData exchanged_data() {
                return std::make_tuple(*send_data_count, *send_data_p);
            }
            inline void set_exchange_data(const InterThreadData & other_data) {
                received_data = other_data;
            }

        private:
            MasterWorker* parent_; 
            const size_t my_id ; 
            std::unique_ptr<std::mutex> slave_mutex;
            std::unique_ptr<condition_variable> slave_condition;
            void main_job();
            void learn();
            void send_message_to_master(SlaveStatus status_code);
            void sample();
            SlaveStatus receive_message_from_master();
            bool paired_update(
                Eigen::MatrixXi & Count_sz, 
                Eigen::MatrixXi & Count_sz_others, 
                Eigen::MatrixXd & P_local
            );

            optional<std::thread> working_thread;
            Eigen::MatrixXi* send_data_count;
            Eigen::MatrixXd* send_data_p; 
            optional<InterThreadData> received_data;
            optional<size_t> pair_id_;
            optional<SlaveStatus> message_from_master;
            map<size_t, Eigen::MatrixXi> previous_encounters_s;
        public:
            optional<ResultData> result; 
    };

    MasterWorker(size_t n_slaves, HyperParams hyper_params,
            InputData && input_data);
    ~MasterWorker();
    void create_slaves();
    OutPutData do_manage();
    void add_message(MasterMessage message);
    inline const map<string, int> & hyper2id () const{
        return std::get<HYPER2ID>(input_data);
    }
    inline const map<string, int> & site2id () const{
        return std::get<SITE2ID>(input_data);
    }

    inline const Eigen::MatrixXd Dts () const {
        return std::get<DTS>(input_data);
    }
    inline const Eigen::MatrixXi Trace () const{
        return std::get<TRACE>(input_data);
    }
    inline const Eigen::MatrixXi & Count_zh () const {
        return std::get<COUNT_ZH>(input_data);
    }
    inline const Eigen::MatrixXi & Count_sz () const {
        return std::get<COUNT_SZ>(input_data);
    } 
    inline const vector<size_t> & trace_hyper_ids () const { 
        return std::get<TRACE_HYPER>(this->input_data);
    }

    inline const IntegerVector & trace_topics () const { 
        return std::get<TRACE_TOPIC>(this->input_data);
    }

    size_t n_slaves () const {
        return n_slaves_;
    }

    const vector<vector<int>> & workloads() {
        return workloads_;
    }

    Eigen::MatrixXd P () {
        return kernel_->get_state();
    }

    private: 
    MasterMessage receive_message();
    size_t n_slaves_;
    public:
    const HyperParams hyper_params; 
    private: 
    InputData input_data; 

    vector<Slave> slaves;
    std::mutex master_mutex;
    std::mutex message_add_mutex;
    std::condition_variable condition;
    queue<MasterMessage> messages_from_slaves;
    std::mt19937 gen;
    vector<vector<int>>workloads_;
    std::unique_ptr<KernelBase> kernel_;
    map<size_t, ResultData> slave_results;
};

OutPutData plearn(string trace_fpath, size_t n_workers, size_t n_topics,
    size_t n_iter, double alpha_zh, double beta_zs, string kernel_name,
    vector<double> residency_priors
    );

void interface_function(
    const string & trace_fpath, size_t n_workers, size_t n_topics,
    size_t n_iter, double alpha_zh, double beta_zs, const string & kernel_name,
    const vector<double> & residency_priors
);

#endif
