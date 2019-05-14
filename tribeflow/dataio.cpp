#include <random>
#include "dataio.hpp"
#include <algorithm>
#include "debug.hpp"


using namespace std;

size_t count_line(string file_path) {
    std::ifstream ifs(file_path); 
    return std::count(std::istreambuf_iterator<char>(ifs), 
            std::istreambuf_iterator<char>(), '\n') + 1;
}

InputData initialize_trace(string trace_fpath, size_t n_topics, size_t num_iter,
        size_t from_, tl::optional<size_t> to, tl::optional<vector<int>> initial_assign,
        int random_seed) {
    std::cout << "start reading the trace" << std::endl;
    std::mt19937 gen(random_seed);
    std::uniform_int_distribution<> dist_topic_assignment(0, n_topics - 1) ;

    map<pair<int, int>, int> count_zh_dict;
    map<pair<int, int>, int> count_sz_dict;
    map<int, int> count_z_dict;
    map<int, int> count_h_dict;
    map<string, int> hyper2id;
    map<string, int> site2id; 
    vector<string> hyper_names, site_names;

    string line;
    ifstream ifs(trace_fpath);
    if (!ifs) {
        throw std::runtime_error("Couldn't read the file ");
    }

    vector<vector<double>> Dts;
    vector<vector<int>> Trace;
    vector<int> trace_hyper_ids_unsorted, trace_topics_vec;

    size_t i = 0;
    size_t l_in_file = 1;
    tl::optional<size_t> mem_size_lazy;

    while (!ifs.eof()) {
        std::getline(ifs, line);
        if (line.empty()) {
            cerr << " Warning: Got empty line at " << l_in_file << endl;
            l_in_file++;
            continue;
        }
        stringstream ls(line);
        string word;

        vector<string> words;
        while (getline(ls, word, '\t'))
            words.push_back(word);

        if (words.size() < 4) {
            cerr << "got following line at " << i << ": " << endl;
            cerr << line << endl;
            throw std::runtime_error("each row of input data must consists of 4 data.");

        }
        if ( (words.size() - 2) % 2 != 0) {
            throw std::runtime_error("each row must contains 2 + 2 * n_path entries."); 
        }
        size_t mem_size = (words.size() -2 ) / 2;
        if (mem_size_lazy) {
            if (mem_size != *mem_size_lazy) {
                cout << "At line " << i <<", mem_size is " << mem_size
                    << ", which is different from previously"
                    "seen value" << *mem_size_lazy << ". " << endl;
                throw runtime_error("inconsistent mem_size among lines.");
            }
        }
        mem_size_lazy = mem_size;
        vector<double> line_dts(mem_size);
        for (size_t j = 0; j < mem_size; j ++ ) {
            line_dts[j] = stod(words[j]);
        }
        Dts.push_back(line_dts);
        string hyper_str = words[mem_size];
        if (hyper2id.find(hyper_str) == hyper2id.end()) {
            hyper2id[hyper_str] = hyper2id.size();
            hyper_names.push_back(hyper_str);
        }

        int z;
        if (!initial_assign) {
            z = dist_topic_assignment(gen);
        } else {
            z = (*initial_assign)[i];
        }

        int h = hyper2id[hyper_str];
        count_zh_dict[std::pair<int, int>{z, h}] += 1;
        count_h_dict[h] += 1;
        vector<int> line_visits_ids;
        for (size_t j = mem_size + 1; j < words.size() ; j++ ) {
            auto site_name = words[j];
            if (site2id.find(site_name) == site2id.end()) {
                site2id[site_name] = site2id.size();
                site_names.push_back(site_name);
            }
            int o = site2id[site_name];
            line_visits_ids.push_back(o);

            count_sz_dict[pair<int, int>{o, z}] += 1;
            count_z_dict[z] += 1;
        }
        //line_visits_ids.push_back(z);
        trace_hyper_ids_unsorted.push_back(h);
        trace_topics_vec.push_back(z);
        Trace.push_back(line_visits_ids);

        i++;
        l_in_file++;
    }
    vector<int> argsort_indices(i); 
    std::iota(argsort_indices.begin(), argsort_indices.end(), 0);

    std::sort(
        argsort_indices.begin(),
        argsort_indices.end(),
        [&Dts](int i, int j) { return Dts[i].back() < Dts[j].back(); }
    );

    //vector_print(argsort_indices) ;
    DoubleMatrix Dts_mat(argsort_indices.size(), *mem_size_lazy);
    IntegerMatrix Trace_mat(argsort_indices.size(), *mem_size_lazy + 1);
    vector<size_t> trace_hyper_ids(trace_hyper_ids_unsorted.size());
    IntegerVector trace_topics(trace_topics_vec.size());

    for (size_t i = 0; i < argsort_indices.size(); i++) {
        auto ind = argsort_indices[i];
        trace_hyper_ids[i] = trace_hyper_ids_unsorted[ind];
        trace_topics(i) = trace_topics_vec.at(ind);
        for (size_t j = 0; j < *mem_size_lazy; j++ ) {
            Dts_mat(i, j) = Dts[ind][j];
        }
        for (size_t j = 0; j < *mem_size_lazy + 1 ; j++) {
            Trace_mat(i, j) = Trace[ind][j];
        }
    }

    size_t nh = hyper2id.size();
    size_t ns = site2id.size(); 
    size_t nz = n_topics;

    StampLists previous_stamps(n_topics);

    const size_t last_ind_dts = Dts_mat.cols() - 1;
    for (size_t i = 0; i < argsort_indices.size(); i++ ) {
        int topic = trace_topics(i); 
        previous_stamps.at(topic).push_back(Dts_mat(i, last_ind_dts));
    }

    IntegerMatrix Count_zh = IntegerMatrix::Zero(nz, nh);
    IntegerMatrix Count_sz = IntegerMatrix::Zero(ns, nz);
    IntegerVector count_h = IntegerVector::Zero(nh);
    IntegerVector count_z = IntegerVector::Zero(nz);

    for (size_t h = 0; h < nh; h++) {
        count_h(h) = count_h_dict[h];
    } 
    for (size_t z = 0; z < nz; z++) {
        count_z(z) = count_z_dict[z];

        for (size_t h = 0; h < nh; h++) {
            Count_zh(z, h) = count_zh_dict[std::pair<int, int>{z, h}];
        }

        for (size_t s = 0; s < ns; s++) {
            Count_sz(s, z) = count_sz_dict[std::pair<int, int>{s, z}];
        }
    }
    if (Count_sz.colwise().sum().transpose() != count_z) {
        throw std::runtime_error("");
    }
    if (Count_zh.colwise().sum().transpose() != count_h) {
        throw std::runtime_error(""); 
    }
    return std::make_tuple(Dts_mat, Trace_mat, trace_hyper_ids,
            trace_topics, previous_stamps, Count_zh, Count_sz,
           hyper2id, site2id, hyper_names, site_names);
}

