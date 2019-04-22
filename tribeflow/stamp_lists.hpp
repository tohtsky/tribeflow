#include<vector>

using namespace std;

using StampLists = vector<vector<double>>;

//struct StampLists: vector<vector<double>> {
//    inline StampLists(size_t num_topics): vector<vector<double>>(num_topics){} 
//
//    inline void append(size_t topic, double value) {
//        this->at(topic).push_back(value);
//    }
//
//    inline double get(size_t topic, size_t idx) {
//        return this->at(topic).at(idx);
//    }
//};
