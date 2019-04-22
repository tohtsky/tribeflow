#ifndef DEBUG_HPP
#define DEBUG_HPP
#include<iostream>
#include<vector>



template<typename T> 
inline void vector_print(vector<T> vals) {
    std::cout << "[";
    for (auto it = vals.begin(); it !=vals.end(); it++ ){
        std::cout << *it << ", ";
    }
    std::cout << "]";
    std::cout << std::endl;
}

#endif
