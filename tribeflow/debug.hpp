#ifndef DEBUG_HPP
#define DEBUG_HPP
#include<iostream>
#include<sstream>
#include<vector>

template<typename ... Args>
inline void __print(std::ostream& os, Args... args); 

template<typename T, typename ... Args>
inline void __print(std::ostream& os, T t, Args... args) { 
    os << t;
    __print(os, args...);
}

template<>
inline void __print(std::ostream& os) {}

template<typename ... Args>
inline std::string str_concat(Args... args) {
    std::ostringstream os;
    __print(os, args...);
    return os.str();
};


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
