find_package(Threads REQUIRED)
set(CMAKE_CXX_FLAGS "-O3 -Wall -std=c++11 -march=native")

include_directories($ENV{HOME}/local/include ) 
include_directories(/Users/tomoki.otsuki/src/eigen)
add_executable(plearn main.cpp plearn.cpp dataio.cpp learn_body.cpp kernels/base.cpp)
target_link_libraries (plearn Threads::Threads) 
