find_package(Threads REQUIRED)
set(CMAKE_CXX_FLAGS "-g -Wall -std=c++17")

include_directories($ENV{HOME}/local/include ) 
include_directories(/Users/tomoki.otsuki/src/eigen)
add_executable(plearn main.cpp plearn.cpp dataio.cpp)
target_link_libraries (plearn Threads::Threads) 
