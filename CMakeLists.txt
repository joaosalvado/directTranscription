cmake_minimum_required(VERSION 3.21)
project(directTranscription)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++17")
set(Python3_EXECUTABLE "/usr/bin/python3")

include_directories(${PROJECT_SOURCE_DIR})



find_package(casadi)
find_package( OpenCV REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )
find_package(Matplot++ REQUIRED)



add_executable(directTranscription main_rh_standard.cpp RK4multipleshooting.cpp RK4multipleshooting.h Plotter_dt.h Plotter_dt.cpp LGLms.cpp LGLms.h CGLms.cpp CGLms.h LLG.cpp LLG.h)

target_link_libraries(directTranscription casadi Matplot++::matplot)