cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(realDataTests)

add_executable(realDataTests main.cpp src/MAMBA_presetting.cpp include/MAMBA_presetting.h src/STRAW_presetting.cpp include/STRAW_presetting.h src/SettingVShape.cpp include/SettingVShape.h src/FittingShape.cpp include/FittingShape.h src/Derivative.cpp include/Derivative.h src/Resolution.cpp include/Resolution.h src/QA.cpp include/QA.h)

set(CMAKE_CXX_STANDARD 14)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT REQUIRED COMPONENTS RIO Net)

include(${ROOT_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

