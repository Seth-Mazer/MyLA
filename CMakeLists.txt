cmake_minimum_required(VERSION 4.2)

project(MyLA VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(MyLA STATIC
        src/MyLA.cpp
)

add_library(MyLA::MyLA ALIAS MyLA)

target_include_directories(MyLA
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_features(MyLA PUBLIC cxx_std_20)