cmake_minimum_required(VERSION 4.1)
project(Tesi0)

set(CMAKE_CXX_STANDARD 20)

find_package(glm REQUIRED)

add_executable(Tesi0
        src/main.cpp
        src/Mesh.cpp
        src/Face.cpp
        src/ExperimentConfig.cpp
)
target_include_directories(Tesi0 PRIVATE src include)
target_link_libraries(Tesi0 PRIVATE glm::glm)
