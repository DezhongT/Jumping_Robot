cmake_minimum_required(VERSION 3.11)
project(JumpingRobot)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-Wno-deprecated-declarations -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

find_package(Eigen3 3.4 REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(MKL REQUIRED)

include_directories(
        ${EIGEN_INCLUDE_DIRS}
        ${OPENGL_INCLUDE_DIRS}
        ${GLUT_INCLUDE_DIRS}
        ${MKL_INCLUDE_DIR})

link_directories(${MKL_LIBRARY_DIR})

add_executable(simDER
        src/main.cpp
        src/world.cpp
        src/setInput.cpp
        src/timeStepper.cpp
        src/inertialForce.cpp
        src/externalGravityForce.cpp
        src/dampingForce.cpp
        src/elasticPlate.cpp
        src/elasticStretchingForce.cpp
        src/elasticBendingForce.cpp
        src/externalContactForce.cpp
)

target_link_libraries(simDER
        ${OPENGL_LIBRARIES}
        ${GLUT_LIBRARIES}
        ${MKL_LIBRARIES}
        Eigen3::Eigen
        lapack
        gfortran
        pthread
        rt
        m)
