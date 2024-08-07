# Jumping Robot Inverse Design

## Overview

This project involves an inverse design process for optimizing parameters of a jumping robot to achieve specific target outputs. The script `inverse_design.py` uses a trained model to predict outcomes and optimize design parameters through simulations.

## Prerequisites

- Python 3.x
- Required libraries (install via `requirements.txt`)
- C++ dependencies

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/DezhongT/Jumping_Robot.git
   cd Jumping_Robot

2. Install dependencies:
   ```bash
   pip install -r requirements.txt

## Usage



There are some dependencies required prior to compilation.
Instructions for macOS and Ubuntu are similar (presented below).
For other operating systems you should be able to modify the commands below appropriately.


- **Note**: Some of these packages are installed to the system library for convenience. You may want to install locally to e.g., `~/.local` to avoid conflicts with system libraries. Add the `cmake` flag: `-D CMAKE_INSTALL_PREFIX=~/.local`. Then `sudo` is not required to install. You'll need to ensure subsequent builds know where to find the build libraries.

- X11
  - An X11 (xorg) server is necessary to use the `freeglut` library. This exists already on Linux.
  - **macOS**: This can be installed with MacPorts: `sudo port install xorg-server`. Then log out and back in.
- [Eigen 3.4.0](http://eigen.tuxfamily.org/index.php?title=Main_Page)
  - Eigen is used for various linear algebra operations.
  - **macOS**: You can install this version with MacPorts: `sudo port install eigen3`. Otherwise, build instructions are below.
  - DisMech is built with Eigen version 3.4.0 which can be downloaded [here](https://gitlab.com/libeigen/eigen/-/releases/3.4.0). After downloading the source code, install through cmake as follows.
    ```bash
    cd eigen-3.4.0 && mkdir build && cd build
    cmake ..
    sudo make install
    ```
    
- [Flexible Collision Library (FCL)](https://github.com/flexible-collision-library/fcl)
  - The FCL library is used to perform both broadphase and narrowphase collision detection with each discrete rod represented as a chain of cylinders.
  - FCL depends on both Eigen (instructions above) and [libccd](https://github.com/danfis/libccd). Install [libccd](https://github.com/danfis/libccd) with the following commands, making sure to build shared libraries:
     ```bash
    git clone https://github.com/danfis/libccd
    cd libccd
    cmake -G "Unix Makefiles" -DBUILD_SHARED_LIBS=ON ..
    make -j4
    sudo make install
     ```
  - Next, install FCL from source using the following commands:
    ```bash
    git clone https://github.com/flexible-collision-library/fcl
    cd fcl && mkdir build && cd build
    cmake ..
    make -j4
    sudo make install
    ```
- [SymEngine](https://github.com/symengine/symengine)
  - SymEngine is used for symbolic differentiation and function generation.
  - **macOS**: SymEngine with LLVM can be installed with MacPorts: `sudo port install symengine`.
  - Before installing SymEngine, LLVM is required which can be installed most easily via a package manager:
    - **Ubuntu**: `sudo apt-get install llvm`
    - **macOS**: `sudo port install llvm-15`
  - Afterwards, install SymEngine from source using the following commands:
    ```bash
    git clone https://github.com/symengine/symengine
    cd symengine && mkdir build && cd build
    cmake -D WITH_LLVM=on -D BUILD_BENCHMARKS=off -D BUILD_TESTS=off ..
    make -j4
    sudo make install
    ```
  - **macOS**: You'll need to provide the LLVM root to the build with `-D CMAKE_PREFIX_PATH=/opt/local/libexec/llvm-15` (if installed via MacPorts).
- [Intel oneAPI Math Kernel Library (oneMKL)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html?operatingsystem=linux&distributions=webdownload&options=online)
  - Necessary for access to Pardiso, which is used as a sparse matrix solver.
  - Intel MKL is also used as the BLAS / LAPACK backend for Eigen.
  - **macOS**: Download from [Intel](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html) and use the install script.
  - **Ubuntu**: Follow the below steps.
    ```bash
    cd /tmp
    wget https://registrationcenter-download.intel.com/akdlm/irc_nas/18483/l_onemkl_p_2022.0.2.136.sh

    # This runs an installer, simply follow the instructions.
    sudo sh ./l_onemkl_p_2022.0.2.136.sh
    ```
  - Add one of the following to your .bashrc so that cmake can find the MKL library. Change the directory accordingly if your MKL version is different. 
   Note that older versions require setting `MKLROOT` while newer versions require `MKL_DIR`.
   You can find out which one from the cmake error message.
    ```bash
    export MKLROOT=/opt/intel/oneapi/mkl/2022.0.2   # for older versions
    export MKL_DIR=/opt/intel/oneapi/mkl/2024.2     # for newer versions
    ```

- [OpenGL / GLUT](https://www.opengl.org/)
  - OpenGL / GLUT is used for rendering the knot through a simple graphic.
  - Simply install through apt package manager:
    - **Ubuntu**: `sudo apt-get install libglu1-mesa-dev freeglut3-dev mesa-common-dev`
    - **macOS**: `sudo port install freeglut pkgconfig` (Note: `pkgconfig` is necessary to avoid finding system GLUT instead of `freeglut`.)

- Lapack (*included in MKL*)
