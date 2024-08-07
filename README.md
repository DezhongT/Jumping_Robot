# Jumping Robot Inverse Design

## Overview

This project involves an inverse design process for optimizing parameters of a jumping robot to achieve specific target outputs. The script `inverse_design.py` uses a trained model to predict outcomes and optimize design parameters through simulations.

## Prerequisites

- Ubuntu 18.04 or above
- Python 3.x
- Python libraries (install via `requirements.txt`)
- C++ dependencies

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/DezhongT/Jumping_Robot.git
   cd Jumping_Robot
   ```

2. Install python libraries:
   ```bash
   pip install -r requirements.txt
   ```
   
3. Install C++ dependencies

   - **Note**: Some of these packages are installed to the system library for convenience. You may want to install locally to e.g., `~/.local` to avoid conflicts with system libraries. Add the `cmake` flag: `-D CMAKE_INSTALL_PREFIX=~/.local`. Then `sudo` is not required to install. You'll need to ensure subsequent builds know where to find the build libraries.

- [Eigen 3.4.0](http://eigen.tuxfamily.org/index.php?title=Main_Page)
  - Eigen is used for various linear algebra operations.
  - The project is built with Eigen version 3.4.0 which can be downloaded [here](https://gitlab.com/libeigen/eigen/-/releases/3.4.0). After downloading the source code, install through cmake as follows.
    ```bash
    cd eigen-3.4.0 && mkdir build && cd build
    cmake ..
    sudo make install
    ```
- [SymEngine](https://github.com/symengine/symengine)
  - SymEngine is used for symbolic differentiation and function generation.
  - Before installing SymEngine, LLVM is required which can be installed most easily via a package manager:
    - **Ubuntu**: `sudo apt-get install llvm`
  - Afterwards, install SymEngine from source using the following commands:
    ```bash
    git clone https://github.com/symengine/symengine
    cd symengine && mkdir build && cd build
    cmake -D WITH_LLVM=on -D BUILD_BENCHMARKS=off -D BUILD_TESTS=off ..
    make -j4
    sudo make install
    ```

- [Intel oneAPI Math Kernel Library (oneMKL)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-download.html?operatingsystem=linux&distributions=webdownload&options=online)
  - Necessary for access to Pardiso, which is used as a sparse matrix solver.
  - Intel MKL is also used as the BLAS / LAPACK backend for Eigen.
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

- Lapack (*included in MKL*)


## Usage

1. Configure the simulation engine
   ```bash
   cd simulations
   mkdir build && cd build
   cmake ..
   make -j4
   cd ../..
   ```

2. (Optional) To simulate the jumping robot with customized setting parameters, simply run
   ```bash
   ./simulations/simDER ./simulations/option.txt
   ```
   The parameters are specified in the ```option.txt``` with specifications as follows (we use SI units):
   - ```render (0 or 1) ```- Flag indicating whether OpenGL visualization should be rendered.
   - ```saveData (0 or 1)``` - Flag indicating whether positions should be recorded.
   - ```YoungM``` - Young's modulus.
   - ```totalTime``` - Total simulation time.
   - ```deltaTime``` - Time step size.
   - ```rodRadius``` - Cross-sectional radius of the beam.
   - ```deltaLength``` - Discretized segment length.
   - ```density``` - Material density.
   - ```stol``` - A small number used in solving the linear system.
   - ```tol``` - Force tolerance.
   - 
   - ```helixradius``` - Radius of the helix.
   - ```helixpitch``` - Pitch of the helix.
   - ```density``` - Mass per unit volume.
   
   - ```Poisson``` - Poisson ratio.
   - ```tol``` and ```stol``` - Small numbers used in solving the linear system. Fraction of a percent, e.g. 1.0e-3, is often a good choice.
   - ```maxIter``` - Maximum number of iterations allowed before the solver quits. 
   - ```gVector``` - 3x1 vector specifying acceleration due to gravity.
   - ```viscosity``` - Viscosity for applying damping forces.
   
   
   - ```recordNodes (0 or 1)``` - Flag indicating whether nodal positions will be recorded.
   - ```dataResolution``` - Rate of data recording in seconds. Applies to both ```saveData``` and ```recordNodes```.
   - ```waitTime``` - Initial wait period duration.
   - ```pullTime``` - Duration to pull for (*starts after ```waitTime``` is done*).
   - ```releaseTime``` - Duration to loosen for (*starts after ```waitTime``` + ```pullTime``` is done*).
   - ```pullSpeed``` - Speed at which to pull and/or loosen each end.
   
   - ```colLimit``` - Distance limit for inclusion in contact candidate set (*colLimit must be > delta*).
   - ```delta``` - Distance tolerance for contact.
   - ```kScaler``` - Constant scaling factor for contact stiffness.
   - ```mu``` - Friction coefficient. A value of zero turns friction off.
   - ```nu``` - Slipping tolerance for friction.
   - ```lineSearch (0 or 1)``` - Flag indicating whether line search will be used.
   - ```knotConfig``` - File name for the initial knot configuration. Should be a txt file located in ```knot_configurations``` directory. Note that overhand knot configurations for ```n1, n2, n3, n4``` are provided with a discretization of 301 nodes.
