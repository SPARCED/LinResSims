#!/bin/bash

#--apt-package dependencies install script --#
# For Linux/Ubuntu 22.04-24.02 Desktop distros only #

# Add environment variables to .bashrc if not already present
grep -qxF 'export BLAS_LIBS=-lopenblas' ~/.bashrc || echo 'export BLAS_LIBS=-lopenblas' >> ~/.bashrc
grep -qxF 'export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH' ~/.bashrc || echo 'export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH' >> ~/.bashrc
source ~/.bashrc

# Function to check for root access
check_root() {
    if [ "$EUID" -ne 0 ]; then
        echo "Root access is required to install apt packages. Proceeding without installation..."
        return 1
    fi
    return 0
}

# Update and upgrade the system only if root access is available
if command -v apt >/dev/null 2>&1; then
    echo "Checking for root access to install apt packages..."
    if check_root; then
        echo "Updating system and installing necessary dependencies..."
        sudo apt-get update && sudo apt-get upgrade -y

        # Install necessary dependencies
        sudo     # Set shell and environment
    export DEBIAN_FRONTEND=noninteractive
    export SHELL=/bin/bash

    # Export variables for the build process
    export BLAS_LIBS=-lopenblas
    export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH

    # Update and install basic dependencies
    sudo apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        openmpi-bin \
        openmpi-common \
        libopenmpi-dev \
        libssl-dev \
        libmunge-dev \
        libmunge2 \
        munge \
        bash \
        wget \
        python3-pip \
        python3-venv \
        python3-dev \
        libopenblas-dev \
        gcc \
        g++ \
        git \
        gfortran \
        cmake \
        make \
        file \
        libgfortran5 \
        libatomic1 \
        swig \
        libhdf5-dev \
        libsbml-dev \
        && apt-get clean && rm -rf /var/lib/apt/lists/*

    fi
else
    echo "apt command not found. Please ensure you are using a compatible Linux distribution."
fi

ln -s /usr/bin/python3 /usr/bin/python

# Install Python dependencies using pip
pip install --upgrade pip
if pip install dist/linressim-1.0-py3-none-any.whl --verbose --force; then
    echo "LinResSims was successfully installed!"
else
    echo "LinResSims installation failed. Please check your system's compatibility and dependencies."
fi
