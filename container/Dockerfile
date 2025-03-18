# Must be built from the root directory of the repo
# BUILD: docker buildx build -t birtwistlelab/linressim:latest -f container/Dockerfile .
# TEST LOCAL (optional): docker run --rm -it -v <path/to/LinResSims>:/LinResSim birtwistlelab/linressim:latest
# TAG: docker tag linressim birtwistlelab/linressim:latest
# PUSH: docker push birtwistlelab/linressim:latest

FROM ubuntu:22.04

# Copy SPARCED files (ensure build context matches path)
RUN mkdir -p /LinResSims
COPY . /LinResSims/

# Set working directory
WORKDIR /LinResSims

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    SHELL=/bin/bash \
    BLAS_LIBS=-lopenblas \
    LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH \
    PATH=/usr/local/bin:/usr/bin:/bin:$PATH \
    OMPI_ALLOW_RUN_AS_ROOT=1 \
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 

# For reader's reference: 
#          /bin/bash points to the shell executable
#          -lopenblas is a required flag for AMICI to find the OpenBLAS package
#          /usr/lib/x86_64-linux-gnu points AMICI to the C-related solvers
#          python's package manager pip installs binaries at /usr/local/bin
#          python 3 is installed in /usr/bin/
#          OMPI_ALLOW_RUN_AS_ROOT is an OpenMPI flag saying within the contiainer,
#              its alright to run as root user
#          OMPI_ALLOW_RUN_AS_ROOT_CONFIRM makes sure it gets the hint ;)

# Update and install basic dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        curl \
        git \
        python3-pip \
        python3-venv \
        python3-dev \
        rsync \
        openmpi-bin \
        openmpi-common \
        libopenmpi-dev \
        libopenblas-dev \
        liblapack-dev\
        gfortran \
        gcc \
        g++ \
        cmake \
        make \
        libgfortran5 \
        libatomic1 \
        swig \
        libhdf5-dev \
        libsbml-dev \
        && apt-get clean && rm -rf /var/lib/apt/lists/*

# Ubuntu 22.04 doesn't symlink python to python3 by default:
RUN ln -s /usr/bin/python3 /usr/bin/python

# Install Python dependencies using pip
RUN pip install --upgrade pip
RUN pip install dist/linressim-1.0-py3-none-any.whl --verbose --force

# Set default shell
SHELL ["/bin/bash", "-c"]
