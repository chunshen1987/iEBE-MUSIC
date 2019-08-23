FROM ubuntu:bionic

MAINTAINER Chun Shen <chunshen@wayne.edu>

RUN export DEBIAN_FRONTEND=noninteractive DEBCONF_NONINTERACTIVE_SEEN=true && \
    apt-get update && apt-get install -y --no-install-recommends \
        cmake \
        g++ \
        gcc \
        gfortran \
        git \
        hdf5-tools \
        libfftw3-dev \
        libgsl-dev \
        libhdf5-dev \
        libz-dev \
        make \
        python3 \
        python3-pip \
        time \
        vim \
        wget \
        && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install pre-reqs python packages
RUN pip3 install --upgrade && \
    pip3 install numpy h5py scipy

ENV HOME /home/

# download and compile the iEBE-MUSIC framework
RUN cd ${HOME} \
&& git clone --depth=1 https://github.com/chunshen1987/iEBE-MUSIC \
&& cd iEBE-MUSIC/codes && ./get_code_packages.sh \
&& ./compile_code_packages.sh

ENV PYTHONIOENCODING utf-8

WORKDIR ${HOME}/iEBE-MUSIC

CMD /bin/bash
