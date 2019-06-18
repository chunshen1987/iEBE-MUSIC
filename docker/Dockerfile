# Build from the official docker python base image, based on Debian
FROM python:3.6-stretch

MAINTAINER Chun Shen <chunshen@wayne.edu>

# Install pre-reqs (commented ones are already in base image)
RUN apt-get update && apt-get install -y \
cmake \
doxygen \
g++ \
gfortran \
gsl-bin \
hdf5-tools \
less \
libxpm-dev \
libgsl-dev \
libhdf5-dev \
rsync \
vim \
zlib1g-dev \
&& rm -rf /var/lib/apt/lists/*

# Install pre-reqs python packages
RUN pip install numpy h5py

# Install packages
ARG software_path=/home/

# 3dMCGlauber
RUN cd ${software_path} \
&& git clone --depth=1 https://github.com/chunshen1987/iEBE-MUSIC \
&& cd iEBE-MUSIC/codes && ./get_code_packages.sh \
&& ./compile_code_packages.sh

# Set up a user group
ARG username=iEBE-MUSIC-user
ARG id=1234
RUN groupadd -g ${id} ${username} \
&& useradd --create-home --home-dir /home/${username} -u ${id} -g ${username} ${username}
USER ${username}
ENV HOME /home/${username}
WORKDIR ${HOME}

ENTRYPOINT /bin/bash
