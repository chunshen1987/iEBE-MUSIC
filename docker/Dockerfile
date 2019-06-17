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
ARG software_path=/home/Software
RUN mkdir ${software_path}

# 3dMCGlauber
RUN cd ${software_path} \
&& git clone --depth=1 https://github.com/chunshen1987/3dMCGlauber 3dMCGlauber \
&& cd 3dMCGlauber && ./get_LHAPDF.sh \
&& mkdir build && cd build \
&& cmake .. && make -j4 && make install \
&& export PATH=${software_path}/3dMCGlauber:$PATH

# MUSIC
RUN cd ${software_path} \
&& git clone --depth=1 https://github.com/MUSIC-fluid/MUSIC MUSIC \
&& cd MUSIC && mkdir build && cd build \
&& cmake .. && make -j4 && make install \
&& export PATH=${software_path}/MUSIC:$PATH

# iSS
RUN cd ${software_path} \
&& git clone --depth=1 https://github.com/chunshen1987/iSS iSS \
&& cd iSS && mkdir build && cd build \
&& cmake .. && make -j4 && make install \
&& export PATH=${software_path}/iSS:$PATH

# UrQMD afterburner
RUN cd ${software_path} \
&& git clone --depth=1 https://Chunshen1987@bitbucket.org/Chunshen1987/urqmd_afterburner.git urqmd_code \
&& cd urqmd_code && make && cd .. \
&& mkdir -p osc2u && cp urqmd_code/osc2u/osc2u.e osc2u \
&& mkdir -p urqmd && cp urqmd_code/urqmd/runqmd.sh urqmd/ \
&& cp urqmd_code/urqmd/urqmd.e urqmd/ && cp urqmd_code/urqmd/uqmd.burner urqmd/ \
&& export PATH=${software_path}/osc2u:${software_path}/urqmd:$PATH

# analysis package
RUN cd ${software_path} \
&& git clone --depth=1 https://github.com/chunshen1987/hadronic_afterburner_toolkit hadronic_afterburner_toolkit \
&& cd hadronic_afterburner_toolkit && mkdir build && cd build \
&& cmake .. && make -j4 && make install && cd .. \
&& cd ebe_scripts && g++ convert_to_binary.cpp -lz -o convert_to_binary.e \
&& mv convert_to_binary.e ../ \
&& g++ concatenate_binary_files.cpp -lz -o concatenate_binary_files.e \
&& mv concatenate_binary_files.e ../ && cd .. \
&& export PATH=${software_path}/hadronic_afterburner_toolkit:$PATH

# Set up a user group
ARG username=iEBE-MUSIC-user
ARG id=1234
RUN groupadd -g ${id} ${username} \
&& useradd --create-home --home-dir /home/${username} -u ${id} -g ${username} ${username}
USER ${username}
ENV HOME /home/${username}
WORKDIR ${HOME}

ENTRYPOINT /bin/bash
