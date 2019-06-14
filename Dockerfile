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
python-h5py \
rsync \
vim \
zlib1g-dev \
&& rm -rf /var/lib/apt/lists/*

# Set up a user group
ARG username=iEBE-MUSIC-user
ARG id=1234
RUN groupadd -g ${id} ${username} \
&& useradd --create-home --home-dir /home/${username} -u ${id} -g ${username} ${username}
USER ${username}
ENV HOME /home/${username}
WORKDIR ${HOME}

ENTRYPOINT /bin/bash
