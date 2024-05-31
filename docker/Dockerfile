FROM centos:7

MAINTAINER Chun Shen <chunshen@wayne.edu>

RUN yum -y upgrade

RUN yum -y install epel-release centos-release-scl

# osg repo
#RUN yum -y install http://repo.opensciencegrid.org/osg/3.5/osg-3.5-el7-release-latest.rpm

## well rounded basic system to support a wide range of user jobs
#RUN yum -y grouplist \
#    && yum -y groupinstall "Additional Development" \
#                           "Compatibility Libraries" \
#                           "Console Internet Tools" \
#                           "Development Tools" \
#                           "Internet Applications" \
#                           "Networking Tools" \
#                           "Scientific Support"

RUN yum -y install \
        cmake \
        cmake3 \
        git \
        eigen3-devel \
        gcc \
        gcc-c++ \
        gcc-gfortran \
        devtoolset-7-gcc \
        devtoolset-7-gcc-c++ \
        devtoolset-7-gcc-gfortran \
        fftw \
        fftw-devel \
        hdf5-devel \
        zlib-devel \
        make \
        openmpi-devel \
        python3 \
        python3-pip \
        ssh \
        time \
        vim \
        wget \
        && \
    yum clean all

RUN alternatives --install /usr/local/bin/cmake cmake /usr/bin/cmake 10 \
    --slave /usr/local/bin/ctest ctest /usr/bin/ctest \
    --slave /usr/local/bin/cpack cpack /usr/bin/cpack \
    --slave /usr/local/bin/ccmake ccmake /usr/bin/ccmake \
    --family cmake

RUN  alternatives --install /usr/local/bin/cmake cmake /usr/bin/cmake3 20 \
--slave /usr/local/bin/ctest ctest /usr/bin/ctest3 \
--slave /usr/local/bin/cpack cpack /usr/bin/cpack3 \
--slave /usr/local/bin/ccmake ccmake /usr/bin/ccmake3 \
--family cmake

SHELL ["/usr/bin/scl", "enable", "devtoolset-7"]

RUN yum -y install http://springdale.princeton.edu/data/springdale/7/x86_64/os/Computational/gsl25-2.5-3.sdl7.2.x86_64.rpm

RUN yum -y install http://springdale.princeton.edu/data/springdale/7/x86_64/os/Computational/gsl25-devel-2.5-3.sdl7.2.x86_64.rpm

# osg
#RUN yum -y install osg-ca-certs osg-wn-client && \
#        rm -f /etc/grid-security/certificates/*.r0

# Install pre-reqs python packages
RUN pip3 install --upgrade && \
    pip3 install numpy h5py scipy

ENV HOME /opt/
ENV PATH="${PATH}:/usr/lib64/openmpi/bin:/usr/local/gsl/2.5/x86_64/bin"
ENV LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/usr/local/lib:/usr/local/gsl/2.5/x86_64/lib64"
ENV GSL_ROOT_DIR="/usr/local/gsl/2.5/x86_64"

RUN echo 'export PATH="${PATH}:/usr/lib64/openmpi/bin:/usr/local/gsl/2.5/x86_64/bin"' >> /root/.bashrc
RUN echo 'export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/usr/local/lib:/usr/local/gsl/2.5/x86_64/lib64"' >> /root/.bashrc
RUN echo 'export GSL_ROOT_DIR="/usr/local/gsl/2.5/x86_64"' >> /root/.bashrc

# download and compile the iEBE-MUSIC framework
RUN cd ${HOME} && \
    git clone --depth=1 https://github.com/chunshen1987/iEBE-MUSIC -b 4DEoS && \
    cd iEBE-MUSIC/codes && \
    ./get_code_packages.sh && \
    ./compile_code_packages.sh

# Cleaning caches to reduce size of image
RUN yum clean all

# required directories
RUN for MNTPOINT in \
        /cvmfs \
        /hadoop \
        /hdfs \
        /lizard \
        /mnt/hadoop \
        /mnt/hdfs \
        /xenon \
        /scratch \
        /spt \
        /stash2 \
    ; do \
        mkdir -p $MNTPOINT ; \
    done

# some extra singularity stuff
COPY .singularity.d /.singularity.d
RUN cd / && \
    ln -s .singularity.d/actions/exec .exec && \
    ln -s .singularity.d/actions/run .run && \
    ln -s .singularity.d/actions/test .shell && \
    ln -s .singularity.d/runscript singularity

ENV PYTHONIOENCODING utf-8
WORKDIR ${HOME}/

CMD scl enable devtoolset-7 bash
