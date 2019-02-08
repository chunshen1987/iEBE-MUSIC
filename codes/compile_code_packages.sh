#!/usr/bin/env bash

# compile MUSIC
echo "compile MUSIC ... "
(
    cd MUSIC
    mkdir build
    cd build
    cmake ..
    make -j 4
    make install
)


# download iSS particle sampler
echo "compile iSS ... "
(
    cd iSS
    mkdir build
    cd build
    cmake ..
    make -j 4
    make install
)

# download UrQMD afterburner
echo "compile UrQMD ... "
(
    cd urqmd_afterburner
    make
)

# download hadronic afterner
echo "compile hadronic afterburner toolkit ... "
(
    cd hadronic_afterburner_toolkit
    mkdir build
    cd build
    cmake ..
    make -j 4
    make install
)
