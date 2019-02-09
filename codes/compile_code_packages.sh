#!/usr/bin/env bash

# compile MUSIC
echo "compile MUSIC ... "
(
    cd MUSIC_code
    mkdir -p build
    cd build
    cmake ..
    make -j 4
    make install
)
mkdir -p MUSIC
cp MUSIC_code/mpihydro MUSIC/
mkdir -p MUSIC/EOS
(
    cd MUSIC/EOS
    ln -s ../../../codes/MUSIC_code/EOS/hotQCD hotQCD
    ln -s ../../../codes/MUSIC_code/EOS/s95p-v1.2 s95p-v1.2
    ln -s ../../../codes/MUSIC_code/EOS/pdg-urqmd_v3.3+.dat pdg-urqmd_v3.3+.dat
)
cp MUSIC_code/example_inputfiles/IPGlasma_2D/music_input_mode_2 MUSIC/
cp MUSIC_code/utilities/sweeper.sh MUSIC/

# download iSS particle sampler
echo "compile iSS ... "
(
    cd iSS_code
    mkdir -p build
    cd build
    cmake ..
    make -j 4
    make install
)
mkdir -p iSS
cp iSS_code/iSS.e iSS/
cp -r iSS_code/iSS_tables iSS/
cp -r iSS_code/iSS_parameters.dat iSS/

# download UrQMD afterburner
echo "compile UrQMD ... "
(
    cd urqmd_code
    make
)
mkdir -p osc2u
cp urqmd_code/osc2u/osc2u.e osc2u/
mkdir -p urqmd
cp urqmd_code/urqmd/runqmd.sh urqmd/
cp urqmd_code/urqmd/urqmd.e urqmd/
cp urqmd_code/urqmd/uqmd.burner urqmd/


# download hadronic afterner
echo "compile hadronic afterburner toolkit ... "
(
    cd hadronic_afterburner_toolkit_code
    mkdir -p build
    cd build
    cmake ..
    make -j 4
    make install
    cd ../ebe_scripts
    g++ convert_to_binary.cpp -lz -o convert_to_binary.e
    mv convert_to_binary.e ../
)
mkdir -p hadronic_afterburner_toolkit
cp hadronic_afterburner_toolkit_code/hadronic_afterburner_tools.e hadronic_afterburner_toolkit/
cp -r hadronic_afterburner_toolkit_code/EOS hadronic_afterburner_toolkit/
cp hadronic_afterburner_toolkit_code/convert_to_binary.e hadronic_afterburner_toolkit/
cp hadronic_afterburner_toolkit_code/parameters.dat hadronic_afterburner_toolkit/
