#!/usr/bin/env bash

# compile 3dMCGlauber
echo "compile 3dMCGlauber ... "
(
    cd 3dMCGlauber_code
    ./get_LHAPDF.sh
    mkdir -p build
    cd build
    cmake ..
    make -j 4
    make install
)
mkdir -p 3dMCGlauber
cp 3dMCGlauber_code/input 3dMCGlauber/

# compile IPGlasma
echo "compile IPGlasma ... "
(
    cd ipglasma_code
    ./compile_IPGlasma.sh
)
mkdir -p ipglasma
cp ipglasma_code/input ipglasma/

# compile KoMPoST
echo "compile KoMPoST ... "
(
    cd kompost_code
    make
)
mkdir -p kompost
cp kompost_code/setup.ini kompost/

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
#cp MUSIC_code/mpihydro MUSIC/
cp MUSIC_code/example_inputfiles/IPGlasma_2D/music_input_mode_2 MUSIC/
cp MUSIC_code/utilities/sweeper.sh MUSIC/
(cd MUSIC; mkdir initial)

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
#cp iSS_code/iSS.e iSS/
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
#cp urqmd_code/urqmd/urqmd.e urqmd/
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
    g++ concatenate_binary_files.cpp -lz -o concatenate_binary_files.e
    mv concatenate_binary_files.e ../
)
mkdir -p hadronic_afterburner_toolkit
#cp hadronic_afterburner_toolkit_code/hadronic_afterburner_tools.e hadronic_afterburner_toolkit/
#cp -r hadronic_afterburner_toolkit_code/EOS hadronic_afterburner_toolkit/
cp hadronic_afterburner_toolkit_code/convert_to_binary.e hadronic_afterburner_toolkit/
cp hadronic_afterburner_toolkit_code/concatenate_binary_files.e hadronic_afterburner_toolkit/
cp hadronic_afterburner_toolkit_code/parameters.dat hadronic_afterburner_toolkit/
