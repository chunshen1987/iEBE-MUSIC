#!/usr/bin/env bash

Green='\033[0;32m'
NC='\033[0m'

machine="$(uname -s)"
case "${machine}" in
    Linux*)     number_of_cores=`nproc --all`;;
    Darwin*)    number_of_cores=`sysctl -n hw.ncpu`;;
    *)          number_of_cores=1;;
esac
number_of_cores_to_compile=$(( ${number_of_cores} > 10 ? 10 : ${number_of_cores} ))

# compile 3dMCGlauber
echo -e "${Green}compile 3dMCGlauber ... ${NC}"
(
    cd 3dMCGlauber_code
    ./get_LHAPDF.sh
    mkdir -p build
    cd build
    cmake .. -Dlink_with_lib=OFF
    make -j${number_of_cores_to_compile}
    make install
)
mkdir -p 3dMCGlauber
cp 3dMCGlauber_code/input 3dMCGlauber/

# compile IPGlasma
echo -e "${Green}compile IPGlasma ... ${NC}"
(
    cd ipglasma_code
    ./compile_IPGlasma.sh
)
mkdir -p ipglasma
cp ipglasma_code/input ipglasma/

# compile KoMPoST
echo -e "${Green}compile KoMPoST ... ${NC}"
(
    cd kompost_code
    make
)
mkdir -p kompost
cp kompost_code/setup.ini kompost/

# compile MUSIC
echo -e "${Green}compile MUSIC ... ${NC}"
(
    cd MUSIC_code
    mkdir -p build
    cd build
    cmake .. -Dlink_with_lib=OFF
    make -j${number_of_cores_to_compile}
    make install
)
mkdir -p MUSIC
#cp MUSIC_code/mpihydro MUSIC/
cp MUSIC_code/example_inputfiles/IPGlasma_2D/music_input_mode_2 MUSIC/
cp MUSIC_code/utilities/sweeper.sh MUSIC/
(cd MUSIC; mkdir initial)

# download iSS particle sampler
echo -e "${Green}compile iSS ... ${NC}"
(
    cd iSS_code
    mkdir -p build
    cd build
    cmake .. -Dlink_with_lib=OFF
    make -j${number_of_cores_to_compile}
    make install
)
mkdir -p iSS
#cp iSS_code/iSS.e iSS/
cp -r iSS_code/iSS_parameters.dat iSS/

# download UrQMD afterburner
echo -e "${Green}compile UrQMD ... ${NC}"
(
    cd urqmd_code
    make -j${number_of_cores_to_compile}
)
mkdir -p osc2u
cp urqmd_code/osc2u/osc2u.e osc2u/
mkdir -p urqmd
cp urqmd_code/urqmd/runqmd.sh urqmd/
#cp urqmd_code/urqmd/urqmd.e urqmd/
cp urqmd_code/urqmd/uqmd.burner urqmd/


# download hadronic afterner
echo -e "${Green}compile hadronic afterburner toolkit ... ${NC}"
(
    cd hadronic_afterburner_toolkit_code
    mkdir -p build
    cd build
    cmake .. -Dlink_with_lib=OFF
    make -j${number_of_cores_to_compile}
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
cp hadronic_afterburner_toolkit_code/ebe_scripts/average_event_HBT_correlation_function.py hadronic_afterburner_toolkit/
