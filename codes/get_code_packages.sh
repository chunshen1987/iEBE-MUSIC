#!/usr/bin/env bash

# download the code package

# download 3DMCGlauber
rm -fr 3dMCGlauber_code
git clone --depth=1 https://github.com/chunshen1987/3dMCGlauber -b electriCharge 3dMCGlauber_code
(cd 3dMCGlauber_code; git checkout ed18f599bbfba89cde2e13ee31c7804317460757)
rm -fr 3dMCGlauber_code/.git

# download IPGlasma
rm -fr ipglasma_code
git clone --depth=1 https://github.com/chunshen1987/ipglasma ipglasma_code
(cd ipglasma_code; git checkout 7195a564e131c3a19d6d1c174ce83aa0d92a2bff)
rm -fr ipglasma_code/.git

# download KoMPoST
rm -fr kompost_code
git clone --depth=1 https://github.com/chunshen1987/KoMPoST kompost_code
(cd kompost_code; git checkout 3a8f873bcf8e20ec8522cb851d20ae5e66610085)
rm -fr kompost_code/.git

# download MUSIC
rm -fr MUSIC_code
git clone --depth=1 https://github.com/MUSIC-fluid/MUSIC -b greg_dev MUSIC_code
(cd MUSIC_code; git checkout 13e6c583ba1c73bd96b8e715d138545124cc7a70)
rm -fr MUSIC_code/.git

# download iSS particle sampler
rm -fr iSS_code
git clone --depth=1 https://github.com/chunshen1987/iSS -b 4DEoS iSS_code
(cd iSS_code; git checkout a3fab9b4c232d6fb01a96e95bb1671c296bc57c7)
rm -fr iSS_code/.git

# download photonEmission wrapper
rm -fr photonEmission_hydroInterface_code
git clone --depth=1 https://github.com/chunshen1987/photonEmission_hydroInterface photonEmission_hydroInterface_code
(cd photonEmission_hydroInterface_code; git checkout a8a9f14c98a3e40519d9704090c3b42eba0107be)
rm -fr photonEmission_hydroInterface_code/.git

# download UrQMD afterburner
rm -fr urqmd_code
git clone --depth=1 https://Chunshen1987@bitbucket.org/Chunshen1987/urqmd_afterburner.git urqmd_code
(cd urqmd_code; git checkout 704c886)
rm -fr urqmd_code/.git

# download hadronic afterner
rm -fr hadronic_afterburner_toolkit_code
git clone --depth=5 https://github.com/chunshen1987/hadronic_afterburner_toolkit -b rapQn hadronic_afterburner_toolkit_code
(cd hadronic_afterburner_toolkit_code; git checkout f78e71ee68a5cc12b7c1d856c34cc445410d64eb)
rm -fr hadronic_afterburner_toolkit_code/.git

# dowload 4D NEoS fomr bitbucket https://bitbucket.org/wayne_state_nuclear_theory/neos/src/main/EoS_UrQMD/

# Bitbucket repository details
USERNAME="wayne_state_nuclear_theory"
REPO_SLUG="neos"
BRANCH="main"

# Bitbucket API base url
API_BASE_URL="https://api.bitbucket.org/2.0/"

# NEOS version and considered quantities
EOS_TYPE="urqmd"
FILE_NAME_LIST=("cs" "mub" "muq" "mus" "p" "t")

for name in "${FILE_NAME_LIST[@]}"; do
    # Define file name
    FILE_PATH="EoS_UrQMD/neos4d_${EOS_TYPE}_${name}_b.dat" 
    # Set the name of URL
    API_FILE_URL="${API_BASE_URL}repositories/${USERNAME}/${REPO_SLUG}/src/${BRANCH}/${FILE_PATH}"

    LOCAL_DESTINATION="neos4d_${EOS_TYPE}_${name}_b.dat"
    curl -L -o "$LOCAL_DESTINATION" "$API_FILE_URL"
done
