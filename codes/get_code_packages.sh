#!/usr/bin/env bash

# download the code packages at the pinned commits; stop at the first error
set -e

# get_pinned <repository> <commit> <folder>: download exactly <commit>. It is
# fetched by its hash, so it is found however far the branch has moved on
# (a shallow clone of the branch loses it after a few new commits).
get_pinned() {
    rm -fr "$3"
    git init -q "$3"
    git -C "$3" fetch -q --depth=1 "$1" "$2"
    git -C "$3" checkout -q FETCH_HEAD
    echo "$3: $(git -C "$3" log -1 --format='%h %s')"
    rm -fr "$3/.git"
}

# download 3DMCGlauber
get_pinned https://github.com/chunshen1987/3dMCGlauber \
    57f25a066bb37cb5b29ac3d912f33fb6a5085ed0 3dMCGlauber_code

# download IPGlasma (devel branch)
get_pinned https://github.com/schenke/ipglasma \
    5b1a4c6ec8149be5c1992a9c1169da9743a3fe28 ipglasma_code

# download KoMPoST
get_pinned https://github.com/chunshen1987/KoMPoST \
    ad5fe9d3b26434bb1d5c29820499ef26808b5a47 kompost_code

# download MUSIC (main branch)
get_pinned https://github.com/MUSIC-fluid/MUSIC \
    febb2ab687dab86a8d2593daad9760560a27376b MUSIC_code

# download iSS particle sampler (dev branch)
get_pinned https://github.com/chunshen1987/iSS \
    7e1a7a698967b7824b06788c82753f43d44431bb iSS_code

# download photonEmission wrapper
get_pinned https://github.com/chunshen1987/photonEmission_hydroInterface \
    aceacef05f610691e80ecfaf2a8a823d0323917b photonEmission_hydroInterface_code

# download UrQMD afterburner
get_pinned https://Chunshen1987@bitbucket.org/Chunshen1987/urqmd_afterburner.git \
    09eeac28b5861d68d166c1f89b2e97e0a0ebfe8f urqmd_code

# download hadronic afterner (main branch)
get_pinned https://github.com/chunshen1987/hadronic_afterburner_toolkit \
    2326ba21a76d7bb533aa1d39a50f3dca298ef9c3 hadronic_afterburner_toolkit_code

# download nucleus configurations for 3D-Glauber
(cd 3dMCGlauber_code/tables; bash download_nucleusTables.sh;)
# download nucleus configurations for IP-Glasma
(cd ipglasma_code/nucleusConfigurations; bash download_nucleusTables.sh;)
# download essential EOS files for hydro simulations
(cd MUSIC_code/EOS; bash download_hotQCD.sh; bash download_Neos2D.sh bqs;)
