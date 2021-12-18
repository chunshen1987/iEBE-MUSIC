#!/usr/bin/env bash

# download the code package

# download 3DMCGlauber-UPC
rm -fr 3dMCGlauber_code
git clone --depth=1 https://github.com/chunshen1987/3dMCGlauber.git --branch UPC 3dMCGlauber_code
rm -fr 3dMCGlauber_code/.git

# download IPGlasma
rm -fr ipglasma_code
git clone --depth=1 https://github.com/chunshen1987/ipglasma ipglasma_code
rm -fr ipglasma_code/.git

# download KoMPoST
rm -fr kompost_code
git clone --depth=1 https://github.com/j-f-paquet/kompost kompost_code
rm -fr kompost_code/.git

# download MUSIC
rm -fr MUSIC_code
git clone --depth=1 https://github.com/MUSIC-fluid/MUSIC -b wenbin-dev MUSIC_code
rm -fr MUSIC_code/.git

# download iSS particle sampler
rm -fr iSS_code
git clone --depth=1 https://github.com/chunshen1987/iSS -b momentumSampler iSS_code
rm -fr iSS_code/.git

# download UrQMD afterburner
rm -fr urqmd_code
git clone --depth=1 https://Chunshen1987@bitbucket.org/Chunshen1987/urqmd_afterburner.git urqmd_code
rm -fr urqmd_code/.git


# download hadronic afterner
rm -fr hadronic_afterburner_toolkit_code
git clone --depth=1 https://github.com/chunshen1987/hadronic_afterburner_toolkit --branch UPC hadronic_afterburner_toolkit_code
rm -fr hadronic_afterburner_toolkit_code/.git
 
