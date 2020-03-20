#!/usr/bin/env bash

# download the code package

# download 3DMCGlauber
git clone --depth=1 https://github.com/chunshen1987/3dMCGlauber -b remnant 3dMCGlauber_code

# download MUSIC
git clone https://Chunshen1987@bitbucket.org/Chunshen1987/music_dev.git -b new_string_format MUSIC_code

# download iSS particle sampler
git clone --depth=1 https://github.com/chunshen1987/iSS iSS_code

# download UrQMD afterburner
git clone --depth=1 https://Chunshen1987@bitbucket.org/Chunshen1987/urqmd_afterburner.git urqmd_code

# download hadronic afterner
git clone --depth=1 https://github.com/chunshen1987/hadronic_afterburner_toolkit hadronic_afterburner_toolkit_code
