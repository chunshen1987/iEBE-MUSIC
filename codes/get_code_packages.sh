#!/usr/bin/env bash

# download the code package

# download 3DMCGlauber
git clone --depth=1 https://github.com/chunshen1987/3dMCGlauber 3dMCGlauber_code
# download MUSIC
git clone --depth=1 https://github.com/MUSIC-fluid/MUSIC MUSIC_code

# download iSS particle sampler
git clone --depth=1 https://github.com/chunshen1987/iSS iSS_code

# download UrQMD afterburner
git clone --depth=1 https://Chunshen1987@bitbucket.org/Chunshen1987/urqmd_afterburner.git urqmd_code

# download hadronic afterner
git clone --depth=1 https://github.com/chunshen1987/hadronic_afterburner_toolkit hadronic_afterburner_toolkit_code

# download neos_bqs
wget https://s3.amazonaws.com/qcdneos/neosBQS-v0.11.tar.gz && \
tar -zxf neosBQS-v0.11.tar.gz -C MUSIC_code/EOS/ && \
mv MUSIC_code/EOS/neosBQS-v0.11 MUSIC_code/EOS/neos_bqs && \
rm neosBQS-v0.11.tar.gz
