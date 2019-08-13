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
neos_bqs_name="neosBQS-v0.11"
neos_bqs_destination="MUSIC_code/EOS/"
wget https://s3.amazonaws.com/qcdneos/${neos_bqs_name}.tar.gz && \
tar -zxf ${neos_bqs_name}.tar.gz -C ${neos_bqs_destination} && \
mv ${neos_bqs_destination}/${neos_bqs_name} ${neos_bqs_destination}/neos_bqs && \
rm ${neos_bqs_name}.tar.gz
