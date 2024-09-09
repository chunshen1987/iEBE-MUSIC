#!/usr/bin/env bash

# download the code package

# download 3DMCGlauber
rm -fr 3dMCGlauber_code
#git clone --depth=1 https://github.com/wenbin1501110084/3dMCGlauber 3dMCGlauber_code
git clone --depth=5 https://github.com/chunshen1987/3dMCGlauber -b UPC 3dMCGlauber_code
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
#rm -fr MUSIC_code
#git clone --depth=3 https://github.com/MUSIC-fluid/MUSIC -b greg_dev MUSIC_code
#(cd MUSIC_code; git checkout d513ec4617634f61052ae3d743278fdf87b1d9d2)
#(cd MUSIC_code/EOS; bash download_Neos4D.sh)
#rm -fr MUSIC_code/.git

# download MUSIC
rm -fr MUSIC_code
git clone --depth=1 https://github.com/MUSIC-fluid/MUSIC -b chun_dev MUSIC_code
(cd MUSIC_code; git checkout 062762b8a15b487259571f35517de9283af0a7ef)
(cd MUSIC_code/EOS; bash download_Neos4D.sh)
(cd MUSIC_code/EOS; bash download_Neos2D.sh)
rm -fr MUSIC_code/.git


# download iSS particle sampler
#rm -fr iSS_code
#git clone --depth=3 https://github.com/chunshen1987/iSS -b 4DEoS iSS_code
#(cd iSS_code; git checkout 3b1f8b71200984afbd8b3e86c6754423caa45e41)
#rm -fr iSS_code/.git

# download iSS particle sampler
rm -fr iSS_code
git clone --depth=1 https://github.com/chunshen1987/iSS -b dev iSS_code
(cd iSS_code; git checkout b19766ec566b278308654d37462cbade9b941f0e)
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
git clone --depth=1 https://github.com/chunshen1987/hadronic_afterburner_toolkit --branch UPC hadronic_afterburner_toolkit_code
rm -fr hadronic_afterburner_toolkit_code/.git


