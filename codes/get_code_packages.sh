#!/usr/bin/env bash

# download the code package

# download 3DMCGlauber
rm -fr 3dMCGlauber_code
git clone --depth=3 https://github.com/chunshen1987/3dMCGlauber -b electriCharge 3dMCGlauber_code
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
git clone --depth=3 https://github.com/MUSIC-fluid/MUSIC -b greg_dev MUSIC_code
(cd MUSIC_code; git checkout bb7a9a957978b576c77adbe1323d7193e574a49e)
rm -fr MUSIC_code/.git

# download iSS particle sampler
rm -fr iSS_code
git clone --depth=3 https://github.com/chunshen1987/iSS -b 4DEoS iSS_code
(cd iSS_code; git checkout 3b1f8b71200984afbd8b3e86c6754423caa45e41)
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
git clone --depth=3 https://github.com/chunshen1987/hadronic_afterburner_toolkit -b rapQn hadronic_afterburner_toolkit_code
(cd hadronic_afterburner_toolkit_code; git checkout d07edc2741ebeabca45b7504d3ad02f067f4bc11)
rm -fr hadronic_afterburner_toolkit_code/.git

# dowload 4D NEoS fomr bitbucket https://bitbucket.org/wayne_state_nuclear_theory/neos/src/main/EoS_UrQMD/
(cd MUSIC_code/EOS; bash download_Neos4D.sh pdg;)
