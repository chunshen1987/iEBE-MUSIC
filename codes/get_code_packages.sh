#!/usr/bin/env bash

# download the code package

# download 3DMCGlauber
rm -fr 3dMCGlauber_code
git clone --depth=1 https://github.com/chunshen1987/3dMCGlauber -b electriCharge 3dMCGlauber_code
(cd 3dMCGlauber_code; git checkout 931ab6bc0fff76e6dc9bb492aa5d45f8b3514fdd)
rm -fr 3dMCGlauber_code/.git

# download IPGlasma
rm -fr ipglasma_code
git clone --depth=1 https://github.com/chunshen1987/ipglasma ipglasma_code
(cd ipglasma_code; git checkout 202773b1260ac883e32f802a0aa07cc6340a4ca7)
rm -fr ipglasma_code/.git

# download KoMPoST
rm -fr kompost_code
git clone --depth=1 https://github.com/chunshen1987/KoMPoST kompost_code
(cd kompost_code; git checkout 3a8f873bcf8e20ec8522cb851d20ae5e66610085)
rm -fr kompost_code/.git

# download MUSIC
rm -fr MUSIC_code
git clone --depth=1 https://github.com/MUSIC-fluid/MUSIC -b greg_dev MUSIC_code
(cd MUSIC_code; git checkout e36929e9ba550d9e307d395852c84db985fa5d4e)
rm -fr MUSIC_code/.git

# download iSS particle sampler
rm -fr iSS_code
git clone --depth=1 https://github.com/chunshen1987/iSS -b 4DEoS iSS_code
(cd iSS_code; git checkout a3fab9b4c232d6fb01a96e95bb1671c296bc57c7)
rm -fr iSS_code/.git

# download photonEmission wrapper
rm -fr photonEmission_hydroInterface_code
git clone --depth=1 https://github.com/chunshen1987/photonEmission_hydroInterface photonEmission_hydroInterface_code
(cd photonEmission_hydroInterface_code; git checkout 282397c2ca423886d806c755f120ea7b16572e03)
rm -fr photonEmission_hydroInterface_code/.git

# download UrQMD afterburner
rm -fr urqmd_code
git clone --depth=1 https://Chunshen1987@bitbucket.org/Chunshen1987/urqmd_afterburner.git urqmd_code
(cd urqmd_code; git checkout 704c886)
rm -fr urqmd_code/.git

# download hadronic afterner
rm -fr hadronic_afterburner_toolkit_code
git clone --depth=1 https://github.com/chunshen1987/hadronic_afterburner_toolkit -b rapQn hadronic_afterburner_toolkit_code
(cd hadronic_afterburner_toolkit_code; git checkout fd2c639b8ec76a9e6fc496d3fbd073f93ce4d745)
rm -fr hadronic_afterburner_toolkit_code/.git

