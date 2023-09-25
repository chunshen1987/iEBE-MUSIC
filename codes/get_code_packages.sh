#!/usr/bin/env bash

# download the code package

# download 3DMCGlauber
rm -fr 3dMCGlauber_code
git clone --depth=1 https://github.com/chunshen1987/3dMCGlauber 3dMCGlauber_code
(cd 3dMCGlauber_code; git checkout 15e1bb43399b34f3a52e993e9d03577cf1511a80)
rm -fr 3dMCGlauber_code/.git

# download IPGlasma
rm -fr ipglasma_code
git clone --depth=1 https://github.com/chunshen1987/ipglasma ipglasma_code
(cd ipglasma_code; git checkout 65cdb1e29ebcfb1868742ffe8ac6d1bff78012a8)
rm -fr ipglasma_code/.git

# download KoMPoST
rm -fr kompost_code
git clone --depth=1 https://github.com/chunshen1987/KoMPoST kompost_code
(cd kompost_code; git checkout 3a8f873bcf8e20ec8522cb851d20ae5e66610085)
rm -fr kompost_code/.git

# download MUSIC
rm -fr MUSIC_code
git clone --depth=1 https://github.com/MUSIC-fluid/MUSIC -b chun_dev MUSIC_code
(cd MUSIC_code; git checkout 062762b8a15b487259571f35517de9283af0a7ef)
rm -fr MUSIC_code/.git

# download iSS particle sampler
rm -fr iSS_code
git clone --depth=1 https://github.com/chunshen1987/iSS -b dev iSS_code
(cd iSS_code; git checkout b19766ec566b278308654d37462cbade9b941f0e)
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
(cd hadronic_afterburner_toolkit_code; git checkout a5e0901b3bd2b575b01630cd7d29948ee94b35e5)
rm -fr hadronic_afterburner_toolkit_code/.git

