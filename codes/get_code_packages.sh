#!/usr/bin/env bash

# download the code package

# download 3DMCGlauber
rm -fr 3dMCGlauber_code
git clone --depth=5 https://github.com/chunshen1987/3dMCGlauber 3dMCGlauber_code
(cd 3dMCGlauber_code; git checkout 9af251c283e700a2abb9eb844a91da861136920f)
rm -fr 3dMCGlauber_code/.git

# download IPGlasma
rm -fr ipglasma_code
git clone --depth=1 https://github.com/chunshen1987/ipglasma ipglasma_code
(cd ipglasma_code; git checkout 6cab80d2ad30df2c724b991828202ecf2dffe6af)
rm -fr ipglasma_code/.git

# download KoMPoST
rm -fr kompost_code
git clone --depth=1 https://github.com/chunshen1987/KoMPoST kompost_code
(cd kompost_code; git checkout 3a8f873bcf8e20ec8522cb851d20ae5e66610085)
rm -fr kompost_code/.git

# download MUSIC
rm -fr MUSIC_code
git clone --depth=3 https://github.com/MUSIC-fluid/MUSIC -b eosGenerator MUSIC_code
(cd MUSIC_code; git checkout 4c9c970e791646a2217d71bbd838a42e94318c0e)
rm -fr MUSIC_code/.git

# download iSS particle sampler
rm -fr iSS_code
git clone --depth=3 https://github.com/chunshen1987/iSS -b dev iSS_code
(cd iSS_code; git checkout 3006e5fe0c22c9cc94c6aeb0afa865ecc2563771)
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
(cd hadronic_afterburner_toolkit_code; git checkout 4c976c272b53db80246b3f53725496c0037fafef)
rm -fr hadronic_afterburner_toolkit_code/.git

