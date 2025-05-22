#!/usr/bin/env bash

# download the code package

# download 3DMCGlauber
rm -fr 3dMCGlauber_code
git clone --depth=5 https://github.com/chunshen1987/3dMCGlauber 3dMCGlauber_code
(cd 3dMCGlauber_code; git checkout 22955a36969429e99a85957eb6b3c231870a2015)
rm -fr 3dMCGlauber_code/.git

# download IPGlasma
rm -fr ipglasma_code
git clone --depth=1 https://github.com/chunshen1987/ipglasma -b ipglasma_jimwlk ipglasma_code
(cd ipglasma_code; git checkout 4aca36710e1b210f522d2961de82c7a4299246c9)
rm -fr ipglasma_code/.git

# download KoMPoST
rm -fr kompost_code
git clone --depth=1 https://github.com/chunshen1987/KoMPoST kompost_code
(cd kompost_code; git checkout ad5fe9d3b26434bb1d5c29820499ef26808b5a47)
rm -fr kompost_code/.git

# download MUSIC
rm -fr MUSIC_code
git clone --depth=3 https://github.com/MUSIC-fluid/MUSIC -b chun_dev MUSIC_code
(cd MUSIC_code; git checkout e898cc826efc0303f2dd688602e6c604c9ef0f91)
rm -fr MUSIC_code/.git

# download iSS particle sampler
rm -fr iSS_code
git clone --depth=3 https://github.com/chunshen1987/iSS -b dev iSS_code
(cd iSS_code; git checkout b612a8e425d3e1dfc2d2b71cd208df6810c783be)
rm -fr iSS_code/.git

# download photonEmission wrapper
rm -fr photonEmission_hydroInterface_code
git clone --depth=1 https://github.com/chunshen1987/photonEmission_hydroInterface photonEmission_hydroInterface_code
(cd photonEmission_hydroInterface_code; git checkout b80fb78c154cc9131162c8205615faffc86d6a49)
rm -fr photonEmission_hydroInterface_code/.git

# download UrQMD afterburner
rm -fr urqmd_code
git clone --depth=1 https://Chunshen1987@bitbucket.org/Chunshen1987/urqmd_afterburner.git urqmd_code
(cd urqmd_code; git checkout 704c886)
rm -fr urqmd_code/.git

# download hadronic afterner
rm -fr hadronic_afterburner_toolkit_code
git clone --depth=5 https://github.com/chunshen1987/hadronic_afterburner_toolkit -b main hadronic_afterburner_toolkit_code
(cd hadronic_afterburner_toolkit_code; git checkout 4fb0852ab6c1c81f0dba52cbbdc4dcc8a0e92131)
rm -fr hadronic_afterburner_toolkit_code/.git

# download nucleus configurations for 3D-Glauber
(cd 3dMCGlauber_code/tables; bash download_nucleusTables.sh;)
# download nucleus configurations for IP-Glasma
(cd ipglasma_code/nucleusConfigurations; bash download_nucleusTables.sh;)
# download essential EOS files for hydro simulations
(cd MUSIC_code/EOS; bash download_hotQCD.sh; bash download_Neos2D.sh bqs;)
