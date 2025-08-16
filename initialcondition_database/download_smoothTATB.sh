#!/usr/bin/env bash

# AuAu 200 GeV
wget --no-check-certificate 'https://www.dropbox.com/scl/fi/h3wf7u77ffj9w75e6vrq4/MCGlbAuAu200_sigmaNN_gauss_d0.9_withMultFluct.tar.gz?rlkey=wm9tlkj20rjaqbdnh3nj16vr6&dl=0' -O MCGlbAuAu200_sigmaNN_gauss_d0.9_withMultFluct.tar.gz
tar -xf MCGlbAuAu200_sigmaNN_gauss_d0.9_withMultFluct.tar.gz
rm -rf MCGlbAuAu200_sigmaNN_gauss_d0.9_withMultFluct.tar.gz

# PbPb 5020 GeV
wget --no-check-certificate 'https://www.dropbox.com/scl/fi/6hrehbfl2epqsan4ob152/MCGlbPbPb5020_sigmaNN_gauss_d0.9_withMultFluct.tar.gz?rlkey=w0iwzqpyorvdiqu34mkqefe9v&dl=0' -O MCGlbPbPb5020_sigmaNN_gauss_d0.9_withMultFluct.tar.gz
tar -xf MCGlbPbPb5020_sigmaNN_gauss_d0.9_withMultFluct.tar.gz
rm -rf MCGlbPbPb5020_sigmaNN_gauss_d0.9_withMultFluct.tar.gz
