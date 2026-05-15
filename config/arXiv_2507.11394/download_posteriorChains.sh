#!/usr/bin/env bash

curl --output temp.tar.gz "https://zenodo.org/api/records/15920131/files/posterior_chains.tar.gz/content"
tar -xf temp.tar.gz
rm temp.tar.gz
mv PosteriorChains/Posterior1/posteriorChain.pkl Posterior1/
mv PosteriorChains/Posterior2/posteriorChain.pkl Posterior2/
mv PosteriorChains/Posterior3/posteriorChain.pkl Posterior3/
rm -fr PosteriorChains
