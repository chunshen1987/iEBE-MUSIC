#!/usr/env/bash

python3 -m venv venv
source venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r iEBE-MUSIC/Cluster_supports/Anvil/pip_requirements.txt
deactivate
