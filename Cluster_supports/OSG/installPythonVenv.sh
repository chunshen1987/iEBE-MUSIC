#!/usr/env/bash

python3 -m venv venv
source venv/bin/activate
python3 -m pip install -r iEBE-MUSIC/Cluster_supports/OSG/pip_requirements.txt
deactivate
