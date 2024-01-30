#!/usr/env/bash

python3 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -r iEBE-MUSIC/Cluster_supports/Stampede2/pip_requirements.txt
deactivate
