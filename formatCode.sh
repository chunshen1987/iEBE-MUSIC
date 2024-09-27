#!/usr/bin/env bash

if command -v yapf &> /dev/null
then
    echo "formatting code ..."
    find . -iname '*.py' -maxdepth 2 | xargs yapf -i
    find analysisKit -iname '*.py' -maxdepth 2 | xargs yapf -i
else
    echo "yapf not found, skip formatting code"
fi
