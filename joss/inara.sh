#!/bin/bash
# https://github.com/openjournals/inara
# Boot "Docker Desktop" before running this script.
# Run `bash inara.sh` on pwsh or WSL.
docker run --rm -it \
    -v $PWD:/data \
    -u $(id -u):$(id -g) \
    openjournals/inara \
    -o pdf,preprint \
    -p \
    ./paper.md