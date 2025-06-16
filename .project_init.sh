#!/bin/bash

if [[ $PWD == /tmp* ]]; then
    echo "We are /tmp/. This is copier update. Exiting"
    exit 0
fi

if [ -d .git/ ]; then
    echo "This project is already initialized. Exiting"
    exit 0
fi

git init --initial-branch=main
git add .
git commit -a -m '🎉 init: sanger-tol/nf-core-modules'
