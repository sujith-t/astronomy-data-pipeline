#!/bin/bash

# if it's a conda environment enable
conda activate astro

python galaxy_profile_build.py

conda deactivate