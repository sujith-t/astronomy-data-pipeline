#!/bin/bash

# install all prerequisite tools
sudo dnf update
sudo dnf install python -y
sudo dnf install git -y
sudo yum install python3-pip -y

# clone the repo
cd /var
git clone https://github.com/sujith-t/astronomy-data-pipeline.git

# create a virtual environment and install the requirements
cd astronomy-data-pipeline
python -m venv astro
source astro/bin/activate
pip install -r pkg_dependancy.txt
