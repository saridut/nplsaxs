#!/usr/bin/env bash

#Install miniconda
if [ ! -d ~/soft/miniconda3 ]
then
    mkdir -p ~/soft/miniconda3
fi

#Download miniconda3 to ~/soft/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/soft/miniconda3/miniconda.sh

#Run the installer
bash ~/soft/miniconda3/miniconda.sh -b -u -p ~/soft/miniconda3

#Delete the installer
rm ~/soft/miniconda3/miniconda.sh

#Refresh the shell
source ~/miniconda3/bin/activate
source ~/.bashrc

#Conda init
conda init --all
