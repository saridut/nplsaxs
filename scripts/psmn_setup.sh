#!/usr/bin/env bash

#ssh s92node01
#mkdir -p /scratch/Cascade/babecass
#module use /applis/PSMN/debian13/Cascade/modules/all
#module load foss/2025a
#module load CMake/3.31.3-GCCcore-14.2.0

if [ ! -f ~/.profile ]
then
    touch ~/.profile
fi
#Needed in both ~/.profile and ~/.bashrc
echo 'source /usr/share/lmod/lmod/init/bash' >> ~/.profile
echo 'source /usr/share/lmod/lmod/init/bash' >> ~/.bashrc

if [ ! -f ~/.bash_profile ]
then
    touch ~/.bash_profile
fi

cat >> ~/.bash_profile << EOF
#Source the bashrc
if [ -f ~/.bashrc ]; then
    source ~/.bashrc
fi
EOF
