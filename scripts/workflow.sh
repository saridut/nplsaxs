#!/usr/bin/env bash

#Go to Cascade login node
ssh s92node01

#Does ~/.bash_profile exist?
#No: Create it
#Yes: Check contents. The following should be there.
cat >> ~/.bash_profile << EOF
#Source the bashrc
if [ -f ~/.bashrc ]; then
    source ~/.bashrc
fi
EOF

#Does ~/.profile exist?
#No:
touch ~/.profile
#Yes: Check contents
#Add these lines
echo 'source /usr/share/lmod/lmod/init/bash' >> ~/.profile
#As well as to ~/.bashrc
echo 'source /usr/share/lmod/lmod/init/bash' >> ~/.bashrc


#Make scratch directory
mkdir -p /scratch/Cascade/babecass
#Check:
ls -lahv /scratch/Cascade | grep babecass
#Create symbolic link in home directory
ln -s /scratch/Cascade/babecass ${HOME}/scratch_cascade

#LOGOUT & LOGIN again to the login node

#Copy load_modules.sh to the home directory
source load_modules.sh

#Install miniconda by script install_conda.sh

#LOGOUT & LOGIN


source load_modules.sh
which g++
which ccmake

#Create new conda environment
conda create -n saxs python=3.14 pip
conda activate saxs
pip install "nplsaxs @ git+https://github.com/saridut/nplsaxs.git"
conda clean --all

#Run job in scratch
mkdir -p ~/scratch_cascade/test_run
cp js_psmn.sh ~/scratch_cascade/test_run
sinfo

sbatch js_psmn.sh

squeue --me

saxsds collate "outdir/*.npz" --coutfile "comb.npz"
