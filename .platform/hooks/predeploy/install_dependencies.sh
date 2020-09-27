#!/bin/bash
#sudo yum install -y epel-release
#sudo yum install openbabel
cd /
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x "./Miniconda3-latest-Linux-x86_64.sh"
./Miniconda3-latest-Linux-x86_64.sh -b -p /
cd /
wget http://mgltools.scripps.edu/downloads/downloads/tars/releases/REL1.5.6/mgltools_x86_64Linux2_1.5.6.tar.gz
tar -xvzf mgltools_x86_64Linux2_1.5.6.tar.gz
cd mgltools_x86_64Linux2_1.5.6
./install.sh
source ./initMGLtools.sh
conda install -c conda-forge rdkit
conda install numpy
conda install scipy
#mgltools_directory=$(which mgltools)
#sudo yum install jq
#jsonStr=$(cat ./autogrow4-4.0.2/sample_submit_autogrow.json)
#jq 'del(.mgltools_directory)' <<<"$jsonStr"
#temp='. + { "mgltools_directory": "'
#temp+=mgltools_directory
#temp+='" }'
#jq temp <<<"$jsonStr" > tmp.$$.json && mv tmp.$$.json ./autogrow4-4.0.2/sample_submit_autogrow.json