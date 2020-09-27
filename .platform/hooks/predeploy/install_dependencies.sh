#!/bin/bash
#sudo yum install -y epel-release
#sudo yum install openbabel
cd /
if [ ! -f "mgltools_x86_64Linux2_1.5.6.tar.gz" ]; then
    echo "mgltools_x86_64Linux2_1.5.6.tar.gz does not exist."
    wget http://mgltools.scripps.edu/downloads/downloads/tars/releases/REL1.5.6/mgltools_x86_64Linux2_1.5.6.tar.gz
    tar -xvzf mgltools_x86_64Linux2_1.5.6.tar.gz
    cd mgltools_x86_64Linux2_1.5.6
    ./install.sh
    source ./initMGLtools.sh
fi
cd /
if [ ! -f "Miniconda3-latest-Linux-x86_64.sh" ]; then
    echo "Miniconda3-latest-Linux-x86_64.sh does not exist."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    chmod +x "./Miniconda3-latest-Linux-x86_64.sh"
    ./Miniconda3-latest-Linux-x86_64.sh -b -p /conda
fi
echo "exporting conda path"
export PATH=~/miniconda3/bin:$PATH
conda list
if [ ! $? -eq 0 ];
then
    'alias conda="/conda/bin/conda"' >> ~/.bashrc
    source ~/.bashrc
fi
echo "ls root"
ls /
echo "ls /conda"
ls /conda
echo "ls /conda/bin"
ls /conda/bin
echo "bashrc"
cat ~/.bashrc
echo "bash_profile"
cat ~/.bash_profile
conda list
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