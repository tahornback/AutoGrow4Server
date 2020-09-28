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
echo "conda before"
conda list
if [ ! $? -eq 0 ];
then
    echo "adding path to bashrc"
    export PATH=/conda/bin:$PATH
    echo "export PATH=/conda/bin:$PATH" >> /root/.bashrc
    echo "alias python=/conda/bin/python" >> /root/.bashrc
    echo "alias conda=/conda/bin/conda" >> /root/.bashrc
    echo "alias pip=/conda/bin/pip" >> /root/.bashrc
    echo "sourcing bashrc"
    source ~/.bashrc
fi
echo "conda after"
conda list
#echo "ls root"
#ls /
#echo "ls /conda"
#ls /conda
#echo "ls /conda/bin"
#ls /conda/bin
echo "bashrc"
cat ~/.bashrc
echo "bash_profile"
cat ~/.bash_profile
conda list
conda install -c conda-forge rdkit
conda install numpy
conda install scipy
/conda/bin/pip install matplotlib==3.2.1
/conda/bin/pip install func_timeout==4.3.5
conda list
#mgltools_directory=$(which mgltools)
#sudo yum install jq
#jsonStr=$(cat ./autogrow4-4.0.2/sample_submit_autogrow.json)
#jq 'del(.mgltools_directory)' <<<"$jsonStr"
#temp='. + { "mgltools_directory": "'
#temp+=mgltools_directory
#temp+='" }'
#jq temp <<<"$jsonStr" > tmp.$$.json && mv tmp.$$.json ./autogrow4-4.0.2/sample_submit_autogrow.json