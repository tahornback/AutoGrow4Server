#!/bin/bash
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
    source /root/.bashrc
fi
echo "export PYTHONPATH=/conda/bin/python" >> /root/.bashrc
export PYTHONPATH=/conda/bin/python
conda install -c conda-forge rdkit
conda install numpy
conda install scipy
conda install -c anaconda gunicorn
conda install -c anaconda flask
/conda/bin/pip install matplotlib==3.2.1
/conda/bin/pip install func_timeout==4.3.5
