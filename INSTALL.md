# 1) 3ML and LiFF installation (only for Linux, not Mac yet):
#    3ML is a multi-wavelength/multi-messenger analysis framework that provides
#    a unified interface to software from many different instruments. It is 
#    available as open-source code under: https://github.com/giacomov/3ML . 
#    LiFF is the HAWC analysis software (compatible with 3ML) that is used to 
#    study HAWC sky maps. It is available as open source code under:
#    https://github.com/rjlauer/aerie-liff .
#    Both frameworks can be installed via conda:

# Install miniconda
wget -q https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O ~/miniconda.sh

# This will install conda in ~/miniconda . If you want to put it somewhere else change accordingly these 3 lines
bash ~/miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"  
source ~/miniconda/bin/activate

# Create a test environment (to avoid playing with the root conda installation)
conda create --name hawc_test -y python=2.7

# Activate environment
source activate hawc_test

# install liff and threeml
conda install -c giacomov -c hawc -y liff threeml boost=1.61

#install mpmath package 
#(for convenient definition of inverse gamma function, used in cutoff-powerlaw integration)
conda install mpmath

#for convenient environment loading, add to your login script:
export PATH="$HOME/software/miniconda/bin:$PATH"
source ~/software/miniconda/bin/activate
source activate hawc_test
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

#test:
ipython
import threeML
# It should say:
#   Configuration read from /PATH/TO/.threeml/threeML_config.yml
#   INFO [LoadProject.cc, load_project:69]: liff loaded


