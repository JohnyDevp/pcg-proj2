export NVARCH=`uname -s`_`uname -m`
export NVCOMPILERS=/opt/nvidia/hpc_sdk
export MANPATH=$MANPATH:$NVCOMPILERS/$NVARCH/24.11/compilers/man
export PATH=$NVCOMPILERS/$NVARCH/24.11/cmake:$NVCOMPILERS/$NVARCH/24.11/compilers/bin:$PATH
