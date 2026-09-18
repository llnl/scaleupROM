#!/usr/bin/bash
check_result () {
  # $1: Result output of the previous command ($?)
  # $2: Name of the previous command
  if [ $1 -eq 0 ]; then
      echo "$2 succeeded"
  else
      echo "$2 failed"
      exit -1
  fi
}

module load gcc/11.2.1
module load git/2.36.1
module load mvapich2
module load cmake/3.26.3
module load hdf5-parallel/1.14.0

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
INSTALL_HELPER=$SCRIPT_DIR/../install-helper
LIB_DIR=$SCRIPT_DIR/../dependencies
mkdir -p $LIB_DIR

export CFLAGS="-fPIC ${CFLAGS}"
export CPPFLAGS="-fPIC ${CPPFLAGS}"
export CXXFLAGS="-fPIC ${CXXFLAGS}"

# get libROM first
cd $LIB_DIR
git clone https://github.com/LLNL/libROM.git

# scalapack
if [ -f "$LIB_DIR/lib/libscalapack.a" ]; then
    echo "Using $LIB_DIR/lib/libscalapack.a"
else
    tar -zxvf ./libROM/dependencies/scalapack-2.2.0.tar.gz
    cd scalapack-2.2.0
    cp $INSTALL_HELPER/scalapack/CMakeLists.txt ./
    mkdir build && cd build
    cmake .. -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_INSTALL_PREFIX=$LIB_DIR
    check_result $? ScaLAPACK-config
    
    make -j 16
    make install
    check_result $? ScaLAPACK-installation
fi

# hypre
cd $LIB_DIR
if [ ! -d "hypre" ]; then
  wget https://github.com/hypre-space/hypre/archive/v2.28.0.tar.gz
  tar -zxvf v2.28.0.tar.gz
  mv hypre-2.28.0 hypre
  cd hypre/src
  ./configure --disable-fortran
  make -j
  check_result $? hypre-installation
fi

cd $LIB_DIR
if [ ! -d "gslib" ]; then
    wget -O gslib-1.0.7.tar.gz https://github.com/gslib/gslib/archive/v1.0.7.tar.gz
    export MG=YES
    tar -zxvf gslib-1.0.7.tar.gz
    mv gslib-1.0.7 gslib
    cd ./gslib
    make CC=mpicc -j
    check_result $? gslib-installation
fi

cd $LIB_DIR
if [ ! -d "parmetis-4.0.3" ]; then
    wget -O parmetis-4.0.3.tar.gz https://github.com/LLNL/libROM/raw/master/dependencies/parmetis-4.0.3.tar.gz
    tar -zxvf parmetis-4.0.3.tar.gz
    cd ./parmetis-4.0.3
    make config prefix=$LIB_DIR
    make
    make install
    check_result $? parmetis-installation

    cd ./metis
    make config prefix=$LIB_DIR
    make
    make install
    check_result $? metis-installation
fi

##MUMPS
cd $LIB_DIR
if [ ! -d "mumps" ]; then
    git clone https://github.com/scivision/mumps.git
    check_result $? mumps-clone
    cd ./mumps
    git checkout v5.6.2.1
    export tempMKLROOT=$MKLROOT
    unset MKLROOT
    cp $INSTALL_HELPER/mumps/* ./cmake/
    cmake -B build -Dparmetis=YES -Dparallel=YES -DCMAKE_TLS_VERIFY=OFF -DCMAKE_INSTALL_PREFIX=$LIB_DIR -DCMAKE_INCLUDE_PATH=$LIB_DIR/include -DCMAKE_LIBRARY_PATH=$LIB_DIR/lib -DCMAKE_Fortran_COMPILER=mpif90
    check_result $? mumps-config

    cmake --build build
    cmake --install build
    check_result $? mumps-install
    export MKLROOT=$tempMKLROOT
fi

cd $LIB_DIR
if [ ! -d "mfem" ]; then
    git clone https://github.com/mfem/mfem.git
    check_result $? mfem-clone
    cd mfem
    git checkout v4.6
    mkdir -p ./build && cd ./build
    cmake .. -DBUILD_SHARED_LIBS=YES -DMFEM_USE_MPI=YES -DMFEM_USE_GSLIB=YES -DMFEM_USE_METIS=YES -DMFEM_USE_METIS_5=YES -DMFEM_USE_MUMPS=YES -DGSLIB_DIR="$LIB_DIR/gslib/build" -DMUMPS_DIR="$LIB_DIR/mumps/build/local" -DCMAKE_INCLUDE_PATH=$LIB_DIR/include -DCMAKE_LIBRARY_PATH=$LIB_DIR/lib -DScaLAPACK_DIR=$LIB_DIR/scalapack-2.2.0/build
    check_result $? mfem-config

    make -j 16
    check_result $? mfem-build
    ln -s . include && ln -s . lib
    check_result $? mfem-install
fi

cd $LIB_DIR
if [ ! -d "yaml-cpp" ]; then
    git clone https://github.com/jbeder/yaml-cpp.git
    check_result $? yaml-cpp-clone
    mkdir -p ./yaml-cpp/lib && cd ./yaml-cpp/lib
    cmake .. -DYAML_BUILD_SHARED_LIBS=on -DCMAKE_INSTALL_PREFIX=$LIB_DIR
    check_result $? yaml-config
    make
    check_result $? yaml-build
    make install
    check_result $? yaml-install
fi

cd $LIB_DIR
if [ ! -d "googletest" ]; then
    git clone https://github.com/google/googletest
    check_result $? googletest-clone
    cd ./googletest
    git checkout tags/release-1.12.1 -b v1.12.1
    mkdir ./build && cd ./build
    cmake .. -DCMAKE_INSTALL_PREFIX=$LIB_DIR
    check_result $? gtest-config
    make
    check_result $? gtest-build
    make install
    check_result $? gtest-install
fi

cd $LIB_DIR/libROM
if [ ! -f "build/lib/libROM.so" ]; then
    cp $INSTALL_HELPER/libROM/CMakeLists.txt ./lib/
    mkdir ./build && cd ./build
    cmake .. -DCMAKE_TOOLCHAIN_FILE=../cmake/toolchains/simple.cmake -DCMAKE_BUILD_TYPE=${BUILD_TYPE} -DUSE_MFEM=OFF -DCMAKE_INCLUDE_PATH=$LIB_DIR/include -DCMAKE_LIBRARY_PATH=$LIB_DIR/lib -DHDF5_IS_PARALLEL=On
    check_result $? librom-config
    make -j 16
    check_result $? librom-build
fi

cd $SCRIPT_DIR/../
rm -rf build && mkdir build && cd build
cmake .. -DLIB_DIR=$LIB_DIR -DGTest_ROOT=$LIB_DIR -DYAML_DIR=$LIB_DIR/lib64/cmake/yaml-cpp -DCMAKE_LIBRARY_PATH=$LIB_DIR/lib -DCMAKE_INCLUDE_PATH=$LIB_DIR/include
check_result $? scaleuprom-config
make -j 16
check_result $? scaleuprom-build
