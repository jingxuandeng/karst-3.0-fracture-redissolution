#! /bin/sh


dir="build"

if [ ! -d $dir ]; then
    mkdir $dir
fi

cd build || exit

if [ "$1" = "clean" ]; then
    make clean
fi

cmake ../CMakeLists.txt
if ! make; then
  exit 1
fi

cp ./karst ../tests