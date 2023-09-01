#! /bin/sh

if [ ! -d bin ]; then
    mkdir bin
fi

# Get the system name using uname
system_name=$(uname -s)

# Check if the system is Linux
# shellcheck disable=SC2039
if [ "$system_name" = "Linux" ]; then

    cd algorithms || exit
    make clean
    make
    cd ../src || exit
    make clean
    make

elif [ "$system_name" = "Darwin" ]; then

    echo "This is macOS."
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
else
    echo "Unsupported system: $system_name"
fi

