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
    make all
    cd ../src || exit
    make clean
    make all

elif [ "$system_name" = "Darwin" ]; then

    echo "This is macOS."
    cd algorithms || exit
    make clean
     make all
     cd ../src || exit
     make clean
     make all

#    dir="build"
#
#    if [ ! -d $dir ]; then
#        mkdir $dir
#    fi
#
##    # Configure build system in 'build' directory
##    cmake -S . -B build
#
#    cd build || exit
#
#    if [ "$1" = "clean" ]; then
#        make clean
#    fi
#
#    make all
##    cmake ..
##    if ! make; then
##      exit 1
##    fi

    cp ../bin/karst ../build/.
else
    echo "Unsupported system: $system_name"
fi

