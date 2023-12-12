#! /bin/sh

printf "Preparing the simulation...\n\n"
cd ../..
bash ./build.sh
if ! bash ./build.sh; then
    echo "Problem with compilation."
    exit 1
fi
cd tests/basic || exit
#
#if [ ! -d "test_1" ]; then
#    mkdir "test_1"
#    echo "Utworzono katalog test_1"
#fi

# Creating proper directory
current_date_time=$(date +test_%Y_%m_%d_%H_%M)
mkdir "$current_date_time"
if [ -d "$current_date_time" ]; then
  echo "Directory '$current_date_time' created successfully."
else
  echo "Failed to create directory."
  exit 1
fi

cd "$current_date_time" || exit
rm *
cp ../config.txt .

#
#gamma=1
#Da=0.25
#kappa=0.2
#
#d0=0.1
#d_min=0.05
#{
#  echo gamma = $gamma
#  echo kappa = $kappa
#  echo Da    = $Da
#  echo d0    = $d0
#  echo d_min = $d_min
#} >> config.txt

printf "Running the simulation...\n\n"
../../karst config.txt #>wyjscie.out 2>bledy.out &

cd ..

