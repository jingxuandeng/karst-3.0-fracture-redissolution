#! /bin/sh

printf "Preparing the simulation...\n\n"

bash ./build.sh
if ! bash ./build.sh; then
    echo "Problem with compilation."
    exit 1
fi

cd ../DATA/1D || exit

# Creating proper directory
current_date_time=$(date +test_%Y_%m_%d_%H_%M)
printf $current_date_time
mkdir "$current_date_time"
if [ -d "$current_date_time" ]; then
  echo "Directory '$current_date_time' created successfully."
else
  echo "Failed to create directory."
  exit 1
fi

cd "$current_date_time" || exit

pwd
cp ../../../karst_3.0/simulation_setups/1D/config.txt . || exit


printf "Running the simulation...\n\n"

Da=0.02
d0=0.3


for kappa in 1 0.5
do
for gamma in 1  0.5 #0.1 0.5 0.7 0.9 0.99 1.0 1.01 1.1 1.5
do
for dmin in  0.0001
do
  (
                param=Da-$Da-gamma-$gamma-kappa-$kappa-d0-$d0-dmin-$dmin
                printf "Creating variant: %s\n" "$param"
                mkdir $param
                cd    $param || exit
                cp ../config.txt .

                {
                  echo gamma = $gamma
                  echo kappa = $kappa
                  echo Da    = $Da
                  echo d0    = $d0
                  echo d_min = $dmin
                } >> config.txt

                ../../../../karst_3.0/build/karst config.txt >wyjscie.out 2>bledy.out &
             )
done
done
done