#! /bin/sh

printf "Preparing the simulation...\n\n"

bash ./build.sh
if ! bash ./build.sh; then
    echo "Problem with compilation."
    exit 1
fi

cd ../DATA/k2_dyn || exit

# Creating proper directory
current_date_time=$(date +small_%Y_%m_%d_%H_%M)

mkdir "$current_date_time"
if [ -d "$current_date_time" ]; then
  echo "Directory '$current_date_time' created successfully."
else
  echo "Failed to create directory."
  exit 1
fi

cd "$current_date_time" || exit
cp ../../../karst_3.0/simulation_setups/k2_dyn/config_small.txt ./config.txt || exit


printf "Running the simulation...\n\n"


dmin=0.001
cut=true
los=107

Da=0.5
d0=0.1

for kappa in 1 #0.1 0.
do
  for gamma in   1 1.25   #0.25 0.3 0.8
  do
  (
                param=Da-$Da-d0-$d0-gamma-$gamma-kappa-$kappa
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
                  echo if_cut_d_min = $cut
                  echo random_seed = $los

                } >> config.txt

                ../../../../karst_3.0/build/karst config.txt >wyjscie.out 2>bledy.out&

             )
done
done


