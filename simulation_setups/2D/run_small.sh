#! /bin/sh

printf "Preparing the simulation...\n\n"

bash ./build.sh
if ! bash ./build.sh; then
    echo "Problem with compilation."
    exit 1
fi

cd ../DATA/2D || exit

# Creating proper directory
current_date_time=$(date +small_%Y_%m_%d_%H_%M)
printf $current_date_time
mkdir "$current_date_time"
if [ -d "$current_date_time" ]; then
  echo "Directory '$current_date_time' created successfully."
else
  echo "Failed to create directory."
  exit 1
fi

cd "$current_date_time" || exit
cp ../../../karst_3.0/simulation_setups/2D/config_small.txt ./config.txt || exit


printf "Running the simulation...\n\n"

Da=0.5
gamma=0.9
kappa=1
d0=0.3
dmin=0.0001
cut=true
los=105

for gamma in 1.5 #0.1 0.2 0.5 0.75 1 1.25 2 10
do
  for kappa in  20                                #0.01 0.1 0.2 0.5 1 2 5 10 100
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

                ../../../../karst_3.0/build/karst config.txt # >wyjscie.out 2>bledy.out &

             )
done
done


