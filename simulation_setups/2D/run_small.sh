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

Da=0.1
gamma=1
kappa=5
d0=0.1
dmin=0.001

for gamma in  1.5 1 0.75 #0.1 0.5 0.7 0.9 0.99 1.0 1.01 1.1 1.5
do
for kappa in    0.1  1 10 #0.001 0.2 0.3  0.25 0.35
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

                ../../../../karst_3.0/build/karst config.txt # >wyjscie.out 2>bledy.out &

             )
done
done



