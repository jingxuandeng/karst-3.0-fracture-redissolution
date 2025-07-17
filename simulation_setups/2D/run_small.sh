#! /bin/sh

printf "Preparing the simulation...\n\n"

bash ./build.sh
if ! bash ./build.sh; then
    echo "Problem with compilation."
    exit 1
fi

cd ~/Desktop/KARST/DATA/2D || exit

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
cp ~/Desktop/KARST/karst_3.0/simulation_setups/2D/config_small.txt ./config.txt || exit


printf "Running the simulation...\n\n"

Da=0.5
gamma=1
kappa=1000
dmin=0.001
cut=true


Da=0.1
d0=0.1
los=13
kappa=1
for N in 10 #20  50  100  200 # `seq 1 1 20`
do
  for gamma in  0.5 #0.25 0.3 0.8 0.75 1.25   #0.01 0.1 0.2 0.5 1 2 5 10 100
  do
  (
                param=Da-$Da-d0-$d0-gamma-$gamma-kappa-$kappa-N-$N
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
                  echo N_x = $N
                  echo N_y = $N


                } >> config.txt

               {
                 time ~/Desktop/KARST/karst_3.0/build/karst config.txt  >wyjscie.out 2>bledy.out
                 } 2>czas.tmp  &

             )
done
done


