#! /bin/sh

printf "Preparing the simulation...\n\n"

bash ./build.sh
if ! bash ./build.sh; then
    echo "Problem with compilation."
    exit 1
fi

cd ../DATA/well_with_crack/100x100 || exit

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
cp ../../../../karst_3.0/simulation_setups/well_with_crack/config_small.txt ./config.txt || exit


printf "Running the simulation...\n\n"

Da=0.5
gamma=0.0000001
kappa=1000
dmin=0.001
cut=true
los=107
K_f0=10
K_f1=5
K_goal=1

Da=0.01
d0=0.02
for if_dynamic_k2 in "true" # "false"
do
for inlet_cut_factor in 2 #3 4 5
do
for kappa in 100 #0.1 0.
do
  for gamma in   1.01  #1 1.1 1.05  #2 1 1.5   #0.01 0.1 0.2 0.5 1 2 5 10 100
  do
  (
                param=Da-$Da-d0-$d0-gamma-$gamma-kappa-$kappa-cut_factor-$inlet_cut_factor-dyn-$if_dynamic_k2
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
                  echo inlet_cut_factor = $inlet_cut_factor
                  echo if_dynamic_k2 = $if_dynamic_k2
                  echo K_f0 = $K_f0
                  echo K_f1 = $K_f1
                  echo K_goal = $K_goal


                } >> config.txt

                ../../../../../karst_3.0/build/karst config.txt    >wyjscie.out 2>bledy.out &

             )
done
done
done
done

