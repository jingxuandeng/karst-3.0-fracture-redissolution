#! /bin/sh

printf "Preparing the simulation...\n\n"

bash ~/Desktop/KARST/karst_3.0/build.sh
if ! bash ~/Desktop/KARST/karst_3.0/build.sh; then
    echo "Problem with compilation."
    exit 1
fi

cd ~/Desktop/KARST/DATA/fracture/1D || exit

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
#mkdir debuging_tmp
#cd debuging_tmp || exit

cp ~/Desktop/KARST/karst_3.0/simulation_setups/fracture/config_1D.txt ./config.txt || exit


printf "Running the simulation...\n\n"

Da=0.5
gamma=0.0000001
kappa=1000
dmin=0.001
los=210
#K_f0=1
#K_f1=50
#K_goal=1

Da=0.02
d0=0.3
C_eq=0.2
dyn=1
for C_eq in 0 0.3 # 0.1 0 0.3
do
for d0 in 0.4  # 0.2 0.3
do
for inlet_cut_factor in 1 #5 #3 4 5
do
for kappa in 2 0.5  # 1000  #0.1 0.
do
  for gamma in  1 #1.5  #1 1.1 1.05  #2 1 1.5   #0.01 0.1 0.2 0.5 1 2 5 10 100
  do
  (
                param=Da-$Da-d0-$d0-gamma-$gamma-kappa-$kappa-cut_factor-$inlet_cut_factor-dyn-$C_eq
                printf "Creating variant: %s\n" "$param"
                pwd
                mkdir $param
                cd    $param || exit
                rm *.gz *.pdf
                cp ../config.txt .

                {
                  echo gamma = $gamma
                  echo kappa = $kappa
                  echo Da    = $Da
                  echo d0    = $d0
                  echo random_seed = $los
                  echo inlet_cut_factor = $inlet_cut_factor
                  echo C_eq = $C_eq
                  echo if_dynamic_k2 = true
                  echo dyn_k2_c0 = 1

                } >> config.txt

                ../../../../../karst_3.0/build/karst config.txt   >wyjscie.out 2>bledy.out &

             )
done
done
done
done
done

