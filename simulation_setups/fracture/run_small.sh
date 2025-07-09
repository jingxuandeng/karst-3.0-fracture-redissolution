#! /bin/sh

printf "Preparing the simulation...\n\n"

bash ~/Desktop/KARST/karst_3.0/build.sh
if ! bash ~/Desktop/KARST/karst_3.0/build.sh; then
    echo "Problem with compilation."
    exit 1
fi

cd ~/Desktop/KARST/DATA/fracture/50x50 || exit

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

cp ~/Desktop/KARST/karst_3.0/simulation_setups/fracture/config_small.txt ./config.txt || exit


printf "Running the simulation...\n\n"

Da=0.5
gamma=0.0000001
kappa=1000
dmin=0.001
los=210

Da=0.1
d0=0.3

dyn=1


sandwich_pores=true
no_max_z=false

C_eq=0.5

inlet_cut_factor=1
Da=0.2
gamma=0.5
d0=0.24
kappa=1
dyn=1
for d0 in 0.24  # 0.2 0.3
do
for inlet_cut_factor in 1 3.5

  do
  (
                param=Da-$Da-d0-$d0-gamma-$gamma-kappa-$kappa-cut_factor-$inlet_cut_factor-dyn-$dyn
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
                  echo sandwich_pores = $sandwich_pores
                  echo no_max_z = $no_max_z
                  echo los = $los
                  echo if_dynamic_k2 = true
                  echo dyn_k2_c0 = 1
                  echo inlet_cut_factor = $inlet_cut_factor
#                  echo K_f0 = $K_f0
#                  echo K_f1 = $K_f1
                  echo K_goal = $K_goal




                } >> config.txt

                ~/Desktop/KARST/karst_3.0/build/karst config.txt   >wyjscie.out 2>bledy.out &

             )
done
done




