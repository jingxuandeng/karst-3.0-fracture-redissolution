#! /bin/sh

printf "Preparing the simulation...\n\n"

#bash ./build.sh
#if ! bash ./build.sh; then # change the path for build.sh
#    echo "Problem with compilation."
#    exit 1
#fi

cd /Users/jingxuandeng/phd/KRG/research/diss_pre_ML/DATA/2D_test/100x200 || exit
#cd /Users/jingxuandeng/phd/KRG/research/diss_pre_ML/DATA/2D_test/radial || exit

# Creating proper directory
current_date_time=$(date +small_%Y_%m_%d_%H_%M_%S)

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

cp /Users/jingxuandeng/phd/KRG/research/diss_pre_ML/fracture/simulation_setups/redissolution/2D/config_small.txt ./config.txt || exit
#cp /Users/jingxuandeng/phd/KRG/research/diss_pre_ML/fracture/simulation_setups/redissolution/2D/config_pH_threshold.txt ./config.txt || exit
#cp /Users/jingxuandeng/phd/KRG/research/diss_pre_ML/fracture/simulation_setups/redissolution/2D/config_radial.txt ./config.txt || exit


printf "Running the simulation...\n\n"

dmin=0.001 # was 0.01 (the larger d_min the easier to get clog)
cut=true
los=107
d0=0.24

redissolution=true
if_randomness_in_regular_net=false
gauss_sigma_d=-0.5

for Da in 0.1 #0.1 0.5
do
for kappa2 in 0.1 1 10
do
for kappa in 0.1 10
do
  for gamma in 3
  do
  (
                param=Da-$Da-d0-$d0-gamma-$gamma-kappa-$kappa-cut_factor-$kappa2-dyn-1
                printf "Creating variant: %s\n" "$param"
                pwd
                mkdir $param
                cd    $param || exit
                rm *.gz *.pdf
                cp ../config.txt .

                {
                  echo gamma = $gamma
                   echo kappa = $kappa
                   echo kappa2 = $kappa2
                   echo Da    = $Da
                   echo d0    = $d0
                   echo d_min = $dmin
                   echo if_cut_d_min = $cut
                   echo random_seed = $los
                   echo if_redissolution = $redissolution
                   echo if_randomness_in_regular_net = $if_randomness_in_regular_net
                   echo gauss_sigma_d = $gauss_sigma_d

                } >> config.txt

                /Users/jingxuandeng/phd/KRG/research/diss_pre_ML/fracture/bin/karst config.txt  >run_output.txt 2>run_errors.txt&
#                 /Users/jingxuandeng/phd/KRG/research/diss_pre_ML/fracture/bin/karst config.txt

             )
done
done
done
done

