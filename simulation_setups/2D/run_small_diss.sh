#! /bin/sh

printf "Preparing the simulation...\n\n"

#bash ./build.sh
#if ! bash ./build.sh; then
#    echo "Problem with compilation."
#    exit 1
#fi

cd /Users/jingxuandeng/phd/KRG/research/diss_pre_ML/DATA/3Mineral/2D_test || exit

# Creating proper directory
#current_date_time=$(date +small_%Y_%m_%d_%H_%M)
current_date_time=$(date +%Y_%m_%d)

printf $current_date_time
mkdir "$current_date_time"
if [ -d "$current_date_time" ]; then
  echo "Directory '$current_date_time' created successfully."
else
  echo "Failed to create directory."
  exit 1
fi

cd "$current_date_time" || exit
cp /Users/jingxuandeng/phd/KRG/research/diss_pre_ML/fracture_claude/simulation_setups/2D/config_small_diss.txt ./config.txt || exit


printf "Running the simulation...\n\n"

d0=0.1
dmin=0.001
cut=true
los=106
Va1_perc=0.5
G4=1
if_save_vtk=false
d0=0.1

for Da in 0.1 #0.1 0.2 0.5 0.75 1 1.25 2 10
do
for gamma_a1 in 0.1 1 2 #0.1 0.2 0.5 0.75 1 1.25 2 10
do
  for kappa3 in  0.1 #0.3 0.2                             #0.01 0.1 0.2 0.5 1 2 5 10 100
  do
  (
                param=Da-$Da-d0-$d0-perc-$Va1_perc-kappa3-$kappa3-gamma-$gamma_a1
                printf "Creating variant: %s\n" "$param"
                mkdir $param
                cd    $param || exit
                cp ../config.txt .

                {
#                  echo gamma = $gamma
#                  echo kappa = $kappa
                  echo Da    = $Da
                  echo d0    = $d0
                  echo d_min = $dmin
                  echo if_cut_d_min = $cut
                  echo random_seed = $los
                  echo Va1_perc = $Va1_perc
                  echo kappa3 = $kappa3
                  echo G4 = $G4
                  echo if_save_vtk = $if_save_vtk
                  echo gamma_a1 = $gamma_a1
                } >> config.txt

                /Users/jingxuandeng/phd/KRG/research/diss_pre_ML/fracture_claude/bin/karst config.txt  >run_output.txt 2>run_errors.txt&

             )
done
done
done

