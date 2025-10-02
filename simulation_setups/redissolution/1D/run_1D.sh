#! /bin/sh

printf "Preparing the simulation...\n\n"

#bash ./build.sh
#if ! bash ./build.sh; then
#    echo "Problem with compilation."
#    exit 1
#fi

# Navigate to base directory
cd /Users/jingxuandeng/phd/KRG/research/diss_pre_ML/DATA/1D_test/ || exit 1

# Generate timestamped folder name
# Use hourly/minute if you want unique runs, otherwise daily
 current_date_time=$(date +"test_%Y_%m_%d_%H_%M")
#current_date_time=$(date +"test_%Y_%m_%d")

echo "Target directory: $current_date_time"

# Create directory if it does not exist
if [ ! -d "$current_date_time" ]; then
    mkdir -p "$current_date_time" || { echo "Failed to create directory."; exit 1; }
    echo "Directory '$current_date_time' created successfully."
else
    echo "Directory '$current_date_time' already exists."
fi

# Enter the directory
cd "$current_date_time" || { echo "Failed to enter directory."; exit 1; }
echo "Now in $(pwd)"

# Copy config file
cp "/Users/jingxuandeng/phd/KRG/research/diss_pre_ML/fracture/simulation_setups/redissolution/1D/config_small.txt" "./config.txt" || { echo "Failed to copy config file."; exit 1; }

echo "Config file copied as ./config.txt"

printf "Running the simulation...\n\n"

Da=0.02
d0=0.6
C_eq=0
inlet_cut_factor=1

for kappa2 in 0.1 #0.5 1 #0.0001
do
for kappa in 0.75
do
for gamma in 0.05 #0.001  #0.1 0.5 0.7 0.9 0.99 1.0 1.01 1.1 1.5
do
for dmin in  0.0001
do
  (
                #param=Da-$Da-d0-$d0-gamma-$gamma-kappa-$kappa-cut_factor-1-dyn--0.3-dmin-$dmin-kappa2-$kappa2
                param=kappa2-${kappa2}-Da-$Da-d0-$d0-gamma-$gamma-kappa-$kappa-cut_factor-${inlet_cut_factor}-dyn-${C_eq} #folder name matching Mathematica notebook
                printf "Creating variant: %s\n" "$param"
                mkdir $param
                cd    $param || exit
                cp ../config.txt .

                {
                  echo gamma = $gamma
                  echo kappa = $kappa
                  echo kappa2 = $kappa2
                  echo Da    = $Da
                  echo d0    = $d0
                  echo d_min = $dmin
                  echo "C_eq = $C_eq"
                  echo "inlet_cut_factor = $inlet_cut_factor"
                } >> config.txt

                /Users/jingxuandeng/phd/KRG/research/diss_pre_ML/fracture/bin/karst config.txt  >run_output.txt 2>run_errors.txt&
#                /Users/jingxuandeng/phd/KRG/research/diss_pre_ML/fracture/bin/karst config.txt > >(tee run_output.txt) 2> >(tee run_errors.txt >&2) &
             )
done
done
done
done