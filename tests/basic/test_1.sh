#! /bin/sh



G1=1

Da=0.25
kappa=0.2

d0=0.5
d_min=0.05


if [ ! -d "test_1" ]; then
    mkdir "test_1"
    echo "Utworzono katalog test_1"
fi
cd test_1
cp ../config.txt .

echo gamma = $gamma >>config.txt 
echo kappa = $kappa >>config.txt 
echo Da    = $Da    >>config.txt 
echo d0    = $d0    >>config.txt 
echo d_min = $d_min >>config.txt 

../../karst config.txt #>wyjscie.out 2>bledy.out &

cd ..
