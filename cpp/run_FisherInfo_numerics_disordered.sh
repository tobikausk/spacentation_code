#!/bin/bash
./generate_seeds_and_folder
for i in {0..99}
do
   echo "This is the $i 'th execution"
   ./Fisher_information_attractor $i &
done
