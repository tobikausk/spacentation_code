#!/bin/bash
./generate_xis_and_folder
for i in {0..1199}
do
   echo "This is the $i 'th execution"
   ./Generate_raw_states_attractor $i &
done
