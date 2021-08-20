#!/bin/bash


HOSTS=("yossarian" "zaphod" "billy" "walter" "milo" "jennifer" "raskolnikov" "gatsby" "cathy" "winston")

FOLDERNAMES_0=( 10 11 12 13 14 15 16 17 18 19)

#INPUT (Output block must be commented out while running this)

for i in "${!FOLDERNAMES_0[@]}"; do
  ssh ${HOSTS[i]} "pkill screen; cd Desktop/; rm -r Koustav; mkdir Koustav"
  scp -r ${FOLDERNAMES_0[i]}/* ${HOSTS[i]}:~/Desktop/Koustav
  ssh ${HOSTS[i]} "cd Desktop/Koustav; g++ MultiTransformersDP.cpp cluster_dynamics.cpp -fopenmp; screen -S percolation_stuff -d -m ./a.out"
done
