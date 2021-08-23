#!/bin/bash


HOSTS=("jane" "ignatius" "harry" "elliot" "ford" "aziraphale" "crowley" "sidhartha" "arcadio")

FOLDERNAMES=(0 1 2 3 4 5 6 7 8) 

#INPUT (Output block must be commented out while running this)

for i in "${!FOLDERNAMES[@]}"; do
  ssh ${HOSTS[i]} "pkill screen; cd Desktop/; rm -r Koustav; mkdir Koustav" 
  scp -r ${FOLDERNAMES[i]}/* ${HOSTS[i]}:~/Desktop/Koustav
  ssh ${HOSTS[i]} "cd Desktop/Koustav; g++ finding_pc_dp.cpp cluster_dynamics.cpp -fopenmp; screen -S percolation_stuff -d -m ./a.out"
done

#OUTPUT (Input block has to be commented out while running this)

for i in "${!FOLDERNAMES[@]}"; do
  rm -r ${FOLDERNAMES[i]}
  mkdir ${FOLDERNAMES[i]}
  scp -r ${HOSTS[i]}:~/Desktop/Koustav/* ${FOLDERNAMES[i]}
done