#!/bin/bash


HOSTS=("jane" "ignatius" "harry" "elliot" "ford" "aziraphale" "crowley" "sidhartha" "arcadio" "aureliano")

FOLDERNAMES_0=(0 1 2 3 4 5 6 7 8 9 )

#OUTPUT (Input block has to be commented out while running this)

for i in "${!FOLDERNAMES_0[@]}"; do
  rm -r ${FOLDERNAMES_0[i]}
  mkdir ${FOLDERNAMES_0[i]}
  scp -r ${HOSTS[i]}:~/Desktop/Koustav/* ${FOLDERNAMES_0[i]}
done
