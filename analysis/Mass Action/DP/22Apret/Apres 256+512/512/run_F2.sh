#!/bin/bash


HOSTS=("yossarian" "zaphod" "billy" "walter" "milo" "jennifer" "raskolnikov" "gatsby" "cathy" "winston")

FOLDERNAMES_0=(0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19)

#OUTPUT (Input block has to be commented out while running this)

for i in "${!FOLDERNAMES_0[@]}"; do
  rm -r ${FOLDERNAMES_0[i]}
  mkdir ${FOLDERNAMES_0[i]}
  scp -r ${HOSTS[i]}:~/Desktop/Koustav/* ${FOLDERNAMES_0[i]}
done
