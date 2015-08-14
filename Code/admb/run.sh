#!/bin/bash
### note this is for models that have been completed
### ./jjm -nox -ind ${1}.ctl -ainp arc/${1}.par -iprint 100
### note this is for models run for the first time (and normal situations)
./jjm -nox -ind ${1}.ctl -iprint 100 ${2}
cp jjm.par arc/${1}.par
cp jjm.rep arc/${1}.rep
cp jjm.std arc/${1}.std
cp $1.prj arc/${1}.prj
### cp jjm.cor arc/${1}.cor
cp Fprof.yld arc/${1}.yld
count=0
countnew=0
for filename in *.rep; do 
  if [ ${filename} == *For_R.rep* ] ; then
    count=$((count+1))
  fi
  if [[ ${filename} == *For_R_?.rep* ]] ; then
    countnew=$((countnew+1))
  fi
done
if [ "$count" == "1" ] ; then
  cp For_R.rep arc/${1}_R.rep 
fi
if [ $countnew -gt 0 ]; then
  for (( i=1; i<=countnew; i++)); do
    cp For_R_${i}.rep arc/${1}_${i}_R.rep
  done
fi
sh cleanad.sh
