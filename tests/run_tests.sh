#!/bin/bash
succ=0
total=0
for d in */; do
  d=`echo $d | head -c -2`
  t=1
  pat=".pat"
  exp=".exp"
  fa=".fa"
  while [ -f "$d/t$t$pat" ]
  do
    $(../src/scan_for_matches ./$d/t$t$pat < ./$d/t$t$fa > ./tmpout)
    outNew=""
    lines=$(sed -n '$=' ./tmpout)
    i=1
    diff=$(diff -w ./tmpout ./$d/t$t$exp)
    succ=$((succ+1))
    total=$((total+1))
    if [ "$diff" != "" ]
    then
      succ=$((succ-1))
      echo "Pattern failed: `cat $d/t$t$pat`"
    fi
    t=$((t+1))
  done
done
echo "Tests passed: $succ/$total"
rm ./tmpout
