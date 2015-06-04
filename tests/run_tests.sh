#!/bin/bash
for d in */; do
  d=`echo $d | head -c -2`
  echo "\\subsubsection{$d}"
  echo '\begin{center}'
  echo '  \begin{tabular}{| c | c | p{34mm} | p{34mm} | p{34mm} |}'
  echo '    \hline'
  echo '    Pattern & Sequence & Expected & Old result & New result\\'
  t=1
  pat=".pat"
  exp=".exp"
  fa=".fa"
  while [ -f "$d/t$t$pat" ]
  do
    pattern=$(sed "1 q;d" ./$d/t$t$pat)
    pattern=$(echo $pattern | sed -e 's/\$/\\$/g')
    pattern=$(echo $pattern | sed -e 's/{/\\{/g')
    pattern=$(echo $pattern | sed -e 's/}/\\}/g')
    pattern=$(echo $pattern | sed -e 's/\^/\\^{}/g')
    sequence=$(sed "2 q;d" ./$d/t$t$fa)
    expected=""
    lines=$(sed -n '$=' ./$d/t$t$exp)
    i=1
    while [[ "$i" -le "$lines" ]]
    do
      line=$(sed "$i q;d" ./$d/t$t$exp)
      if [ "$i" -gt "1" ]
      then
        expected=$expected" \\newline "$line
      else
        expected=$line
      fi
      i=$((i+1))
    done
    $(../../vanilla/scan_for_matches/scan_for_matches ./$d/t$t$pat < ./$d/t$t$fa > ./tmpout)
    outOld=""
    lines=$(sed -n '$=' ./tmpout)
    i=1
    while [[ "$i" -le "$lines" ]]
    do
      line=$(sed "$i q;d" ./tmpout)
      if [ "$i" -gt "1" ]
      then
        outOld=$outOld" \\newline "$line
      else
        outOld=$line
      fi
      i=$((i+1))
    done
    colourOld="green"
    diff=$(diff -w ./tmpout ./$d/t$t$exp)
    if [ "$diff" != "" ]
    then
      colourOld="red"
    fi

    $(../src/scan_for_matches ./$d/t$t$pat < ./$d/t$t$fa > ./tmpout)
    outNew=""
    lines=$(sed -n '$=' ./tmpout)
    i=1
    while [[ "$i" -le "$lines" ]]
    do
      line=$(sed "$i q;d" ./tmpout)
      if [ "$i" -gt "1" ]
      then
        outNew=$outNew" \\newline "$line
      else
        outNew=$line
      fi
      i=$((i+1))
    done
    colourNew="green"
    diff=$(diff -w ./tmpout ./$d/t$t$exp)
    if [ "$diff" != "" ]
    then
      colourNew="red"
    fi
    echo '    \hline'
    echo "    {\\scriptsize $pattern} & {\\scriptsize $sequence} & {\\scriptsize $expected} & {\\scriptsize \\cellcolor{$colourOld!25} $outOld } & {\\scriptsize \\cellcolor{$colourNew!25} $outNew}\\\\"
    t=$((t+1))
  done
  echo '    \hline'
  echo '  \end{tabular}'
  echo '\end{center}'
done
rm ./tmpout
