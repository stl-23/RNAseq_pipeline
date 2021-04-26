#!/usr/bin/env bash
ref=$1
prefix=$2
awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0"\n":$0.upper() }' $ref > $prefix.fa
####statistics####
genome_size=`awk '$0 !~ /^>/{print length}' $prefix.fa |awk '{s+=$1}END{print s}'`
A=`grep -v '>' $prefix.fa |grep -o 'A' | wc -l`
T=`grep -v '>' $prefix.fa |grep -o 'T' | wc -l`
C=`grep -v '>' $prefix.fa |grep -o 'C' | wc -l`
G=`grep -v '>' $prefix.fa |grep -o 'G' | wc -l`
N=`grep -v '>' $prefix.fa |grep -o 'N' | wc -l`
GC=${G}+${C}
atcg=${A}+${T}+${C}+${G}
GC_content=`echo "scale=2; ${GC}*100/${atcg}" |bc`
echo -e "A:${A}\nT:${T}\nC:${C}\nG:${G}\nGC content(No N's):${GC_content}%\nTotal:${genome_size}" > $prefix.fa.stat.txt
#### build chr.list ####
grep '>' $prefix.fa |sed 's/>//g' |awk '{split($0,a," "}' > $prefix"_chr.list"
####split genome####
mkdir -p split_fa/  && cd split_fa/
for i in `grep '>' ../$prefix.fa |sed 's/>//g'`;do echo $i > $i.txt && awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' $i.txt ../$prefix.fa > $i.fa && rm *.txt;done




    



