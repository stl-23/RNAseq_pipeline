#!/usr/bin/bash
ref=$1
prefix=$2
new_file=$prefix.fa
if [[ ! -e "$new_file" ]];then
    awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0"\n":toupper($0) }' $ref > $new_file
fi
####statistics####
{
awk '$0 !~ /^>/{print length}' $new_file |awk '{s+=$1}END{print s}' >genome_size.tmp
}&
{
grep -v '>' $new_file |grep -o 'A' | wc -l > A.tmp
}&
{
grep -v '>' $new_file |grep -o 'T' | wc -l > T.tmp
}&
{
grep -v '>' $new_file |grep -o 'C' | wc -l > C.tmp
}&
{
grep -v '>' $new_file |grep -o 'G' | wc -l > G.tmp
}&
{
grep -v '>' $new_file |grep -o 'N' | wc -l > N.tmp
}
wait
genome_size=`cat genome_size.tmp`
A=`cat A.tmp`
T=`cat T.tmp`
C=`cat C.tmp`
G=`cat G.tmp`
N=`cat N.tmp`
GC=$((${G}+${C}))
atcg=$((${A}+${T}+${C}+${G}))
GC_content=`echo "scale=2; ${GC}*100/${atcg}" |bc`
echo -e "A:${A}\nT:${T}\nC:${C}\nG:${G}\nN:${N}\nGC content(No N's):${GC_content}%\nTotal:${genome_size}" > $prefix.fa.stat.txt
rm genome_size.tmp A.tmp T.tmp C.tmp G.tmp N.tmp
#### build chr.list ####
grep '>' $new_file |sed 's/>//g' |awk '{split($0,a," ");print a[1]}' > $prefix"_chr.list"
####split genome####
mkdir -p ./split_fa/  && cd ./split_fa/
for i in `grep '>' ../$new_file |sed 's/>//g'`
do 
        echo $i > $i.txt && awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' $i.txt ../$new_file > $i.fa && rm $i.txt
done
