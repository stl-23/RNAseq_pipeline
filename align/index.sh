#!/usr/bin/env bash
$ref=$1
$tool_name=$2
$prefix=$3
soft=../softwares.config
tool=`grep -w -F $tool_name"=" $soft |cut -d '=' -f 2`
awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0"\n":$0.upper() }' $ref > $prefix.fa
####statistics####
genome_size=`awk '$0 !~ /^>/{print length}' $prefix.fa |awk '{s+=$1}END{print s}'`

####index####
if [[ $tool_name = 'BWA' && $genome_size -lt 4000000000 ]] ; then  ## <4G small genome
  $tool index $prefix.fa -p $prefix
elif [[ $tool_name = 'BWA' && $genome_size -ge 4000000000 ]] ; then ## >4G large genome
  $tool index -a bwtsw -p $prefix
elif [[ $tool_name = 'samtools' ]] ;then
  $tool faidx $prefix.fa
elif [[ $tool_name = 'hisat2' ]];then
  $tool"-build" -p 8 $prefix.fa $prefix
elif [[ $tool_name = 'gatk4' ]];then
  $tool CreateSequenceDictionary -R $prefix.fa -O $prefix.dict
elif [[ $tool_name = 'fato2bit' ]];then
  $tool $prefix.fa $prefix.fa.2bit

else
  echo "Error,no such tool in softwares.config or input wrong name, please add software in software.config or use correct name"
fi



    



