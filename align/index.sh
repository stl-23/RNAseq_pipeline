#!/usr/bin/bash
ref=$1
tool_name=$2
prefix=$3
soft=$4
tool=`grep $tool_name"=" $soft |cut -d '=' -f 2`
new_file=$prefix.aline.fa
if [[ ! -e "$new_file" ]];then
    awk '/^>/&&NR>1{print "";}{ printf "%s",/^>/ ? $0"\n":toupper($0) }' $ref > $new_file
fi
####index####
if [[ $tool_name = 'bwa' ]];then
    genome_size=`awk '$0 !~ /^>/{print length}' $new_file |awk '{s+=$1}END{print s}'`
    if [[ $genome_size -lt 4000000000 ]];then  ## <4G small genome
        $tool index $new_file -p $prefix
    else [[ $genome_size -ge 4000000000 ]]     ## >4G large genome
        $tool index $new_file -a bwtsw -p $prefix
    fi
elif [[ $tool_name = 'samtools' ]]; then
    $tool faidx $new_file
elif [[ $tool_name = 'hisat2' ]]; then
    ${tool}"-build" -p 4 $new_file $prefix
elif [[ $tool_name = 'gatk4' ]];then
    $tool --java-options "-Xmx8G" CreateSequenceDictionary -R $new_file -O $prefix.dict
elif [[ $tool_name = 'fato2bit' ]]; then
    $tool $new_file $prefix.fa.2bit

else
    echo "Error,no such tool in softwares.config or input wrong name, please add software in software.config or use correct name"
fi
