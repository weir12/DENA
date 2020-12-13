#!/bin/bash
# @Author: liangou
# @Date:   2020-11-23 14:55:45
# @Last Modified by:   liangou
# @Last Modified time: 2020-12-09 15:47:28

#config info
#Plz place the config file in the same directory
source ./config
set -e
#defien functions

CheckFolder(){
	dir=$1
    if [ "`ls -A $dir`" = "" ]; then
    	echo -e "The folder used to save the result is empty,nothing needed to do"
    else
    	while [ 1 ]; do
    	echo -e "Tip: The folder already exists"
    	read -p "Clear folder(c) or auto overwrite(w):" cho
    	case $cho in
    		"c")
			rm -rf ${dir}/* && echo "Clear Done" && break
			;;
			"w")
			echo -e "Folder will be automatically overwritten" && break
			;;
			*)
			echo -e "Please enter c(clean) or w(auto overwrite)"
			;;
		esac
	done
fi
return 0
}

#Extract the packages containing fast5 files and classify them according to the level of m6A
unzip_data(){
	declare -A dic
	dic=(['vir1']="non-m6a" [VIRc]="m6a" )
	for condition in $(echo ${!dic[*]})
	do
		for (( i = 1; i <= 4; i++ ))
		do
			echo $condition
			echo $i
			mkdir -p /home/weir/tair_rawdata/elife_rawdata/${dic[${condition}]}/${i}/
			(tar -zxf ${condition}_nanopore_drs_${i}.tar.gz  -C /home/weir/tair_rawdata/elife_rawdata/${dic[${condition}]}/${i}/ && \
			echo -e `date "+%Y-%m-%d %H:%M:%S"` ${condition}"\t"${i} has successfully done >> unzip.log ) &
		done
	done
}

#quality control of FASTQ files
read_QC(){
	QC_result="/home/weir/m6a_model/result/qc_result/"
	mkdir -p ${QC_result}{"VIRc","vir1"}
	CheckFolder $QC_result
	NanoPlot -t ${threads} --fastq ${vir1_fastq}  --maxlength 40000 --plots hex dot --N50 --drop_outliers -o ${QC_result}vir1 &
	NanoPlot -t ${thread} --fastq ${VIRc_fastq}  --maxlength 40000 --plots hex dot --N50 --drop_outliers -o ${QC_result}VIRc &
}


#mapping reads to ref(tair10) and sort+index bam files
mapping(){
	align_result="${output_wd}/align_result/"
	CheckFolder $align_result
	mkdir -p ${align_result}{${neg_id},${pos_id}}
	cat ${neg_fastq} |minimap2 -t ${threads}  -ax map-ont  --secondary=no ${ref} -  |samtools sort -@ {threads} -T bob -O bam -o  ${align_result}/$neg_id/$neg_id.sort.bam -  && samtools index -@ {threads} ${align_result}/$neg_id/$neg_id.sort.bam
	cat ${pos_fastq} |minimap2 -t ${threads}  -ax map-ont  --secondary=no ${ref} -  |samtools sort -@ {threads} -T bob -O bam -o  ${align_result}/$pos_id/$pos_id.sort.bam -  && samtools index -@ {threads} ${align_result}/$pos_id/$pos_id.sort.bam
}

#The read_ids of positive and negative samples are classified from candidate sites
label(){
	labeling_result="/home/weir/m6a_model/result/labeling_result/"
	mkdir -p ${labeling_result}
	CheckFolder $labeling_result
	py ./label_reads.py --bam $1 --label $2 --sites_file $3 --outFolder $labeling_result --plain_text_output
}


if (whiptail --title "Welcome using preprocessing workflow for train nano m6a model" --fb --backtitle "auhtor:liangou@ips.ac.cn" --yes-button "Continue" --no-button "Exit" --yesno \
	"This script will guide you generate feature table from raw_data." 10 120);then
	echo "begin to analyze..."
else
	exit 1
fi

OPTION=$(whiptail --title "The steps you want to complete" --fb --checklist "Select the options you want to do(Allow multiple)" 15 120 5 \
"unzip_data" "Unzip the zipped package containing the fast5 file" OFF \
"read_QC" "Implement quality control on input fastq files based on NanoPlot" OFF \
"mapping" "The sequence was aligned to the reference transcriptome" OFF \
"bam_QC" "evaluation of alignment result(bam files)"  OFF \
"label" "Label the reads that mapping the known m6A locus" OFF 3>&1 1>&2 2>&3 ) 

exitstatus=$?
if [ $exitstatus = 0 ]; then
    echo "Next we will complete the following steps:" $OPTION
   	for step in "${!OPTION[@]}";do
   		eval ${OPTION[$step]}
   	done

else
    echo "You chose Cancel."
fi
