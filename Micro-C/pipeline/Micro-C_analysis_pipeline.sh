#!/usr/bin/env bash
##-------------------------
# @Author: Siling Hu
# @Date: 2023/01/01
##-------------------------

function usage {
    echo -e "usage : Micro-C_analysis_pipeline.sh -i prefixes -p thread -w workdir [-h]"
	#echo -e "for example : Micro-C_analysis_pipeline.sh -i "K562_MicroC_DMSO_R1 K562_MicroC_DMSO_R2" -p 8 -w ~/data/" 
    echo -e "Use option -h|--help for more information"
}

function help {
    usage;
    echo 
    echo "Analyze Micro-C data. Firstly, you need to create a folder to use as your working path and create a \"1rawdata\" folder inside it to put in your fq.gz file."
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -i|--input prefixes : prefixes of Micro-C samples"
    echo "   -p|--thread : number of threads"
    echo "   -w|--workdir : workdir"
    echo "   [-h|--help] : help"
    exit;
}

if [ $# -lt 1 ]
then
    usage
    exit
fi

while getopts ":i:p:w:h" OPT
do
    case $OPT in
		
		i) prefixes=$OPTARG;;
        p) thread=$OPTARG;;
        w) workdir=$OPTARG;;
        h) help ;;		
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            exit 1
            ;;
    esac
done

# check parameter
if [ -z "${prefixes}" -o -z "${thread}" -o -z "${workdir}" ]; then
	printf "[ERROR] `date '+%F %T'` following parameters is empty:\n-i=${prefixes}\n-p=${thread}\n-w=${workdir}\n"
	exit 1
fi

echo -e "input prefixes: ${prefixes}\nthread: ${thread}\nworkdir: ${workdir}\n============================================"

echo -e "Micro-C data analysis, start!\n"

cd ${workdir}/


for prefix in ${prefixes}
do
	
	cd ${workdir}
	mkdir -p ${workdir}/2fastqc ${workdir}/3trim_adapter ${workdir}/4hicpro/fastqc 
    echo ${prefix}
    #qc before trim
    fastqc -o ${workdir}/2fastqc -q -f fastq ${workdir}/1rawdata/${prefix}*fq.gz -t ${thread}

    #trim adapter
    trim_galore -q 10 --phred33 --stringency 1 --length 10 -e 0.1 -j ${thread} --paired ${workdir}/1rawdata/${prefix}_1.fq.gz ${workdir}/1rawdata/${prefix}_2.fq.gz --gzip -o ${workdir}/3trim_adapter > ${workdir}/3trim_adapter/${prefix}_trimgalore.log 2>&1

    #trim linker
    cd ${workdir}/4hicpro
    mkdir -p ${prefix}/${prefix}/
    cutadapt -j ${thread} -q 15 -O 10 -m 10 -e 0.1 -a file:/mnt/disk3/husiling/private/pipeline/nasMicroC/MicroC_UMIlinker.fa -A file:/mnt/disk3/husiling/private/pipeline/nasMicroC/MicroC_UMIlinker.fa -o ${prefix}/${prefix}/tmp.${prefix}_1.fq.gz -p ${prefix}/${prefix}/tmp.${prefix}_2.fq.gz ${workdir}/3trim_adapter/${prefix}_1_val_1.fq.gz ${workdir}/3trim_adapter/${prefix}_2_val_2.fq.gz --info ${prefix}_trim.info --rename='{id};{r1.adapter_name}{r2.adapter_name} {comment}' --json=${prefix}_trim.json --times 10 > ${prefix}_trim.log

    cutadapt -j ${thread} -q 15 -O 10 -m 10 -e 0.1 -a file:/mnt/disk3/husiling/private/pipeline/nasMicroC/MicroC_UMIlinker_noAT.fa -A file:/mnt/disk3/husiling/private/pipeline/nasMicroC/MicroC_UMIlinker_noAT.fa -o ${prefix}/${prefix}/${prefix}_1.fq.gz -p ${prefix}/${prefix}/${prefix}_2.fq.gz ${prefix}/${prefix}/tmp.${prefix}_1.fq.gz ${prefix}/${prefix}/tmp.${prefix}_2.fq.gz --info ${prefix}_trim_noAT.info --rename='{id};{r1.adapter_name}{r2.adapter_name} {comment}' --json=${prefix}_trim_noAT.json --times 10 > ${prefix}_trim_noAT.log

    rm ${prefix}/${prefix}/tmp.*.fq.gz

    #qc after trim linker
    fastqc -o ${workdir}/4hicpro/fastqc -q -f fastq ${workdir}/4hicpro/*/*/${prefix}*fq.gz -t ${thread}

    #hicpro
    cd ${workdir}/4hicpro
    HiC-Pro -c /mnt/disk3/husiling/private/pipeline/nasMicroC/Micro-C_hicpro_hg38XX.config -i ${prefix} -o ${prefix}_hicpro

    echo ${prefix}  
	#Remove duplicates based on UMI and coordinates	
    LANG=en; awk -F "\t|;" '{$10=$2","$3;print $1";"$2","$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}' ${workdir}/4hicpro/${prefix}_hicpro/hic_results/data/${prefix}/${prefix}_hg38XX.bwt2pairs.validPairs | sort -T /mnt/disk3/husiling/private/data/tmp/ -S 50% -k2,2V -k3,3n -k5,5V -k6,6n -k8,8n -m | awk -F "\t|;" 'BEGIN{c1=0;c2=0;s1=0;s2=0;a1=0}(c1!=$3 || c2!=$6 || s1!=$4 || s2!=$7 || a1!= $2){print;c1=$3;c2=$6;s1=$4;s2=$7;a1=$2}' > ${workdir}/4hicpro/${prefix}_hicpro/hic_results/data/${prefix}/${prefix}.UMI.allValidPairs

	#statistics
    echo -e "${prefix}\t1rawdata:\t`zcat ${workdir}/1rawdata/${prefix}_1.fq.gz|awk 'END{print NR/4}'`" > ${workdir}/Statistics_${prefix}.txt
    echo -e "${prefix}\t2rm_adapter:\t`zcat ${workdir}/3trim_adapter/${prefix}_1_val_1.fq.gz|awk 'END{print NR/4}'`" >> ${workdir}/Statistics_${prefix}.txt
    echo -e "${prefix}\t3trim_linker:\t`zcat ${workdir}/4hicpro/${prefix}/${prefix}/${prefix}_1.fq.gz|awk 'END{print NR/4}'`" >> ${workdir}/Statistics_${prefix}.txt

    Unmapped_pairs=`grep Unmapped_pairs ${workdir}/4hicpro/${prefix}_hicpro/hic_results/stats/${prefix}/${prefix}.mpairstat|cut -f 2`
    Total_pairs_processed=`grep Total_pairs_processed ${workdir}/4hicpro/${prefix}_hicpro/hic_results/stats/${prefix}/${prefix}.mpairstat|cut -f 2`
    echo -e "${prefix}\t4Mapping:\t`echo -e "$Total_pairs_processed - $Unmapped_pairs"|bc`" >> ${workdir}/Statistics_${prefix}.txt

    Multiple_pairs_alignments=`grep Multiple_pairs_alignments ${workdir}/4hicpro/${prefix}_hicpro/hic_results/stats/${prefix}/${prefix}.mpairstat|cut -f 2`
    echo -e "${prefix}\t5rm_Multiple_pairs_alignments:\t`echo -e "$Total_pairs_processed - $Unmapped_pairs - $Multiple_pairs_alignments"|bc`" >> ${workdir}/Statistics_${prefix}.txt

    Pairs_with_singleton=`grep Pairs_with_singleton ${workdir}/4hicpro/${prefix}_hicpro/hic_results/stats/${prefix}/${prefix}.mpairstat|cut -f 2`
    echo -e "${prefix}\t6rm_Pairs_with_singleton:\t`echo -e "$Total_pairs_processed - $Unmapped_pairs - $Multiple_pairs_alignments -$Pairs_with_singleton"|bc`" >> ${workdir}/Statistics_${prefix}.txt

    Low_qual_pairs=`grep Low_qual_pairs ${workdir}/4hicpro/${prefix}_hicpro/hic_results/stats/${prefix}/${prefix}.mpairstat|cut -f 2`
    echo -e "${prefix}\t7rm_Low_qual_pairs:\t`echo -e "$Total_pairs_processed - $Unmapped_pairs - $Multiple_pairs_alignments -$Pairs_with_singleton -$Low_qual_pairs"|bc`" >> ${workdir}/Statistics_${prefix}.txt

    Unique_paired_alignments=`grep Unique_paired_alignments ${workdir}/4hicpro/${prefix}_hicpro/hic_results/stats/${prefix}/${prefix}.mpairstat|cut -f 2`
    echo -e "${prefix}\t7Unique_paired_alignments:\t$Unique_paired_alignments" >> ${workdir}/Statistics_${prefix}.txt

    Deduplication_by_coor=`wc -l ${workdir}/4hicpro/${prefix}_hicpro/hic_results/data/${prefix}/${prefix}.allValidPairs|cut -f 1 -d " "`
    echo -e "${prefix}\t8Deduplication_by_coor:\t$Deduplication_by_coor" >> ${workdir}/Statistics_${prefix}.txt

    Deduplication_by_UMI_coor=`wc -l ${workdir}/4hicpro/${prefix}_hicpro/hic_results/data/${prefix}/${prefix}.UMI.allValidPairs|cut -f 1 -d " "`
    echo -e "${prefix}\t9Deduplication_by_UMI_coor:\t$Deduplication_by_UMI_coor" >> ${workdir}/Statistics_${prefix}.txt

    awk -v prefix=$prefix '$4=="+" && $7=="+" {sum1+=1}$4=="+" && $7=="-" {sum2+=1}$4=="-" && $7=="+" {sum3+=1} $4=="-" && $7=="-" {sum4+=1}END{sum=sum1+sum2+sum3+sum4;print prefix"\t10valid_pairs:\t"sum"\t"sum/sum"\n"prefix"\t11FF:\t"sum1"\t"sum1/sum"\n"prefix"\t12FR:\t"sum2"\t"sum2/sum"\n"prefix"\t13RF:\t"sum3"\t"sum3/sum"\n"prefix"\t14RR:\t"sum4"\t"sum4/sum}' ${workdir}/4hicpro/${prefix}_hicpro/hic_results/data/${prefix}/${prefix}.UMI.allValidPairs >> ${workdir}/Statistics_${prefix}.txt

    awk -v prefix=$prefix '$2==$5 && $6-$3<=20000{cis_short+=1}$2==$5 && $6-$3>20000{cis_long+=1}$2!=$5{trans+=1}END{sum=cis_short+cis_long+trans;print prefix"\t15cis_short:\t"cis_short"\t"cis_short/sum"\n"prefix"\t16cis_long:\t"cis_long"\t"cis_long/sum"\n"prefix"\t17trans:\t"trans"\t"trans/sum}' ${workdir}/4hicpro/${prefix}_hicpro/hic_results/data/${prefix}/${prefix}.UMI.allValidPairs >> ${workdir}/Statistics_${prefix}.txt


    samtools sort -@ 80 ${workdir}/4hicpro/${prefix}_hicpro/bowtie_results/bwt2/${prefix}/${prefix}_hg38XX.bwt2pairs.bam -o ${workdir}/4hicpro/${prefix}_hicpro/bowtie_results/bwt2/${prefix}/tmp.${prefix}_hg38XX.bwt2pairs.bam

    /mnt/disk3/husiling/private/miniconda3/bin/picard MarkDuplicates -REMOVE_DUPLICATES true --VALIDATION_STRINGENCY SILENT --I ${workdir}/4hicpro/${prefix}_hicpro/bowtie_results/bwt2/${prefix}/tmp.${prefix}_hg38XX.bwt2pairs.bam --O ${workdir}/4hicpro/${prefix}_hicpro/bowtie_results/bwt2/${prefix}/${prefix}_hg38XX.bwt2pairs_rmdup.bam --M ${workdir}/4hicpro/${prefix}_hicpro/bowtie_results/bwt2/${prefix}/${prefix}_hg38XX.bwt2pairs_rmdup.metrics --OPTICAL_DUPLICATE_PIXEL_DISTANCE 10000 --MAX_OPTICAL_DUPLICATE_SET_SIZE 300000 --TAGGING_POLICY All

    nQ10=`samtools view -@ 80 -c ${workdir}/4hicpro/${prefix}_hicpro/bowtie_results/bwt2/${prefix}/${prefix}_hg38XX.bwt2pairs.bam | awk '{print $1/2}'`

    index=`less ${workdir}/1rawdata/${prefix}_1.fq.gz |head -n 1|cut -f 2 -d " "`
    READ_PAIR_DUPLICATES=`cat ${workdir}/4hicpro/${prefix}_hicpro/bowtie_results/bwt2/${prefix}/${prefix}_hg38XX.bwt2pairs_rmdup.metrics|grep -A 2 "LIBRARY"|cut -f 7-8|awk 'NR==2{print $1}'`

    READ_PAIR_OPTICAL_DUPLICATES=`cat ${workdir}/4hicpro/${prefix}_hicpro/bowtie_results/bwt2/${prefix}/${prefix}_hg38XX.bwt2pairs_rmdup.metrics|grep -A 2 "LIBRARY"|cut -f 7-8|awk 'NR==2{print $2}'`

    PCR_dup_rate=`echo "scale=2; 100*($READ_PAIR_DUPLICATES-$READ_PAIR_OPTICAL_DUPLICATES)/$nQ10" | bc`

    optical_dup_rate=`echo "scale=2; 100*$READ_PAIR_OPTICAL_DUPLICATES/$nQ10" | bc`

    echo -e "${prefix}\t18PCR_dup_rate:\t$PCR_dup_rate" >> ${workdir}/Statistics_${prefix}.txt

    rm ${workdir}/4hicpro/${prefix}_hicpro/bowtie_results/bwt2/${prefix}/tmp.${prefix}_hg38XX.bwt2pairs.bam

    rm ${workdir}/4hicpro/${prefix}_hicpro/bowtie_results/bwt2/${prefix}/${prefix}_hg38XX.bwt2pairs_rmdup.bam
done
