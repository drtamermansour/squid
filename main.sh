#!/bin/sh

## construction of the basic diretory structure
git clone https://github.com/drtamermansour/squid.git
cd squid
squid=$(pwd)

## create working directory and define paths for raw data and scripts
mkdir -p $squid/{data,QC_raw,adap_remove}
script_path=$squid/scripts
data_path=$squid/data
prepData=$squid/prepdata
QC_raw=$squid/QC_raw
trimmed_data=$squid/adap_remove

## get the data from SRA DB
## http://www.ncbi.nlm.nih.gov/sra/?term=pealei
SRA_URL=$"ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR"
cd $data_path
wget -r $SRA_URL/SRR172/SRR1725172/SRR1725172.sra
wget -r $SRA_URL/SRR172/SRR1725236/SRR1725236.sra
wget -r $SRA_URL/SRR172/SRR1725235/SRR1725235.sra
wget -r $SRA_URL/SRR172/SRR1725213/SRR1725213.sra
wget -r $SRA_URL/SRR172/SRR1725171/SRR1725171.sra
wget -r $SRA_URL/SRR172/SRR1725169/SRR1725169.sra
wget -r $SRA_URL/SRR172/SRR1725167/SRR1725167.sra
wget -r $SRA_URL/SRR172/SRR1725164/SRR1725164.sra
wget -r $SRA_URL/SRR172/SRR1725163/SRR1725163.sra
wget -r $SRA_URL/SRR152/SRR1522988/SRR1522988.sra
wget -r $SRA_URL/SRR152/SRR1522987/SRR1522987.sra
#wget -r $SRA_URL/SRR152/SRR1522984/SRR1522984.sra ## genomic DNA

##Data Striping on Scratch
#lfs setstripe --count -1 $prepData
cd $prepData
for dir in $data_path/$SRA_URL/*/*; do if [ -d $dir ]; then
  echo $dir
  qsub -v inputdir=${dir} ${script_path}/fastq-dump.sh
fi done;

## change the name of the data files to match the standard format
for R1 in *_1.fastq.gz; do
  echo $R1; R2=$(echo $R1 | sed s/_1.fastq.gz/_2.fastq.gz/);
  newR1=$(echo $R1 | sed s/_1.fastq.gz/_R1_001.fastq.gz/);
  newR2=$(echo $R1 | sed s/_1.fastq.gz/_R2_001.fastq.gz/);
  mv $R1 $newR1; mv $R2 $newR2; done;
chmod u-w $prepData/*.fastq.gz

## QC check for the original data files 
cd $QC_raw
for f in ${prepData}/*.gz; do qsub -v INPUT_FILE="$f" ${script_path}/fastqc.sh; done
mv ${prepData}/{*.zip,*.html} $QC_raw/.
#============================================================================================
## adaptor removal only 
cd ${trimmed_data}
bash ${script_path}/run_adapter_trimmer.sh ${prepData} ${trimmed_data} ${script_path}/adapter_removal.sh

## combine single reads 
for f in ${trimmed_data}/*; do if [ -d "${f}" ]; then echo "${f}"; cat "${f}"/s1_se.fq "${f}"/s2_se.fq > "${f}"/$(basename "$f").s_se.fq; fi; done;

## QC check
#for dir in ${trimmed_data}/*; do cd $dir; for f in *_pe.fq *.s_se.fq; do qsub -v INPUT_FILE="$f" ${script_path}/fastqc.sh; done; done;
#============================================================================================
## Abundance estimation
extTrans=$squid/extTrans
mkdir -p $extTrans
cd $extTrans ## create the file "targetTrans.fa"
qsub -v index="salmon_index",transcriptome="$extTrans/targetTrans.fa" ${script_path}/salmonIndex2.sh
for dir in ${trimmed_data}/*; do if [ -d "$dir" ]; then
  cd $dir
  identifier=$(basename $dir).PE
  qsub -v index="$extTrans/salmon_index",identifier=$identifier ${script_path}/salmonQuant_PE.sh  
  identifier=$(basename $dir).SE
  qsub -v index="$extTrans/salmon_index",identifier=$identifier ${script_path}/salmonQuant_SE.sh
fi;done
cd ${trimmed_data}
find */*.quant -name *.sf -exec grep -H "mapping rate" {} \; | sort > salmonQuant_summary.txt
python $script_path/gather-counts.py -i "$(pwd)"
echo "transcript"$'\t'"length" > transcripts.lengthes
sf=$(find */*.quant -name \*.sf | head -n1)
cat $sf | grep -v "^#" | awk -F "\t" -v OFS='\t' '{print $1,$2}' >> transcripts.lengthes
grep "^>" $extTrans/targetTrans.fa | sed 's/>//g' > gene_transcript.map
for dir in $trimmed_data/*;do if [ -d $dir ];then cd $dir; wl=$(cat s1_pe.fq *.s_se.fq | wc -l); c=$(($wl/4)); echo $(basename $dir) $c;fi; done > readCounts
module load R/3.0.1
for target in *;do if [ -d $target ];then
  echo $target
  inputreads=$(grep $target readCounts | awk '{print $2}')
  Rscript ${script_path}/calcTPM_tis_noIso.R "$(pwd)" "$target" "$inputreads" >> targets_list # > /dev/null
fi; done 
for f in *.dataSummary_comp;do echo $f; cat $f;done > dataSummary_comp.Summary
#===========================================================================================
## prepare the files for subsequent analysis
for dir in ${trimmed_data}/*; do if [ -d "${dir}" ]; then echo $dir; cd $dir; qsub -v output=$(basename "$dir").s_pe.fq ${script_path}/interleave.sh; fi; done;

## filter abundance
mkdir $squid/abundFilter
cd $squid/abundFilter
input_files=()
for f in ${trimmed_data}/*/*.s_[ps]e.fq; do input_files+=($f); done
qsub -v graph=$"countgraph_k20.kt",input_files="${input_files[*]}" ${script_path}/load-into-counting.sh
qsub -v input=$"countgraph_k20.kt",files="${input_files[*]}" ${script_path}/filter_abund.sh
for f in filter_abund.e*;do grep -B2 '^output in' "$f";done > filter_abund.summary
## break out the orphaned and still-paired reads and rename files (this step ends with .s_pe.fq and .s_se.fq for each sample)
#for i in *.s_pe.*.abundfilt; do extract-paired-reads.py $i; done
for i in *.s_pe.*.abundfilt; do qsub -v input=$i $script_path/extract-paired-reads.sh; done
##  combine the orphaned reads into a single file & rename pe files
for i in *.s_se.fq.abundfilt; do
   pe_orphans=$(basename $i .s_se.fq.abundfilt).s_pe.fq.abundfilt.se
   cat $i $pe_orphans > $(basename $i .abundfilt)
done
#rm *.abundfilt.se
for i in *.abundfilt.pe; do mv $i $(basename $i .abundfilt.pe); done
## split interleaved files & ## merge the single reads to the end of left reads (this step reform the data into .s_pe1_se.fq & .s_pe2.fq for each sample )
for i in *.s_pe.fq; do qsub -v input=$i $script_path/split-paired-reads.sh; done
for f in *.s_pe.fq; do
  base=$(basename $f);
  file_end=$(tail -n 8 $f | head -n 4);
  file_endA=$(tail -n 4 $base.1);
  if [ "$file_end" != "$file_endA" ];then echo $f.1;fi
  file_end=$(tail -n 4 $f);
  file_endB=$(tail -n 4 $base.2);
  if [ "$file_end" != "$file_endB" ];then echo $f.2;fi
done

for f in *.s_pe.fq.1; do echo $f; base=$(basename $f .s_pe.fq.1); cat $f $base.s_se.fq > $base.s_pe1_se.fq; done
#rm *.s_pe.fq.1
for f in *.s_pe.fq.2; do mv $f $(basename $f .s_pe.fq.2).s_pe2.fq; done

## calculation of filtered read counts
for f in *.s_pe1_se.fq ; do wc -l $f;done > filtered_readCounts
cat filtered_readCounts | awk '{ sum+=$1} END {print sum,sum/4}' ## 454,506,826

## Trinity
lf_files=()
for f in *.s_pe1_se.fq; do lf_files+=($f); done;
rt_files=()
for f in *.s_pe2.fq; do rt_files+=($f); done;
qsub -v lf="${lf_files[*]}",rt="${rt_files[*]}" ${script_path}/run_Trinity.sh
echo $(grep "^>" Trinity.fasta | wc -l) ## 881402

