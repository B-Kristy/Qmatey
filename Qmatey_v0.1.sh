#!/bin/bash
[ $# -ge 1 -a -f "$1" ] && input="$1" || { echo -e "$1 \e[31mnot found, exiting"; exit 1; }
. $input
YELLOW='\033[1;33m'
WHITE='\033[1;37m'
proj_dir=${input%/*}

if java -version; then
	:
else
	echo "java not detected, try 'sudo apt install default-jre'"
	exit 1
fi

if datamash --version; then
	:
else
	echo "datamash not detected, try 'sudo apt install datamash"
	exit 1
fi 


##################################################################################################################
#Check for existence of directories specified in config file
[ -d "$tool_dir" ] || { echo -e "$tool_dir \e[31mnot found, exiting"; exit 1; }
[ -d "$input_dir" ] || { echo -e "$input_dir \e[31mnot found, exiting"; exit 1; }
[ -d "$proj_dir" ] || { echo -e "$proj_dir \e[31mnot found, exiting"; exit 1; }
[ -d "$ref_dir" ] || { echo -e "$ref_dir \e[31mnot found, exiting"; exit 1; }
##################################################################################################################
#Create all necessary subdirectories for downstream processing
cd $proj_dir
mkdir metagenome || { echo -e "\e[31mproject directory can't contain folder named metagenome, exiting"; exit 1; }
cd $proj_dir/metagenome
mkdir haplotig
mkdir alignment
mkdir sighits
mkdir results

#################################################################################################################
#Make dictionary of host reference genome
echo -e "${YELLOW}------------------------------------------------------------------------------ \n \n Qmatey is compiling the host reference genome \n \n------------------------------------------------------------------------------${WHITE}"
cd $ref_dir
cat *.fa* > master_ref.fasta
$tool_dir/bwa-0.7.17/bwa index -a bwtsw master_ref.fasta &>/dev/null
java -jar $tool_dir/picard-2.18.27.jar CreateSequenceDictionary REFERENCE=master_ref.fasta    OUTPUT=master_ref.dict &>/dev/null 
$tool_dir/samtools-1.9/samtools faidx master_ref.fasta 
#################################################################################################################
#Perform alignment of sample reads to host reference genome
echo -e "${YELLOW}------------------------------------------------------------------------------ \n \n Qmatey is aligning raw reads to the reference genome \n \n------------------------------------------------------------------------------${WHITE}"
cd $input_dir
for i in $(ls *.fa*); do
	$tool_dir/bwa-0.7.17/bwa mem -t $threads $ref_dir/master_ref.fasta $i > $proj_dir/metagenome/${i%.fa*}.sam && \
	java -XX:ParallelGCThreads=$threads -jar $tool_dir/picard-2.18.27.jar SortSam I= $proj_dir/metagenome/${i%.fa*}.sam O= $proj_dir/metagenome/${i%.fa*}.bam SORT_ORDER=coordinate && \
	rm $proj_dir/metagenome/${i%.fa*}.sam 
done &>/dev/null
#################################################################################################################
#Extract unmatched (probably mostly microbiome) reads
echo -e "${YELLOW}------------------------------------------------------------------------------ \n \n Qmatey is identifying metagenomic reads \n \n------------------------------------------------------------------------------${WHITE}"
cd $proj_dir/metagenome
for i in $(ls *.bam); do
	$tool_dir/samtools-1.9/samtools view -b -f4 $i > ${i%.*}_metagenome.bam && \
	$tool_dir/samtools-1.9/samtools bam2fq ${i%.*}_metagenome.bam | gzip > ${i%.*}_metagenome.fastq.gz && \
	rm ${i%.*}_metagenome.bam
done &>/dev/null
##################################################################################################################
#Uses word count to determine total amount of starting sequences. This will be used to calculate coverage.
echo -e "${YELLOW}------------------------------------------------------------------------------ \n \n Qmatey is calculating metagenomic coverage \n \n------------------------------------------------------------------------------"

cd $proj_dir/metagenome
for i in $(ls *.fastq.gz); do
	gunzip $i
done

cd $proj_dir/metagenome
for i in $(ls *.fastq); do
	awk '{if(NR%4==2) print $0}' $i > ./results/${i%.*}_seq.txt
done

cd $proj_dir/metagenome
gzip *.fastq

cd $proj_dir/metagenome/results
for i in $(ls *_seq.txt); do
	wc -l $i > ${i%meta*}seqcov.txt
	find . -type f -name '*_seqcov.txt' -exec cat {} + > tempcov.txt
	awk '{sum+=$1}END{print sum}' tempcov.txt > seqcov.txt
done


rm  *_seq.txt && rm *tempcov.txt && rm *_seqcov.txt
##################################################################################################################
#Compares unmatched microbiome (metagenome) reads to the master reference (host genomes) to calculate coverage normalization factor
echo -e "${YELLOW}------------------------------------------------------------------------------ \n \n Qmatey is calculating a host-coverage normalization factor \n \n------------------------------------------------------------------------------"
cd $proj_dir/metagenome
touch host_coverage.txt
touch empty_host_coverage.txt
touch microbiome_coverage.txt
touch empty_microbiome_coverage.txt
for i in $(ls *.bam); do
	$tool_dir/samtools-1.9/samtools view -c -F 4 ${i%.*}.bam > ${i%.*}_host_coverage.txt && \
	grep "" *_host_coverage.txt > appendhost_coverage.txt && \
	rm ${i%.*}_host_coverage.txt && cat host_coverage.txt appendhost_coverage.txt > holdhost_coverage.txt  && \
	cat holdhost_coverage.txt > host_coverage.txt
	$tool_dir/samtools-1.9/samtools view -c -f 4 ${i%.*}.bam > ${i%.*}_microbiome_coverage.txt && \
	grep "" *_microbiome_coverage.txt > appendmicrobiome_coverage.txt && \
	rm ${i%.*}_microbiome_coverage.txt && cat microbiome_coverage.txt appendmicrobiome_coverage.txt > holdmicrobiome_coverage.txt  && \
	cat holdmicrobiome_coverage.txt > microbiome_coverage.txt
done
awk  'gsub("_host_coverage.txt:", "\t", $0)'  holdhost_coverage.txt > host_coverage.txt
awk  'gsub("_microbiome_coverage.txt:", "\t", $0)'  holdmicrobiome_coverage.txt > microbiome_coverage.txt
join host_coverage.txt microbiome_coverage.txt > coverage_normalize.txt
average=$(awk '{ total += $2 } END { print total/NR }' coverage_normalize.txt)
awk -v average=$average '{print $1,$2/average}' coverage_normalize.txt > coverage_normalization_factor.txt
rm empty* && rm hold* && rm append* && rm host* && rm microbiome*
#################################################################################################################
#Processes fastq input (metagenomic reads) files into fasta format for further processing
cd $proj_dir/metagenome
for i in $(ls *_metagenome.fastq.gz); do
	gunzip -k $i && mv ${i%_meta*}_metagenome.fastq ./haplotig/${i%_meta*}_metagenome.fastq && cd ./haplotig
	awk 'NR%2==0' ${i%_meta*}_metagenome.fastq | awk 'NR%2==1' | sort | uniq -c | tr -s ' ' '\t' | awk '$0=">hseqid_"NR$0' | awk '$2>"0"' > ${i%_meta*}_cov_haplotig.txt
	awk '{print $1"\t"$2}' ${i%_meta*}_cov_haplotig.txt | awk  'gsub(">", "", $0)' > ${i%_meta*}_haplocov.txt
	awk '{print $1"\n"$3}' ${i%_meta*}_cov_haplotig.txt > ${i%_meta*}_haplotig.fasta
	rm ${i%_meta*}_metagenome.fastq && rm ${i%_meta*}_cov_haplotig.txt && cd $proj_dir/metagenome
done
##################################################################################################################
#Use standard/host-normalization-factor (coverage normalization factor) and DNA-normalization/dilution-factor (in house spike) to normalize count data
#If user has normalization factor they must provide it, if they don't a normalization factor of 1 will be used
cd $proj_dir/metagenome/haplotig
for i in $(ls *_haplocov.txt);do
	id=${i%_haplocov*}
	if [ -f "$normfactor" ]; then 
		libnormfactor=$(awk -v id=$id '$1 == id {print $2}' $normfactor) 
	else libnormfactor=1; 
	fi
	awk -v libnormfactor=$libnormfactor '{print $1,$2 * libnormfactor}' $i | awk 'gsub(" ", "\t", $0)' > ${i%_haplocov*}_libnormalized.txt
	internalnormfactor=$(awk -v id=$id '$1 == id {print $2}' ../coverage_normalization_factor.txt)
	awk -v internalnormfactor=$internalnormfactor '{print $1,$2 * internalnormfactor}' ${i%_haplocov*}_libnormalized.txt | awk 'gsub(" ", "\t", $0)' > ${i%_haplocov*}_normalized.txt
done
##################################################################################################################
#install blast ### https://www.exoscale.com/syslog/blast/
#permanently set environment’s PATH variable: "export PATH=$PATH:/home/bode/ncbi-blast-2.8.1+/bin/"
#download ncbi database (16SMicrobial, nt, taxdump, taxdb) >>> https://github.com/yeban/ncbi-blast-dbs
#updatde database ### email="my email address here" ncbi-blast-dbs 16SMicrobial nt taxdump
#specify blast output format ### http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
#download new_taxdump.tar.gz https://ncbiinsights.ncbi.nlm.nih.gov/2018/02/22/new-taxonomy-files-available-with-lineage-type-and-host-information/
echo -e "${YELLOW}------------------------------------------------------------------------------ \n \n Qmatey is performing BLAST \n \n------------------------------------------------------------------------------"
cd $proj_dir/metagenome/haplotig
for i in $(ls *_haplotig.fasta);do
	$tool_dir/ncbi-blast-2.8.1+/bin/blastn -task megablast -query $i -db $db_dir -num_threads $threads -evalue 1e-10 -max_target_seqs 5 -outfmt \
	"6 qseqid sseqid length mismatch evalue pident qcovs qseq sseq staxids stitle" \
	-out ../alignment/${i%_haplotig*}_haplotig.megablast
done
##################################################################################################################
#Removes duplicate rows and vector contamination from *_haplotig.megablast
cd $proj_dir/metagenome/alignment
for i in $(ls *_haplotig.megablast);do
	sort -u $i > ${i%_haplotig*}_haplotig_nd.megablast
done
rm *_haplotig.megablast

##################################################################################################################
#Strain-level sighit identification
if [ "$strain_level" == "TRUE" ]; then
cd $proj_dir/metagenome/sighits
mkdir sighits_strain
cd $proj_dir/metagenome/results
mkdir strain_level
cp seqcov.txt ./strain_level

	echo -e "${YELLOW}------------------------------------------------------------------------------ \n \n Qmatey is performing strain-level analysis \n \n------------------------------------------------------------------------------"
	cd $proj_dir/metagenome/alignment
	for i in $(ls *_haplotig_nd.megablast);do
		awk '$6>=100' $i | awk '$7>=100' | awk 'gsub(" ","_",$0)' > ../sighits/sighits_strain/${i%_haplotig*}_filter.txt
		awk 'NR==FNR {h[$1] = $2; next} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,h[$1]}' ../haplotig/${i%_haplotig*}_normalized.txt ../sighits/sighits_strain/${i%_haplotig*}_filter.txt | 
		awk '{print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' | awk 'gsub(" ","\t",$0)' > ../sighits/sighits_strain/${i%_haplotig*}_sighits.txt
	done
	rm ../sighits/sighits_strain/*_filter.txt 
##################################################################################################################
#Removes unique reads that match multiple subjects and adds headers to each column
	cd $proj_dir/metagenome/sighits/sighits_strain
	for i in $(ls *_sighits.txt);do
		awk '{a[$2]++;b[$2]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i > ${i%.txt}_nr_temp.txt
		echo $'abundance\tqseqid\tsseqid\tlength\tmismatch\tevalue\tpident\tqcovs\tqseq\tsseq\tstaxids\tstitle' | \
		cat - ${i%.txt}_nr_temp.txt > ${i%.txt}_nr.txt
	done
	rm *_sighits.txt *_nr_temp.txt
#################################################################################################################
#Combine all taxids for all files/individuals and perform single search against new_taxdump.
	awk '{gsub(/\t\t/,"\tNA\t"); print}' $tool_dir/rankedlineage.dmp | awk '{gsub(/[|]/,""); print}' | awk '{gsub(/\t\t/,"\t"); print}' > $tool_dir/rankedlineage_tabdelimited.dmp
	echo $'tax_id\ttaxname\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tsuperkingdom\t' | \
	cat - $tool_dir/rankedlineage_tabdelimited.dmp > $tool_dir/rankedlineage_edited.dmp
	rm $tool_dir/rankedlineage_tabdelimited.dmp

	cd $proj_dir/metagenome/sighits/sighits_strain
	find . -type f -name '*_sighits_nr.txt' -exec cat {} + > sighits.txt
	awk '{print $11}' sighits.txt | awk '{gsub(";","\n"); print}' | sort -u -n | sed -e '1s/staxids/tax_id/' > taxids_sighits.txt && rm sighits.txt
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}'  $tool_dir/rankedlineage_edited.dmp taxids_sighits.txt | \
	awk '{gsub(/ /,"_"); print }' > rankedlineage_subhits.txt 

	rm taxids_sighits.txt
#################################################################################################################
#Install datamash
#Generate file with mean, number of unique reads per taxID, and standard error
#Now, perform merge subsetted new_taxdump with each file while retaining only taxids in each file
	cd $proj_dir/metagenome/sighits/sighits_strain/
	awk '{print $1}' rankedlineage_subhits.txt > ../../results/strain_level/strain_taxa_mean.txt
	awk '{print $1}' rankedlineage_subhits.txt > ../../results/strain_level/strain_taxa_unique_sequences.txt
	awk '{print $1}' rankedlineage_subhits.txt > ../../results/strain_level/strain_taxa_quantification_accuracy.txt 
	for i in $(ls *_sighits_nr.txt); do
		cut -f 1,11 $i | awk '{print $2,"\t",$1}' | datamash --header-in --sort --group 1 mean 2 sstdev 2 count 2 | \
		awk '{ print $1,"\t",$2,"\t",$4,"\t",((($3/sqrt($4))/$2)*100) }' > stats1.txt
		echo $'tax_id\tmean\tuniq_reads\tstderr' | cat - stats1.txt > stats2.txt
		id=${i%_sighits*}_mean && awk -v id=$id '{gsub(/mean/,id); print }' stats2.txt | awk '{print $1,"\t",$2}' > holdmean.txt
		awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdmean.txt ../../results/strain_level/strain_taxa_mean.txt > holdmean2.txt && cat holdmean2.txt > ../../results/strain_level/strain_taxa_mean.txt
		id=${i%_sighits*}_uniq_reads && awk -v id=$id '{gsub(/uniq_reads/,id); print }' stats2.txt | awk '{print $1,"\t",$3}' > holduniq_reads.txt
		awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holduniq_reads.txt ../../results/strain_level/strain_taxa_unique_sequences.txt > holduniq_reads2.txt && cat holduniq_reads2.txt > ../../results/strain_level/strain_taxa_unique_sequences.txt
		id=${i%_sighits*}_stderr && awk -v id=$id '{gsub(/stderr/,id); print }' stats2.txt | awk '{print $1,"\t",$4}' > holdstderr.txt
		awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdstderr.txt ../../results/strain_level/strain_taxa_quantification_accuracy.txt > holdstderr2.txt && cat holdstderr2.txt > ../../results/strain_level/strain_taxa_quantification_accuracy.txt
		awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $1,"\t",$2,"\t",$3,"\t",$4,"\t", a[$1]}'  rankedlineage_subhits.txt stats2.txt > stats3.txt
		awk '{print $1,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14 }' stats3.txt | awk '{gsub(/ /,"\t"); print }' > ${i%_sighits*}_taxastats.txt
		rm *stats1* *stats2* *stats3* *hold*
	done
#################################################################################################################
#Calculates percent coverage of sighits from unmatched microbiome input files
	cd $proj_dir/metagenome/results/strain_level
	i="_mean$"
	awk -vp="$i" 'NR==1{for(i=1; i<=NF; i++) if ($i~p) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' strain_taxa_mean.txt > temp_mean.txt
	i="uniq_reads$"
	awk -vp="$i" 'NR==1{for(i=1; i<=NF; i++) if ($i~p) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' strain_taxa_unique_sequences.txt > temp_uniq.txt
	paste temp_mean.txt temp_uniq.txt | awk '/^[0-9]/ {for(i=1; i<=NF/2; i++) {s=s OFS $i*$(NF/2+i); }sub(/^ /,x,s);$0=s; s=""} !/[0-9]/{$0=$1;}1' > temp_uniq_mean.txt
	tail -n +2 temp_uniq_mean.txt > temp_uniq_mean_2.txt
	awk '{for (i=1; i<=NF; i++) sum[i]+=$i;}; END{for (i in sum) print sum[i];}' temp_uniq_mean_2.txt > temp_uniq_mean_3.txt
	awk '{sum+=$1}END{print sum}' temp_uniq_mean_3.txt > mean_uniq.txt
	paste mean_uniq.txt seqcov.txt | awk '{print(($1/$2)* 100)}' > strain_percent_coverage.txt
	rm *temp* *mean_uniq*
################################################################################################################
	cd $proj_dir/metagenome/results/strain_level
	for i in {mean,unique_sequences,quantification_accuracy}; do
		awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_strain/rankedlineage_subhits.txt strain_taxa_${i}.txt | \
		awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' > strain_taxainfo_${i}.txt
	done
	rm *_taxa_*
################################################################################################################
#Strain-level visualizations
	cd $proj_dir/metagenome/results/strain_level
	strain_level_mean=strain_taxainfo_mean.txt
	strain_level_uniq=strain_taxainfo_unique_sequences.txt
	strain_level_stderr=strain_taxainfo_quantification_accuracy.txt
	percent_thresh=5
	Rscript $tool_dir/Rscripts/strain_level_corr.R $strain_level_mean $percent_thresh &>/dev/null
	Rscript $tool_dir/Rscripts/strain_level_boxplots.R $strain_level_mean $strain_level_uniq $strain_level_stderr $percent_thresh &>/dev/null
fi
################################################################################################################
#Species-level sighit identification
if [ "$species_level" == "TRUE" ];then
cd $proj_dir/metagenome/sighits
mkdir sighits_species
cd $proj_dir/metagenome/results
mkdir species_level
cp seqcov.txt ./species_level

	echo -e "${YELLOW}------------------------------------------------------------------------------ \n \n Qmatey is performing species-level analysis \n \n------------------------------------------------------------------------------"
	cd $proj_dir/metagenome/alignment
	for i in $(ls *_haplotig_nd.megablast);do
		awk '$6>=95' $i | awk '$7>=95' | awk 'gsub(" ","_",$0)' > ../sighits/sighits_species/${i%_haplotig*}_filter.txt
		awk 'NR==FNR {h[$1] = $2; next} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,h[$1]}' ../haplotig/${i%_haplotig*}_normalized.txt ../sighits/sighits_species/${i%_haplotig*}_filter.txt | 
		awk '{print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' | awk 'gsub(" ","\t",$0)' > ../sighits/sighits_species/${i%_haplotig*}_sighits.txt
	done
	rm ../sighits/sighits_species/*_filter.txt 
################################################################################################################
#Extracts non-redundant hits for unique alignment
	cd $proj_dir/metagenome/sighits/sighits_species
	for i in $(ls *_sighits.txt);do
		awk -F '\t' '{a[$2]++;b[$2]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i > ${i%_sighits*}_species_unique_reads.txt
	done
################################################################################################################
#Extracts duplicate hits for multi-alignment 
	cd $proj_dir/metagenome/sighits/sighits_species
	for i in $(ls *_sighits.txt);do
		awk -F '\t' 'FNR==NR{a[$1,$2]=1; next}  !a[$1,$2]' ${i%_sighits*}_species_unique_reads.txt $i OFS='\t' > ${i%_sighits*}_dup.txt
	done
################################################################################################################
#Appends taxa informations to the duplicate sighits file
	cd $proj_dir/metagenome/sighits/sighits_species
	for i in $(ls *_dup.txt);do
		awk -F '\t' '{print $11}' OFS=';' $i > ${i%_dup*}_taxids_dup_inter.txt
	done

	for i in $(ls *_dup_inter.txt);do
	awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_dup_inter*}_taxids_dup.txt
	done

	rm *_taxids_dup_inter.txt
################################################################################################################
#Reformats taxonomic information to provide species-level taxa in the species column for duplicate reads
	cd $proj_dir/metagenome/sighits/sighits_species
	for i in $(ls *_taxids_dup.txt);do
		awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' $tool_dir/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_dup*}_dup_inter.txt 
	done

	for i in $(ls *_dup_inter.txt);do
		awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_dup_inter*}_species_taxid.txt
	done

	rm *_taxids_dup.txt

	for i in $(ls *_species_taxid.txt);do
		awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_species_taxid*}_species_column.txt
	done

	for i in $(ls *_species_column.txt);do
		paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column*}_species_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_species_taxa.txt
	done

	for i in $(ls *_dup.txt);do
		paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' OFS='\t' ${i%*_dup.txt}_species_taxa.txt) > ${i%_dup*}_species_duplicates_virome.txt
	done

	for i in $(ls *_species_duplicates_virome.txt);do
		awk -F '\t' '{ if ($21!="Viruses") print $0}' $i > ${i%*_species_duplicates_virome*}_species_duplicates.txt
	done

	rm *_species_taxid.txt && rm *_dup_inter.txt && rm *_dup.txt && rm *_species_column.txt && rm *_species_taxa.txt && rm *_species_duplicates_virome.txt
#################################################################################################################
#species-level clustering
	for i in $(ls *_species_duplicates.txt);do
		awk -F '\t' '{print $1, $2"~"$14, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_species_duplicates*}_species_inter.txt
	done

	for i in $(ls *_species_inter.txt);do
		awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_species_inter*}_species_inter2.txt
	done

	for i in $(ls *_species_inter2.txt);do
		awk -F '\t' '{dups[$1]++} END {for (num in dups) {print num}}' $i | sort -k1,1  > ${i%_species_inter2*}_duplicate_count.txt 
	done

	for i in $(ls *_duplicate_count.txt);do
		awk -F '~' '{a[$1]++;b[$1]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i | sort -k1,1 > ${i%_duplicate_count*}_multialign_species_reads.txt
	done

	for i in $(ls *_multialign_species_reads.txt);do
		awk -F '\t'  'FNR==NR {a[$1]; next}; $1 in a' $i ${i%_multialign_species_reads*}_species_inter2.txt | sort -u -k1,1 | awk 'gsub("~","\t",$0)'| awk -F '\t' '{print $3, $1, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' OFS='\t' > ${i%_multialign_species_reads*}_species_OTU.txt
	done

	rm *_species_inter.txt && rm *_species_inter2.txt && rm *_duplicate_count.txt && rm *_multialign_species_reads.txt && rm *_species_duplicates.txt

	for i in $(ls *_species_OTU.txt);do
		cat $i ${i%_species_OTU*}_species_unique_reads.txt > ${i%_species_OTU*}_complete_species_reads.txt
	done 

	for i in $(ls *_complete_species_reads.txt);do
		awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' $i > ${i%_complete_species_reads*}_sighits_temp.txt
		echo $'abundance\tqseqid\tsseqid\tlength\tmismatch\tevalue\tpident\tqcovs\tqseq\tsseq\tstaxids\tstitle' | \
		cat - ${i%_complete_species_reads*}_sighits_temp.txt > ${i%_complete_species_reads*}_sighits.txt
	done

	rm *_complete_species_reads.txt && rm *_sighits_temp.txt && rm *_species_OTU.txt && rm *_unique_reads.txt
#################################################################################################################
#Combine all taxids for all files/individuals and perform single search against new_taxdump.
	cd $proj_dir/metagenome/sighits/sighits_species
	find . -type f -name '*_sighits.txt' -exec cat {} + > sighits.txt
	awk '{print $11}' sighits.txt | awk '{gsub(";","\n"); print}' | sort -u -n | sed -e '1s/staxids/tax_id/' > taxids_sighits.txt && rm sighits.txt
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}'  $tool_dir/rankedlineage_edited.dmp taxids_sighits.txt | \
	awk '{gsub(/ /,"_"); print }' > rankedlineage_subhits.txt 

	rm taxids_sighits.txt
#################################################################################################################
#Install datamash
#Generate file with mean, number of unique reads per taxID, and standard error
#Now, perform merge subsetted new_taxdump with each file while retaining only taxids in each file
	cd $proj_dir/metagenome/sighits/sighits_species/
	awk '{print $1}' rankedlineage_subhits.txt > ../../results/species_level/species_taxa_mean.txt
	awk '{print $1}' rankedlineage_subhits.txt > ../../results/species_level/species_taxa_unique_sequences.txt
	awk '{print $1}' rankedlineage_subhits.txt > ../../results/species_level/species_taxa_quantification_accuracy.txt 
	for i in $(ls *_sighits.txt); do
		cut -f 1,11 $i | awk '{print $2,"\t",$1}' | datamash --header-in --sort --group 1 mean 2 sstdev 2 count 2 | \
		awk '{ print $1,"\t",$2,"\t",$4,"\t",((($3/sqrt($4))/$2)*100) }' > stats1.txt
		echo $'tax_id\tmean\tuniq_reads\tstderr' | cat - stats1.txt > stats2.txt
		id=${i%_sighits*}_mean && awk -v id=$id '{gsub(/mean/,id); print }' stats2.txt | awk '{print $1,"\t",$2}' > holdmean.txt
		awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdmean.txt ../../results/species_level/species_taxa_mean.txt > holdmean2.txt && cat holdmean2.txt > ../../results/species_level/species_taxa_mean.txt
		id=${i%_sighits*}_uniq_reads && awk -v id=$id '{gsub(/uniq_reads/,id); print }' stats2.txt | awk '{print $1,"\t",$3}' > holduniq_reads.txt
		awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holduniq_reads.txt ../../results/species_level/species_taxa_unique_sequences.txt > holduniq_reads2.txt && cat holduniq_reads2.txt > ../../results/species_level/species_taxa_unique_sequences.txt
		id=${i%_sighits*}_stderr && awk -v id=$id '{gsub(/stderr/,id); print }' stats2.txt | awk '{print $1,"\t",$4}' > holdstderr.txt
		awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdstderr.txt ../../results/species_level/species_taxa_quantification_accuracy.txt > holdstderr2.txt && cat holdstderr2.txt > ../../results/species_level/species_taxa_quantification_accuracy.txt
		awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $1,"\t",$2,"\t",$3,"\t",$4,"\t", a[$1]}'  rankedlineage_subhits.txt stats2.txt > stats3.txt
		awk '{print $1,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14 }' stats3.txt | awk '{gsub(/ /,"\t"); print }' > ${i%_sighits*}_taxastats.txt
		rm *stats1* *stats2* *stats3* *hold*
	done
#################################################################################################################
#Calculates percent coverage of sighits from unmatched microbiome input files
	cd $proj_dir/metagenome/results/species_level/
	i="_mean$"
	awk -vp="$i" 'NR==1{for(i=1; i<=NF; i++) if ($i~p) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' species_taxa_mean.txt > temp_mean.txt
	i="uniq_reads$"
	awk -vp="$i" 'NR==1{for(i=1; i<=NF; i++) if ($i~p) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' species_taxa_unique_sequences.txt > temp_uniq.txt
	paste temp_mean.txt temp_uniq.txt | awk '/^[0-9]/ {for(i=1; i<=NF/2; i++) {s=s OFS $i*$(NF/2+i); }sub(/^ /,x,s);$0=s; s=""} !/[0-9]/{$0=$1;}1' > temp_uniq_mean.txt
	tail -n +2 temp_uniq_mean.txt > temp_uniq_mean_2.txt
	awk '{for (i=1; i<=NF; i++) sum[i]+=$i;}; END{for (i in sum) print sum[i];}' temp_uniq_mean_2.txt > temp_uniq_mean_3.txt
	awk '{sum+=$1}END{print sum}' temp_uniq_mean_3.txt > mean_uniq.txt
	paste mean_uniq.txt seqcov.txt | awk '{print(($1/$2)* 100)}' > species_percent_coverage.txt
	rm *temp* *mean_uniq*
################################################################################################################
	cd $proj_dir/metagenome/results/species_level
	for i in {mean,unique_sequences,quantification_accuracy}; do
		awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_species/rankedlineage_subhits.txt species_taxa_${i}.txt | \
		awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' > species_inter_${i}.txt
	done
	rm *_taxa_*
################################################################################################################
#Repopulate NA in species column
	for i in $(ls *_inter_mean.txt);do
		awk -F '\t' 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}}{ print $(f["taxname"])}' $i | awk -F '_' '{print $1, $2}' | sed '1s/taxname/species/g' > ${i%_inter_mean*}_column.txt
	done
	species=$(awk -F '\t' '{for (i=1; i<=NF; i++) if ($i == "species") print i }' species_inter_mean.txt)
	awk -F '\t' -v var="$species" ' FNR==NR {a[NR]=$1;next}{$var=a[FNR]}1' OFS='\t' species_column.txt species_inter_mean.txt > species_taxainfo_mean.txt
	awk -F '\t' -v var="$species" ' FNR==NR {a[NR]=$1;next}{$var=a[FNR]}1' OFS='\t' species_column.txt species_inter_unique_sequences.txt > species_taxainfo_unique_sequences.txt
	awk -F '\t' -v var="$species" ' FNR==NR {a[NR]=$1;next}{$var=a[FNR]}1' OFS='\t' species_column.txt species_inter_quantification_accuracy.txt > species_taxainfo_quantification_accuracy.txt

	rm species_inter_mean.txt && rm species_inter_unique_sequences.txt && rm species_inter_quantification_accuracy.txt

################################################################################################################
#Species-level visualizations
	cd $proj_dir/metagenome/results/species_level
	species_level_input=species_taxainfo_mean.txt
	percent_thresh=5
	Rscript $tool_dir/Rscripts/species_level_corr.R $species_level_input $percent_thresh &>/dev/null
fi
################################################################################################################
#Genus-level sighit identification
if [ "$genus_level" == "TRUE" ]; then
cd $proj_dir/metagenome/sighits
mkdir sighits_genus
cd $proj_dir/metagenome/results
mkdir genus_level
cp seqcov.txt ./genus_level

	echo -e "${YELLOW}------------------------------------------------------------------------------ \n \n Qmatey is performing genus-level analysis \n \n------------------------------------------------------------------------------"
	cd $proj_dir/metagenome/alignment
	for i in $(ls *_haplotig_nd.megablast);do
		awk '$6>=90' $i | awk '$7>=90' | awk 'gsub(" ","_",$0)' > ../sighits/sighits_genus/${i%_haplotig*}_filter.txt
		awk 'NR==FNR {h[$1] = $2; next} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,h[$1]}' ../haplotig/${i%_haplotig*}_normalized.txt ../sighits/sighits_genus/${i%_haplotig*}_filter.txt | 
		awk '{print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' | awk 'gsub(" ","\t",$0)' > ../sighits/sighits_genus/${i%_haplotig*}_sighits.txt
	done
	rm ../sighits/sighits_genus/*_filter.txt 
################################################################################################################
#Extracts non-redundant hits for unique alignment
	cd $proj_dir/metagenome/sighits/sighits_genus/
	for i in $(ls *_sighits.txt);do
		awk -F '\t' '{a[$2]++;b[$2]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i > ${i%_sighits*}_genus_unique_reads.txt
	done
################################################################################################################
#Extracts duplicate hits for multi-alignment 
	cd $proj_dir/metagenome/sighits/sighits_genus
	for i in $(ls *_sighits.txt);do
		awk -F'|' 'FNR==NR{a[$1,$2]=1; next}  !a[$1,$2]' ${i%_sighits.txt}_genus_unique_reads.txt $i > ${i%_sighits*}_dup.txt
	done
################################################################################################################
#Appends taxa informations to the duplicate sighits file
	cd $proj_dir/metagenome/sighits/sighits_genus
	for i in $(ls *_dup.txt);do
		awk -F '\t' '{print $11}' OFS=';' $i > ${i%_dup*}_taxids_dup_inter.txt
	done

	for i in $(ls *_dup_inter.txt);do
		awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_dup_inter*}_taxids_dup.txt
	done

	rm *_taxids_dup_inter.txt
################################################################################################################
#Reformats taxonomic information to provide genus-level taxa in the species column for duplicate reads
	cd $proj_dir/metagenome/sighits/sighits_genus
	for i in $(ls *_taxids_dup.txt);do
		awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' /home/brandon/Desktop/Qmatey/tools/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_dup*}_dup_inter.txt 
	done

	for i in $(ls *_dup_inter.txt);do
		awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_dup_inter*}_species_taxid.txt
	done

	rm *_taxids_dup.txt

	for i in $(ls *_species_taxid.txt);do
		awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_species_taxid*}_species_column.txt
	done

	for i in $(ls *_species_column.txt);do
		paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column*}_species_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_species_taxa.txt
	done

	for i in $(ls *_dup.txt);do
		paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' OFS='\t' ${i%*_dup.txt}_species_taxa.txt) > ${i%_dup*}_genus_duplicates.txt
	done

	rm *_species_taxid.txt && rm *_dup_inter.txt && rm *_dup.txt && rm *_species_column.txt && rm *_species_taxa.txt
################################################################################################################
#Genus-level clustering
	for i in $(ls *_genus_duplicates.txt);do
		awk -F '\t' '{print $1, $2"~"$15, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_genus_duplicates*}_genus_inter.txt
	done

	for i in $(ls *_genus_inter.txt);do
		awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_genus_inter*}_genus_inter2.txt
	done

	for i in $(ls *_genus_inter2.txt);do
		awk -F '\t' '{dups[$1]++} END {for (num in dups) {print num}}' $i | sort -k1,1  > ${i%_genus_inter2*}_duplicate_count.txt 
	done

	for i in $(ls *_duplicate_count.txt);do
		awk -F '~' '{a[$1]++;b[$1]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i | sort -k1,1 > ${i%_duplicate_count*}_multialign_genus_reads.txt
	done

	for i in $(ls *_multialign_genus_reads.txt);do
		awk -F '\t'  'FNR==NR {a[$1]; next}; $1 in a' $i ${i%_multialign_genus_reads*}_genus_inter2.txt | sort -u -k1,1 | awk 'gsub("~","\t",$0)'| awk -F '\t' '{print $3, $1, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' OFS='\t' > ${i%_multialign_genus_reads*}_genus_OTU.txt
	done

	rm *_genus_inter.txt && rm *_genus_inter2.txt && rm *_duplicate_count.txt && rm *_multialign_genus_reads.txt && rm *_genus_duplicates.txt

	for i in $(ls *_genus_OTU.txt);do
		cat $i ${i%_genus_OTU*}_genus_unique_reads.txt > ${i%_genus_OTU*}_complete_genus_reads.txt
	done 

	for i in $(ls *_complete_genus_reads.txt);do
		awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' $i > ${i%_complete_genus_reads*}_sighits_temp.txt
		echo $'abundance\tqseqid\tsseqid\tlength\tmismatch\tevalue\tpident\tqcovs\tqseq\tsseq\tstaxids\tstitle' | \
		cat - ${i%_complete_genus_reads*}_sighits_temp.txt > ${i%_complete_genus_reads*}_sighits.txt
	done

	rm *_complete_genus_reads.txt && rm *_sighits_temp.txt && rm *_genus_OTU.txt && rm *_unique_reads.txt
#################################################################################################################
#Combine all taxids for all files/individuals and perform single search against new_taxdump.
	cd $proj_dir/metagenome/sighits/sighits_genus
	find . -type f -name '*_sighits.txt' -exec cat {} + > sighits.txt
	awk '{print $11}' sighits.txt | awk '{gsub(";","\n"); print}' | sort -u -n | sed -e '1s/staxids/tax_id/' > taxids_sighits.txt && rm sighits.txt
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}'  $tool_dir/rankedlineage_edited.dmp taxids_sighits.txt | \
	awk '{gsub(/ /,"_"); print }' > rankedlineage_subhits.txt 

	rm taxids_sighits.txt
#################################################################################################################
#Install datamash
#Generate file with mean, number of unique reads per taxID, and standard error
#Now, perform merge subsetted new_taxdump with each file while retaining only taxids in each file
	cd $proj_dir/metagenome/sighits/sighits_genus
	awk '{print $1}' rankedlineage_subhits.txt > ../../results/genus_level/genus_taxa_mean.txt
	awk '{print $1}' rankedlineage_subhits.txt > ../../results/genus_level/genus_taxa_unique_sequences.txt
	awk '{print $1}' rankedlineage_subhits.txt > ../../results/genus_level/genus_taxa_quantification_accuracy.txt 
	for i in $(ls *_sighits.txt); do
		cut -f 1,11 $i | awk '{print $2,"\t",$1}' | datamash --header-in --sort --group 1 mean 2 sstdev 2 count 2 | \
		awk '{ print $1,"\t",$2,"\t",$4,"\t",((($3/sqrt($4))/$2)*100) }' > stats1.txt
		echo $'tax_id\tmean\tuniq_reads\tstderr' | cat - stats1.txt > stats2.txt
		id=${i%_sighits*}_mean && awk -v id=$id '{gsub(/mean/,id); print }' stats2.txt | awk '{print $1,"\t",$2}' > holdmean.txt
		awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdmean.txt ../../results/genus_level/genus_taxa_mean.txt > holdmean2.txt && cat holdmean2.txt > ../../results/genus_level/genus_taxa_mean.txt
		id=${i%_sighits*}_uniq_reads && awk -v id=$id '{gsub(/uniq_reads/,id); print }' stats2.txt | awk '{print $1,"\t",$3}' > holduniq_reads.txt
		awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holduniq_reads.txt ../../results/genus_level/genus_taxa_unique_sequences.txt > holduniq_reads2.txt && cat holduniq_reads2.txt > ../../results/genus_level/genus_taxa_unique_sequences.txt
		id=${i%_sighits*}_stderr && awk -v id=$id '{gsub(/stderr/,id); print }' stats2.txt | awk '{print $1,"\t",$4}' > holdstderr.txt
		awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdstderr.txt ../../results/genus_level/genus_taxa_quantification_accuracy.txt > holdstderr2.txt && cat holdstderr2.txt > ../../results/genus_level/genus_taxa_quantification_accuracy.txt
		awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $1,"\t",$2,"\t",$3,"\t",$4,"\t", a[$1]}'  rankedlineage_subhits.txt stats2.txt > stats3.txt
		awk '{print $1,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14 }' stats3.txt | awk '{gsub(/ /,"\t"); print }' > ${i%_sighits*}_taxastats.txt
		rm *stats1* *stats2* *stats3* *hold*
	done
#################################################################################################################
#Calculates percent coverage of sighits from unmatched microbiome input files
	cd $proj_dir/metagenome/results/genus_level
	i="_mean$"
	awk -vp="$i" 'NR==1{for(i=1; i<=NF; i++) if ($i~p) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' genus_taxa_mean.txt > temp_mean.txt
	i="uniq_reads$"
	awk -vp="$i" 'NR==1{for(i=1; i<=NF; i++) if ($i~p) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' genus_taxa_unique_sequences.txt > temp_uniq.txt
	paste temp_mean.txt temp_uniq.txt | awk '/^[0-9]/ {for(i=1; i<=NF/2; i++) {s=s OFS $i*$(NF/2+i); }sub(/^ /,x,s);$0=s; s=""} !/[0-9]/{$0=$1;}1' > temp_uniq_mean.txt
	tail -n +2 temp_uniq_mean.txt > temp_uniq_mean_2.txt
	awk '{for (i=1; i<=NF; i++) sum[i]+=$i;}; END{for (i in sum) print sum[i];}' temp_uniq_mean_2.txt > temp_uniq_mean_3.txt
	awk '{sum+=$1}END{print sum}' temp_uniq_mean_3.txt > mean_uniq.txt
	paste mean_uniq.txt seqcov.txt | awk '{print(($1/$2)* 100)}' > genus_percent_coverage.txt
	rm *temp* *mean_uniq*
################################################################################################################
	cd $proj_dir/metagenome/results/genus_level
	for i in {mean,unique_sequences,quantification_accuracy}; do
		awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_genus/rankedlineage_subhits.txt genus_taxa_${i}.txt | \
		awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' > genus_taxainfo_${i}.txt
	done
	rm *_taxa_*
################################################################################################################
#Genus-level visualizations
	genus_level_input=genus_taxainfo_mean.txt
	percent_thresh=5
	Rscript $tool_dir/Rscripts/genus_level_corr.R $genus_level_input $percent_thresh &>/dev/null
fi
################################################################################################################
#Family-Level sighit identification
if [ "$family_level" == "TRUE" ];then
cd $proj_dir/metagenome/sighits
mkdir sighits_family
cd $proj_dir/metagenome/results
mkdir family_level
cp seqcov.txt ./family_level

	echo -e "${YELLOW}------------------------------------------------------------------------------ \n \n Qmatey is performing family-level analysis \n \n------------------------------------------------------------------------------"
	cd $proj_dir/metagenome/alignment
	for i in $(ls *_haplotig_nd.megablast);do
		awk '$6>=90' $i | awk '$7>=90' | awk 'gsub(" ","_",$0)' > ../sighits/sighits_family/${i%_haplotig*}_filter.txt
		awk 'NR==FNR {h[$1] = $2; next} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,h[$1]}' ../haplotig/${i%_haplotig*}_normalized.txt ../sighits/sighits_family/${i%_haplotig*}_filter.txt | 
		awk '{print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' | awk 'gsub(" ","\t",$0)' > ../sighits/sighits_family/${i%_haplotig*}_sighits.txt
	done
	rm ../sighits/sighits_family/*_filter.txt 
################################################################################################################
#Extracts non-redundant hits for unique alignment
	cd $proj_dir/metagenome/sighits/sighits_family/
	for i in $(ls *_sighits.txt);do
		awk -F '\t' '{a[$2]++;b[$2]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i > ${i%_sighits*}_family_unique_reads.txt
	done
################################################################################################################
#Extracts duplicate hits for multi-alignment 
	cd $proj_dir/metagenome/sighits/sighits_family
	for i in $(ls *_sighits.txt);do
			awk -F'|' 'FNR==NR{a[$1,$2]=1; next}  !a[$1,$2]' ${i%_sighits.txt}_family_unique_reads.txt $i > ${i%_sighits*}_dup.txt
	done
################################################################################################################
#Appends taxa informations to the duplicate sighits file
	cd $proj_dir/metagenome/sighits/sighits_family
	for i in $(ls *_dup.txt);do
		awk -F '\t' '{print $11}' OFS=';' $i > ${i%_dup*}_taxids_dup_inter.txt
	done

	for i in $(ls *_dup_inter.txt);do
		awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_dup_inter*}_taxids_dup.txt
	done

	rm *_taxids_dup_inter.txt
################################################################################################################
#Reformats taxonomic information to provide family-level taxa in the species column for duplicate reads
	cd $proj_dir/metagenome/sighits/sighits_family
	for i in $(ls *_taxids_dup.txt);do
		awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' /home/brandon/Desktop/Qmatey/tools/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_dup*}_dup_inter.txt 
	done

	for i in $(ls *_dup_inter.txt);do
		awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_dup_inter*}_species_taxid.txt
	done

	rm *_taxids_dup.txt

	for i in $(ls *_species_taxid.txt);do
		awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_species_taxid*}_species_column.txt
	done

	for i in $(ls *_species_column.txt);do
		paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column*}_species_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_species_taxa.txt
	done

	for i in $(ls *_dup.txt);do
		paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' OFS='\t' ${i%*_dup.txt}_species_taxa.txt) > ${i%_dup*}_family_duplicates.txt
	done

	rm *_species_taxid.txt && rm *_dup_inter.txt && rm *_dup.txt && rm *_species_column.txt && rm *_species_taxa.txt
################################################################################################################
#Family-level clustering
	for i in $(ls *_family_duplicates.txt);do
		awk -F '\t' '{print $1, $2"~"$16, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_family_duplicates*}_family_inter.txt
	done

	for i in $(ls *_family_inter.txt);do
		awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_family_inter*}_family_inter2.txt
	done

	for i in $(ls *_family_inter2.txt);do
		awk -F '\t' '{dups[$1]++} END {for (num in dups) {print num}}' $i | sort -k1,1  > ${i%_family_inter2*}_duplicate_count.txt 
	done

	for i in $(ls *_duplicate_count.txt);do
		awk -F '~' '{a[$1]++;b[$1]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i | sort -k1,1 > ${i%_duplicate_count*}_multialign_family_reads.txt
	done

	for i in $(ls *_multialign_family_reads.txt);do
		awk -F '\t'  'FNR==NR {a[$1]; next}; $1 in a' $i ${i%_multialign_family_reads*}_family_inter2.txt | sort -u -k1,1 | awk 'gsub("~","\t",$0)'| awk -F '\t' '{print $3, $1, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' OFS='\t' > ${i%_multialign_family_reads*}_family_OTU.txt
	done

	rm *_family_inter.txt && rm *_family_inter2.txt && rm *_duplicate_count.txt && rm *_multialign_family_reads.txt && rm *_family_duplicates.txt
	for i in $(ls *_family_OTU.txt);do
		cat $i ${i%_family_OTU*}_family_unique_reads.txt > ${i%_family_OTU*}_complete_family_reads.txt
	done 

	for i in $(ls *_complete_family_reads.txt);do
		awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' $i > ${i%_complete_family_reads*}_sighits_temp.txt
		echo $'abundance\tqseqid\tsseqid\tlength\tmismatch\tevalue\tpident\tqcovs\tqseq\tsseq\tstaxids\tstitle' | \
		cat - ${i%_complete_family_reads*}_sighits_temp.txt > ${i%_complete_family_reads*}_sighits.txt
	done

	rm *_complete_family_reads.txt && rm *_sighits_temp.txt && rm *_family_OTU.txt && rm *_unique_reads.txt
#################################################################################################################
	#Combine all taxids for all files/individuals and perform single search against new_taxdump.
	cd $proj_dir/metagenome/sighits/sighits_family
	find . -type f -name '*_sighits.txt' -exec cat {} + > sighits.txt
	awk '{print $11}' sighits.txt | awk '{gsub(";","\n"); print}' | sort -u -n | sed -e '1s/staxids/tax_id/' > taxids_sighits.txt && rm sighits.txt
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}'  $tool_dir/rankedlineage_edited.dmp taxids_sighits.txt | \
	awk '{gsub(/ /,"_"); print }' > rankedlineage_subhits.txt 

	rm taxids_sighits.txt
#################################################################################################################
#Install datamash
#Generate file with mean, number of unique reads per taxID, and standard error
#Now, perform merge subsetted new_taxdump with each file while retaining only taxids in each file
	cd $proj_dir/metagenome/sighits/sighits_family
	awk '{print $1}' rankedlineage_subhits.txt > ../../results/family_level/family_taxa_mean.txt
	awk '{print $1}' rankedlineage_subhits.txt > ../../results/family_level/family_taxa_unique_sequences.txt
	awk '{print $1}' rankedlineage_subhits.txt > ../../results/family_level/family_taxa_quantification_accuracy.txt 
	for i in $(ls *_sighits.txt); do
		cut -f 1,11 $i | awk '{print $2,"\t",$1}' | datamash --header-in --sort --group 1 mean 2 sstdev 2 count 2 | \
		awk '{ print $1,"\t",$2,"\t",$4,"\t",((($3/sqrt($4))/$2)*100) }' > stats1.txt
		echo $'tax_id\tmean\tuniq_reads\tstderr' | cat - stats1.txt > stats2.txt
		id=${i%_sighits*}_mean && awk -v id=$id '{gsub(/mean/,id); print }' stats2.txt | awk '{print $1,"\t",$2}' > holdmean.txt
		awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdmean.txt ../../results/family_level/family_taxa_mean.txt > holdmean2.txt && cat holdmean2.txt > ../../results/family_level/family_taxa_mean.txt
		id=${i%_sighits*}_uniq_reads && awk -v id=$id '{gsub(/uniq_reads/,id); print }' stats2.txt | awk '{print $1,"\t",$3}' > holduniq_reads.txt
		awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holduniq_reads.txt ../../results/family_level/family_taxa_unique_sequences.txt > holduniq_reads2.txt && cat holduniq_reads2.txt > ../../results/family_level/family_taxa_unique_sequences.txt
		id=${i%_sighits*}_stderr && awk -v id=$id '{gsub(/stderr/,id); print }' stats2.txt | awk '{print $1,"\t",$4}' > holdstderr.txt
		awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdstderr.txt ../../results/family_level/family_taxa_quantification_accuracy.txt > holdstderr2.txt && cat holdstderr2.txt > ../../results/family_level/family_taxa_quantification_accuracy.txt
		awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $1,"\t",$2,"\t",$3,"\t",$4,"\t", a[$1]}'  rankedlineage_subhits.txt stats2.txt > stats3.txt
		awk '{print $1,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14 }' stats3.txt | awk '{gsub(/ /,"\t"); print }' > ${i%_sighits*}_taxastats.txt
		rm *stats1* *stats2* *stats3* *hold*
	done
#################################################################################################################
#Calculates percent coverage of sighits from unmatched microbiome input files
	cd $proj_dir/metagenome/results/family_level
	i="_mean$"
	awk -vp="$i" 'NR==1{for(i=1; i<=NF; i++) if ($i~p) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' family_taxa_mean.txt > temp_mean.txt
	i="uniq_reads$"
	awk -vp="$i" 'NR==1{for(i=1; i<=NF; i++) if ($i~p) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' family_taxa_unique_sequences.txt > temp_uniq.txt
	paste temp_mean.txt temp_uniq.txt | awk '/^[0-9]/ {for(i=1; i<=NF/2; i++) {s=s OFS $i*$(NF/2+i); }sub(/^ /,x,s);$0=s; s=""} !/[0-9]/{$0=$1;}1' > temp_uniq_mean.txt
	tail -n +2 temp_uniq_mean.txt > temp_uniq_mean_2.txt
	awk '{for (i=1; i<=NF; i++) sum[i]+=$i;}; END{for (i in sum) print sum[i];}' temp_uniq_mean_2.txt > temp_uniq_mean_3.txt
	awk '{sum+=$1}END{print sum}' temp_uniq_mean_3.txt > mean_uniq.txt
	paste mean_uniq.txt seqcov.txt | awk '{print(($1/$2)* 100)}' > family_percent_coverage.txt
	rm *temp* *mean_uniq*
################################################################################################################
	cd $proj_dir/metagenome/results/family_level
	for i in {mean,unique_sequences,quantification_accuracy}; do
		awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_family/rankedlineage_subhits.txt family_taxa_${i}.txt | \
		awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' > family_taxainfo_${i}.txt
	done
	rm *_taxa_*
################################################################################################################
#Family-level visualizations
	family_level_input=family_taxainfo_mean.txt
	percent_thresh=5
	Rscript $tool_dir/Rscripts/family_level_corr.R $family_level_input $percent_thresh &>/dev/null
fi




