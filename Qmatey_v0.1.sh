#!/bin/bash
[ $# -ge 1 -a -f "$1" ] && input="$1" || { echo -e "$1 \e[31mnot found, exiting"; exit 1; }
. $input
BLUE='\e[38;5;210m'
WHITE='\e[97m'
ORANGE='\e[38;5;210m'
proj_dir=${input%/*}
dependency () {
if java -version; then
	:
else
	echo "java not detected, try 'sudo apt install default-jre'"
	exit 1
fi



}
dependency &> $proj_dir/log.out

if [ -z "$threads" ]; then
	threads=$(nproc --all)
fi
if  [[ "$threads" -gt 1 ]]; then
	loopthread=2
	N=$(($threads/2))
else
	N=1 && loopthread=$threads
fi

cd $tool_dir
tar -xf rankedlineage.tar.xz

##################################################################################################################
#Check for existence of directories specified in config file
[ -d "$tool_dir" ] || { echo -e "$tool_dir \e[31mnot found, exiting"; exit 1; }
[ -d "$input_dir" ] || { echo -e "$input_dir \e[31mnot found, exiting"; exit 1; }
[ -d "$proj_dir" ] || { echo -e "$proj_dir \e[31mnot found, exiting"; exit 1; }
##################################################################################################################
#Create all necessary subdirectories for downstream processing
cd $proj_dir
mkdir metagenome || { echo -e "\e[31mproject directory can't contain folder named metagenome, exiting"; exit 1; }
cd $proj_dir/metagenome
mkdir haplotig
mkdir alignment
mkdir sighits
mkdir results
##################################################################################################################
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is installing R dependencies \n\e[97m########################################################\n"
if R --version; then
	:
else
	echo -e "\e[31mR not detected,exiting"
	exit 1
fi 
Rscript $tool_dir/Rscripts/R_packages.R &>/dev/null
#################################################################################################################
ref_norm (){
[ -d "$norm_ref_dir" ] || { echo -e "$norm_ref_dir \e[31mnot found, exiting"; exit 1; }
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is Preforming Normalization \n\e[97m########################################################\n"
cd $norm_ref_dir
echo -e "${YELLOW}- indexing normalization reference genome${WHITE}"
cat *.fa* > master_ref.fasta
$tool_dir/bwa-0.7.17/bwa index -a bwtsw master_ref.fasta
java -jar $tool_dir/picard-2.18.27.jar CreateSequenceDictionary REFERENCE=master_ref.fasta    OUTPUT=master_ref.dict
$tool_dir/samtools-1.9/samtools faidx master_ref.fasta 


cd $input_dir
echo -e "${YELLOW}- aligning input sequences to normalization reference genome${WHITE}"
for i in $(ls *.fa*); do
	$tool_dir/bwa-0.7.17/bwa mem -t $loopthread $norm_ref_dir/master_ref.fasta $i > $proj_dir/metagenome/${i%.fa*}.sam && \
	java -XX:ParallelGCThreads=$loopthread -jar $tool_dir/picard-2.18.27.jar SortSam I= $proj_dir/metagenome/${i%.fa*}.sam O= $proj_dir/metagenome/${i%.fa*}.bam SORT_ORDER=coordinate && \
	rm $proj_dir/metagenome/${i%.fa*}.sam
done

cd $proj_dir/metagenome
echo -e "${YELLOW}- extracting sequence reads that do not match the normalization reference genome (i.e. potential metagenomic data)${WHITE}"
for i in $(ls *.bam); do
	$tool_dir/samtools-1.9/samtools view -b -f4 $i > ${i%.*}_metagenome.bam && \
	$tool_dir/samtools-1.9/samtools bam2fq ${i%.*}_metagenome.bam | gzip > ${i%.*}_metagenome.fastq.gz && \
	rm ${i%.*}_metagenome.bam
done

echo -e "${YELLOW}- calculating metagenomic coverage"

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

echo -e "${YELLOW}- calculating a normalization factor"
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

echo -e "${YELLOW}- reformatting metagenomic reads into .fasta format" 
cd $proj_dir/metagenome
for i in $(ls *_metagenome.fastq.gz); do
	gunzip -k $i && mv ${i%_meta*}_metagenome.fastq ./haplotig/${i%_meta*}_metagenome.fastq && cd ./haplotig
	awk 'NR%2==0' ${i%_meta*}_metagenome.fastq | awk 'NR%2==1' | sort | uniq -c | tr -s ' ' '\t' | awk '$0=">hseqid_"NR$0' | awk '$2>"0"' > ${i%_meta*}_cov_haplotig.txt
	awk '{print $1"\t"$2}' ${i%_meta*}_cov_haplotig.txt | awk  'gsub(">", "", $0)' > ${i%_meta*}_haplocov.txt
	awk '{print $1"\n"$3}' ${i%_meta*}_cov_haplotig.txt > ${i%_meta*}_haplotig.fasta
	rm ${i%_meta*}_metagenome.fastq && rm ${i%_meta*}_cov_haplotig.txt && cd $proj_dir/metagenome
done
cd $proj_dir/metagenome/haplotig
for i in $(ls *_haplocov.txt);do 
	id=${i%_haplocov*}
	normfactor=$(awk -v id=$id '$1 == id {print $2}' ../coverage_normalization_factor.txt)
	awk -v normfactor=$normfactor '{print $1,$2 * normfactor}' $i | awk 'gsub(" ", "\t", $0)' > ${i%_haplocov*}_normalized.txt
done
}
if [ "$normalization" == "TRUE" ]; then
	time ref_norm 2>> $proj_dir/log.out
fi

bypass_norm (){
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is Bypassing Normalization \n\e[97m########################################################\n"
cd $input_dir
if [ "$data_type" == "FASTA" ]; then
	echo -e "${YELLOW}- calculating metagenomic coverage" 
	for i in $(ls *.f*); do 
		grep -o "^>" $i | wc -l > ../metagenome/results/${i%.*}_seq.txt
	done
	cd $proj_dir/metagenome/results
	for i in $(ls *_seq.txt); do
		cat $i >> tempcov.txt | awk '{sum+=$1}END{print sum}' tempcov.txt > seqcov.txt
	done
	rm *_seq.txt & rm *tempcov.txt
	cd $input_dir
	for i in $(ls *.f*);do
		awk < $i '/^>/ { print $0, 1}' | awk  'gsub(">", "", $0)' > ../metagenome/haplotig/${i%.*}_haplocov.txt
	done
	cd $proj_dir/metagenome/haplotig
	for i in $(ls *_haplocov.txt);do
		id=${i%_haplocov*}
		normfactor=1
		awk -v normfactor=$normfactor '{print $1,$2 * normfactor}' $i | awk 'gsub(" ", "\t", $0)' > ${i%_haplocov*}_normalized.txt
	done
else
	echo -e "${YELLOW}- calculatng metagenomic coverage"
	for i in $(ls *.gz);do
		gunzip $i

	done
	for i in $(ls *.f*);do 
		awk '{if(NR%4==2) print $0}' $i > ../metagenome/results/${i%.*}_seq.txt
	done
	cd $proj_dir/metagenome/results
	for i in $(ls *_seq.txt);do
		wc -l $i > ${i%meta*}_seqcov.txt
		find . -type f -name '*_seqcov.txt' -exec cat {} + > tempcov.txt
		awk '{sum+=$1}END{print sum}' tempcov.txt > seqcov.txt
	done
	rm *_seq.txt && rm *tempcov.txt && rm *_seqcov.txt
	cd $input_dir
	echo -e "${YELLOW}- converting files into .fasta format"
	for i in $(ls *.fastq);do
		awk 'NR%2==0' $i | awk 'NR%2==1' | sort | uniq -c | tr -s ' ' '\t' | awk '$0=">hseqid_"NR$0' | awk '$2>"0"' > ../metagenome/haplotig/${i%.fastq*}_cov_haplotig.txt
	done
	cd $proj_dir/metagenome/haplotig
	for i in $(ls *_cov_haplotig.txt);do
		awk '{print $1"\t"$2}' $i | awk  'gsub(">", "", $0)' > ${i%_cov_haplotig*}_haplocov.txt
	done
	for i in $(ls *_cov_haplotig.txt);do
		awk '{print $1"\n"$3}' $i > ${i%_cov_haplotig*}_haplotig.fasta
	done
	rm *_cov_haplotig.txt
	cd $proj_dir/metagenome/haplotig
	for i in $(ls *_haplocov.txt);do
		id=${i%_haplocov*}
		normfactor=1
		awk -v normfactor=$normfactor '{print $1,$2 * normfactor}' $i | awk 'gsub(" ", "\t", $0)' > ${i%_haplocov*}_normalized.txt
	done
fi
}
if [ "$normalization" == "FALSE" ]; then
	time bypass_norm 2>> $proj_dir/log.out
fi

blast(){
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is Prefroming BLAST \n\e[97m########################################################\n"
cd $proj_dir/metagenome/haplotig
if [ "$blast_location" == "LOCAL" ]; then
echo -e "${YELLOW}- preforming a local BLAST"
	for i in $(ls *_haplotig.fasta);do
		$tool_dir/ncbi-blast-2.8.1+/bin/blastn -task megablast -query $i -db $local_db_dir -num_threads $threads -evalue 1e-12 -perc_identity 95 -outfmt "6 qseqid sseqid length mismatch evalue pident qcovs qseq sseq staxids stitle" -out ../alignment/${i%_haplotig*}_haplotig.megablast
	done
	for i in $(ls *_haplotig.megablast);do
		sort -u $i > ${i%_haplotig*}_haplotig_nd.megablast && rm *_haplotig.megablast
	done
fi
if [ "$blast_location" == "REMOTE" ]; then 
echo -e "${YELLOW}- preforming a remote BLAST"
	for i in $(ls *_haplotig.fasta);do
		$tool_dir/ncbi-blast-2.8.1+/bin/blastn -task megablast -query $i -db $remote_db_dir -evalue 1e-12 -perc_identity 95 -outfmt \
		"6 qseqid sseqid length mismatch evalue pident qcovs qseq sseq staxids stitle" \
		-out ../alignment/${i%_haplotig*}_haplotig.megablast -remote
	done
	for i in $(ls *_haplotig.megablast);do
		sort -u $i > ${i%_haplotig*}_haplotig_nd.megablast && rm *_haplotig.megablast
	done
fi
}
time blast 2>> $proj_dir/log.out

##################################################################################################################
#Removes duplicate rows and vector contamination from *_haplotig.megablast
awk '{gsub(/\t\t/,"\tNA\t"); print}' $tool_dir/rankedlineage.dmp | awk '{gsub(/[|]/,""); print}' | awk '{gsub(/\t\t/,"\t"); print}' > $tool_dir/rankedlineage_tabdelimited.dmp
echo $'tax_id\ttaxname\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tsuperkingdom\t' | \
cat - $tool_dir/rankedlineage_tabdelimited.dmp > $tool_dir/rankedlineage_edited.dmp
rm $tool_dir/rankedlineage_tabdelimited.dmp

##################################################################################################################
strain_filter(){
cd $proj_dir/metagenome/sighits
mkdir sighits_strain
cd $proj_dir/metagenome/results
mkdir strain_level
cp seqcov.txt ./strain_level
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is Preforming Strain-Level Classification \n\e[97m########################################################\n"
cd $proj_dir/metagenome/alignment
echo -e "${YELLOW}- preforming exact-matching algorithm"
for i in $(ls *_haplotig_nd.megablast);do
	awk '$6>=100' $i | awk '$7>=100' | awk 'gsub(" ","_",$0)' > ../sighits/sighits_strain/${i%_haplotig*}_filter.txt
	awk 'NR==FNR {h[$1] = $2; next} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,h[$1]}' ../haplotig/${i%_haplotig*}_normalized.txt ../sighits/sighits_strain/${i%_haplotig*}_filter.txt | 
	awk '{print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' | awk 'gsub(" ","\t",$0)' > ../sighits/sighits_strain/${i%_haplotig*}_sighits.txt
done
rm ../sighits/sighits_strain/*_filter.txt 
cd $proj_dir/metagenome/sighits/sighits_strain
for i in $(ls *_sighits.txt);do
	awk '{a[$2]++;b[$2]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i > ${i%.txt}_nr_temp.txt
	echo $'abundance\tqseqid\tsseqid\tlength\tmismatch\tevalue\tpident\tqcovs\tqseq\tsseq\tstaxids\tstitle' | \
	cat - ${i%.txt}_nr_temp.txt > ${i%.txt}_nr.txt
done
rm *_sighits.txt *_nr_temp.txt
for i in $(ls *_sighits_nr.txt);do
	awk -F '\t' 'NR>1{print $0}' $i > ${i%_sighits_nr*}_sighits_nr_temp.txt
done

for i in $(ls *_sighits_nr_temp.txt);do
	awk -F '\t' '{print $11}' OFS=';' $i > ${i%_sighits_nr_temp*}_taxids_uniq_inter.txt
done

for i in $(ls *_uniq_inter.txt);do
	awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_uniq_inter*}_taxids_uniq.txt
done

rm *_uniq_inter.txt
cd $proj_dir/metagenome/sighits/sighits_strain
for i in $(ls *_taxids_uniq.txt);do
	(
	awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' $tool_dir/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_uniq*}_uniq_inter.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi 
done
wait
for i in $(ls *_uniq_inter.txt);do
	(
	awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_uniq_inter*}_species_taxid.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
rm *_taxids_uniq.txt
for i in $(ls *_species_taxid.txt);do
	awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_species_taxid*}_species_column.txt
done
for i in $(ls *_species_column.txt);do
	paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column*}_species_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_species_taxa.txt
done
for i in $(ls *_sighits_nr_temp.txt);do
	paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%*_sighits_nr_temp*}_species_taxa.txt) > ${i%_sighits_nr_temp*}_sighits_uncultured.txt
done
rm *_species_taxid.txt && rm *_uniq_inter.txt && rm *_species_column.txt && rm *_species_taxa.txt && rm *_sighits_nr.txt && rm *_sighits_nr_temp.txt
for i in $(ls *_sighits_uncultured.txt);do
	awk -F '\t' '!/Uncultured/' $i > ${i%_sighits_uncultured*}_sighits_nr.txt
done
rm *_sighits_uncultured.txt
}
if [ "$strain_level" == "TRUE" ]; then
	time strain_filter 2>> $proj_dir/log.out
fi

strain_refseq_filter(){
cd $proj_dir/metagenome/sighits/sighits/strain
mkdir strain_gene_annotation
for i in $(ls *_sighits_nr.txt);do
	cp $i $proj_dir/metagenome/sighits/sighits_strain/strain_gene_annotation
done
cd $proj_dir/metagenome/sighits/sighis_strain/strain_gene_annotation
for i in $(ls *_sighits_nr.txt);do
	awk -F '\t' -v var="${i%*_sighits_nr.txt}" '{print $1"\t"var"_"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20}' $i > ${i%*_sighits_nr.txt}_refseq_sighits.txt
done
rm *_sighits_nr.txt
for i in $(ls *_refseq_sighits.txt);do
	awk 'FNR==NR{ h=$0 }NR>1{ print (!a[$11]++? h ORS $0 : $0) >> $11".txt" }' $i
done
for i in $(ls *.txt);do
	sed -i '/^$/d' $i 
done
for i in $(ls *_refseq_sighits.txt);do
	mv $i $proj_dir/metagenome/sighits/sighits_strain
done
for i in $(ls *.txt);do
	awk -F '\t' '{print ">"$2"\n"$9}' $i > ${i%.txt}.fasta
done
rm *.txt
for i in $(ls *.fasta);do
	$tool_dir/ncbi-blast-2.8.1+/bin/blastn -task megablast -query $i -db $refseq_db_dir -num_threads $threads -evalue 1e-10 -max_target_seqs 1 -outfmt "6 qseqid sseqid length mismatch evalue pident qcovs qseq sseq staxids stitle" -out ${i%*.fasta}_refseq.megablast
done
for i in $(ls *_refseq.megablast);do
	sort -u $i > ${i%*_refseq.megablast}_refseq_nr.txt
done
rm *_refseq.megablast && rm *.fasta
find . -type f -name '*_refseq_nr.txt' -exec cat {} + > refseq_hits.txt
awk -F '\t' '{print $10}' OFS=';' refseq_hits.txt > taxids_refseq_inter.txt
awk -F ';' '{print $1}' OFS='\t' taxids_refseq_inter.txt > taxids_refseq.txt
rm *_refseq_inter.txt
for i in $(ls taxids_refseq.txt);do
	(
	awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' $tool_dir/rankedlineage_edited.dmp OFS='\t' $i> ${i%taxids_refseq*}refseq_inter.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi 
done
paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' OFS='\t' refseq_hits.txt ) <(awk -F '\t' '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' refseq_inter.txt) > refseq_hits_genus.txt
rm refseq_inter.txt && rm taxids_refseq.txt && rm refseq_hits.txt
cd $proj_dir/metagenome/sighits/sighits_strain/
for i in $(ls *_refseq_sighits.txt);do
	mv $i $proj_dir/metagenome/sighits/sighits_strain/strain_gene_annotation
done
cd $proj_dir/metagenome/sighits/sighits_strain/strain_gene_annotation
find . -type f -name '*_refseq_sighits.txt' -exec cat {} + > refseq_prehits_temp.txt
awk -F '\t' '!/abundance/' refseq_prehits_temp.txt > refseq_prehits.txt && rm refseq_prehits_temp.txt
awk -F '\t' '{print $2, $14}' OFS='\t' refseq_prehits.txt > prehits.txt
awk -F '\t' '{print $1, $14}' OFS='\t' refseq_hits_genus.txt > hits.txt
awk -F '\t' 'NR==FNR{c[$1$2]++;next};c[$1$2] > 0' prehits.txt hits.txt > consistent_hits_taxids.txt
awk -F '\t' '{print $1}' consistent_hits_taxids.txt > consistent_hits.txt && rm consistent_hits_taxids.txt
for i in $(ls *_refseq_sighits.txt);do
	awk -F '\t' 'NR==FNR{c[$1]++;next};c[$2] > 0' consistent_hits.txt $i > ${i%*_refseq_sighits.txt}_refseq_filtered.txt
done
rm consistent_hits.txt && rm hits.txt && rm prehits.txt && rm refseq_prehits.txt && rm *_refseq_sighits.txt && rm *refseq_hits_genus
for i in $(ls *_refseq_filtered.txt);do
	mv $i $proj_dir/metagenome/sighits/sighits_strain/
done 
cd $proj_dir/metagenome/sighits/sighits_strain/
for i in $(ls *refseq_filtered.txt);do
	echo $'abundance\tqseqid\tsseqid\tlength\tmismatch\tevalue\tpident\tqcovs\tqseq\tsseq\tstaxids\tstitle\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tdomain' | \
	cat - $i > ${i%*_refseq_filtered.txt}_sighits_nr.txt
done
rm *_refseq_filtered.txt
if [ "$refseq_filter" == "TRUE" ] && [ "$strain_level" == "TRUE" ]; then
	time strain_refseq_filter 2>> $proj_dir/log.out
fi

strain_stats(){
echo -e "${YELLOW}- compiling taxonomic information"
cd $proj_dir/metagenome/sighits/sighits_strain
find . -type f -name '*_sighits_nr.txt' -exec cat {} + > sighits.txt
awk -F '\t' '{print $11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20}' sighits.txt > rankedlineage_subhits.txt && rm sighits.txt
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -u | awk '!/taxname/' > strain_taxa_mean_temp1.txt
echo -e 'taxname' | cat - strain_taxa_mean_temp1.txt > strain_taxa_mean_temp.txt && rm strain_taxa_mean_temp1.txt
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -u | awk '!/taxname/' > strain_taxa_unique_sequences_temp1.txt
echo -e 'taxname' | cat - strain_taxa_unique_sequences_temp1.txt > strain_taxa_unique_sequences_temp.txt && rm strain_taxa_unique_sequences_temp1.txt
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -u | awk '!/taxname/' > strain_taxa_quantification_accuracy_temp1.txt
echo -e 'taxname' | cat - strain_taxa_quantification_accuracy_temp1.txt > strain_taxa_quantification_accuracy_temp.txt && rm strain_taxa_quantification_accuracy_temp1.txt
echo -e "${YELLOW}- quantifying the strain-level taxonomy" 
strain_level=strain
for i in $(ls *_sighits_nr.txt);do
	Rscript $tool_dir/Rscripts/stats_summary.R $i $min_uniq $strain_level
	echo $'tax_id\tmean\tuniq_reads\tstderr' | cat - stats1.txt > stats2.txt 
	id=${i%_sighits*}_mean && awk -v id=$id '{gsub(/mean/,id); print }' stats2.txt | awk '{print $1,"\t",$2}' > holdmean.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdmean.txt strain_taxa_mean_temp.txt > holdmean2.txt && cat holdmean2.txt > strain_taxa_mean_temp.txt
	id=${i%_sighits*}_uniq_reads && awk -v id=$id '{gsub(/uniq_reads/,id); print }' stats2.txt | awk '{print $1,"\t",$3}' > holduniq_reads.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holduniq_reads.txt strain_taxa_unique_sequences_temp.txt > holduniq_reads2.txt && cat holduniq_reads2.txt > strain_taxa_unique_sequences_temp.txt
	id=${i%_sighits*}_stderr && awk -v id=$id '{gsub(/stderr/,id); print }' stats2.txt | awk '{print $1,"\t",$4}' > holdstderr.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdstderr.txt strain_taxa_quantification_accuracy_temp.txt > holdstderr2.txt && cat holdstderr2.txt > strain_taxa_quantification_accuracy_temp.txt
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $1,"\t",$2,"\t",$3,"\t",$4,"\t", a[$1]}'  rankedlineage_subhits.txt stats2.txt > stats3.txt
	awk '{print $1,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14 }' stats3.txt | awk '{gsub(/ /,"\t"); print }' > ${i%_sighits*}_taxastats.txt
	rm *stats1* *stats2* *stats3* *hold*
done
	
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' strain_taxa_mean_temp.txt > strain_taxa_mean_temp2.txt && rm strain_taxa_mean_temp.txt
awk '{gsub(/ /,"\t"); print}' strain_taxa_mean_temp2.txt > ../../results/strain_level/strain_taxa_mean.txt && rm strain_taxa_mean_temp2.txt
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' strain_taxa_unique_sequences_temp.txt > strain_taxa_unique_sequences_temp2.txt && rm strain_taxa_unique_sequences_temp.txt
awk '{gsub(/ /,"\t"); print}' strain_taxa_unique_sequences_temp2.txt > ../../results/strain_level/strain_taxa_unique_sequences.txt && rm strain_taxa_unique_sequences_temp2.txt
awk '{gsub(/ /,"\t"); print}' strain_taxa_quantification_accuracy_temp.txt > strain_taxa_quantification_accuracy_temp2.txt && rm strain_taxa_quantification_accuracy_temp.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/strain_level/strain_taxa_mean.txt strain_taxa_quantification_accuracy_temp2.txt > ../../results/strain_level/strain_taxa_quantification_accuracy.txt && rm strain_taxa_quantification_accuracy_temp2.txt
rm *_temp.txt
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
cd $proj_dir/metagenome/results/strain_level
for i in {mean,unique_sequences,quantification_accuracy}; do
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_strain/rankedlineage_subhits.txt strain_taxa_${i}.txt | \
	awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' > strain_taxainfo_${i}.txt
done
rm *_taxa_*
echo -e "${YELLOW}- creating strain-level visualizations"
cd $proj_dir/metagenome/results/strain_level
strain_level_mean=strain_taxainfo_mean.txt
strain_level_uniq=strain_taxainfo_unique_sequences.txt
strain_level_stderr=strain_taxainfo_quantification_accuracy.txt
percent_thresh=5
Rscript $tool_dir/Rscripts/strain_level_corr.R $strain_level_mean $percent_thresh &>/dev/null
Rscript $tool_dir/Rscripts/strain_level_boxplots.R $strain_level_mean $strain_level_uniq $strain_level_stderr $percent_thresh &>/dev/null
mkdir boxplots
mv *_files $proj_dir/metagenome/results/strain_level/boxplots
mv *.html $proj_dir/metagenome/results/strain_level/boxplots
}
if [ "$strain_level" == "TRUE" ]; then
	time strain_stats 2>> $proj_dir/log.out
fi
################################################################################################################

species_filter(){
cd $proj_dir/metagenome/sighits
mkdir sighits_species
cd $proj_dir/metagenome/results
mkdir species_level
cp seqcov.txt ./species_level
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is Preforming Species-Level Classification \n\e[97m########################################################\n"
cd $proj_dir/metagenome/alignment
echo -e "${YELLOW}- preforming exact-matching algorithm"
for i in $(ls *_haplotig_nd.megablast);do
	awk '$6>=97' $i | awk 'gsub(" ","_",$0)' > ../sighits/sighits_species/${i%_haplotig*}_filter.txt
	awk 'NR==FNR {h[$1] = $2; next} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,h[$1]}' ../haplotig/${i%_haplotig*}_normalized.txt ../sighits/sighits_species/${i%_haplotig*}_filter.txt | 
	awk '{print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' | awk 'gsub(" ","\t",$0)' > ../sighits/sighits_species/${i%_haplotig*}_sighits.txt
done
rm ../sighits/sighits_species/*_filter.txt 
cd $proj_dir/metagenome/sighits/sighits_species
for i in $(ls *_sighits.txt);do
	awk -F '\t' '{a[$2]++;b[$2]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i > ${i%_sighits*}_species_unique_reads.txt
done
echo -e "${YELLOW}- compiling species-level OTUs with multi-alignment algorithm"
cd $proj_dir/metagenome/sighits/sighits_species
for i in $(ls *_sighits.txt);do
	(
	awk -F '\t' 'FNR==NR{a[$1,$2]=1; next}  !a[$1,$2]' ${i%_sighits*}_species_unique_reads.txt $i OFS='\t' > ${i%_sighits*}_dup.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
cd $proj_dir/metagenome/sighits/sighits_species
for i in $(ls *_dup.txt);do
	awk -F '\t' '{print $11}' OFS=';' $i > ${i%_dup*}_taxids_dup_inter.txt
done
for i in $(ls *_dup_inter.txt);do
	(
	awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_dup_inter*}_taxids_dup.txt
	) & 
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
rm *_dup_inter.txt
cd $proj_dir/metagenome/sighits/sighits_species
for i in $(ls -S *_taxids_dup.txt);do
	(
	awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' $tool_dir/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_dup*}_dup_inter.txt 
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
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
	awk -F '\t' '{ if ($21!="Viruses") print $0}' $i | awk -F '\t' '!/Uncultured/' > ${i%*_species_duplicates_virome*}_species_duplicates.txt
done
rm *_species_taxid.txt && rm *_dup_inter.txt && rm *_dup.txt && rm *_species_column.txt && rm *_species_taxa.txt && rm *_species_duplicates_virome.txt
for i in $(ls *_species_duplicates.txt);do
	(
	awk -F '\t' '{print $1, $2"~"$14, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_species_duplicates*}_species_inter.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi	
done
wait
for i in $(ls *_species_inter.txt);do
	(	
	awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_species_inter*}_species_inter2.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
for i in $(ls *_species_inter2.txt);do
	(
	awk -F '\t' '{dups[$1]++} END {for (num in dups) {print num}}' $i | sort -k1,1  > ${i%_species_inter2*}_duplicate_count.txt
	) &
	if  [[ $(jobs -r -p | wc -l) -gt $N ]]; then
		wait
	fi
done
for i in $(ls *_duplicate_count.txt);do
	(
	awk -F '~' '{a[$1]++;b[$1]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i | sort -k1,1 > ${i%_duplicate_count*}_multialign_species_reads.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
		wait
	fi
done

for i in $(ls *_multialign_species_reads.txt);do
	(
	awk -F '\t'  'FNR==NR {a[$1]; next}; $1 in a' $i ${i%_multialign_species_reads*}_species_inter2.txt | sort -u -k1,1 | awk 'gsub("~","\t",$0)'| awk -F '\t' '{print $3, $1, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' > ${i%_multialign_species_reads*}_species_OTU.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
for i in $(ls *_species_unique_reads.txt);do
	awk -F '\t' '{print $11}' OFS=';' $i > ${i%_species_unique_reads*}_taxids_uniq_inter.txt
done

for i in $(ls *_uniq_inter.txt);do
	awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_uniq_inter*}_taxids_uniq.txt
done

rm *_uniq_inter.txt
cd $proj_dir/metagenome/sighits/sighits_species
for i in $(ls *_taxids_uniq.txt);do
	(
	awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' $tool_dir/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_uniq*}_uniq_inter.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi 
done
wait
for i in $(ls *_uniq_inter.txt);do
	(
	awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_uniq_inter*}_species_taxid.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
rm *_taxids_uniq.txt
for i in $(ls *_species_taxid.txt);do
	awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_species_taxid*}_species_column.txt
done
for i in $(ls *_species_column.txt);do
	paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column*}_species_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_species_taxa.txt
done
for i in $(ls *_unique_reads.txt);do
	paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%*_species_unique_reads*}_species_taxa.txt) > ${i%_uniq*}_species_unique_uncultured.txt
done
rm *_species_taxid.txt && rm *_uniq_inter.txt && rm *_species_column.txt && rm *_species_taxa.txt
for i in $(ls *_species_unique_uncultured.txt);do
	awk -F '\t' '!/Uncultured/' $i > ${i%*_species_unique_uncultured*}_species_unique_sequences.txt
done
rm *_species_unique_uncultured.txt && rm *_species_inter.txt && rm *_species_inter2.txt && rm *_duplicate_count.txt && rm *_multialign_species_reads.txt && rm *_species_duplicates.txt
for i in $(ls *_species_OTU.txt);do
	cat $i ${i%_species_OTU*}_species_species_unique_sequences.txt > ${i%_species_OTU*}_complete_species_reads.txt
done 
echo -e "${YELLOW}- compiling taxonomic information"
for i in $(ls *_complete_species_reads.txt);do
	awk -F '\t' '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' $i > ${i%_complete_species_reads*}_sighits_temp.txt
	echo $'abundance\tqseqid\tsseqid\tlength\tmismatch\tevalue\tpident\tqcovs\tqseq\tsseq\tstaxids\tstitle\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tdomain' | \
	cat - ${i%_complete_species_reads*}_sighits_temp.txt > ${i%_complete_species_reads*}_sighits_temp2.txt
done
for i in $(ls *_sighits_temp2.txt);do
	awk -F '\t' '{gsub(/ /,"_");print}' $i > ${i%_sighits_temp2*}_sighits.txt
done

rm *_complete_species_reads.txt && rm *_sighits_temp.txt && rm *_species_OTU.txt && rm *_unique_reads.txt && rm *_sighits_temp2.txt
}
if [ "$species_level" == "TRUE" ]; then
	time species_filter 2>> $proj_dir/log.out
fi

species_refseq_filter(){
cd $proj_dir/metagenome/sighits/sighits/species
mkdir species_gene_annotation
for i in $(ls *_sighits.txt);do
	cp $i $proj_dir/metagenome/sighits/sighits_species/species_gene_annotation
done
cd $proj_dir/metagenome/sighits/sighis_species/species_gene_annotation
for i in $(ls *_sighits.txt);do
	awk -F '\t' -v var="${i%*_sighits.txt}" '{print $1"\t"var"_"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20}' $i > ${i%*_sighits.txt}_refseq_sighits.txt
done
rm *_sighits.txt
for i in $(ls *_refseq_sighits.txt);do
	awk 'FNR==NR{ h=$0 }NR>1{ print (!a[$13]++? h ORS $0 : $0) >> $13".txt" }' $i
done

for i in $(ls *.txt);do
	sed -i '/^$/d' $i 
done
for i in $(ls *_refseq_sighits.txt);do
	mv $i $proj_dir/metagenome/sighits/sighits_species
done
for i in $(ls *.txt);do
	awk -F '\t' '{print ">"$2"\n"$9}' $i > ${i%.txt}.fasta
done
rm *.txt
for i in $(ls *.fasta);do
	$tool_dir/ncbi-blast-2.8.1+/bin/blastn -task megablast -query $i -db $refseq_db_dir -num_threads $threads -evalue 1e-10 -max_target_seqs 1 -outfmt "6 qseqid sseqid length mismatch evalue pident qcovs qseq sseq staxids stitle" -out ${i%*.fasta}_refseq.megablast
done
for i in $(ls *_refseq.megablast);do
	sort -u $i > ${i%*_refseq.megablast}_refseq_nr.txt
done
rm *_refseq.megablast && rm *.fasta
find . -type f -name '*_refseq_nr.txt' -exec cat {} + > refseq_hits.txt
awk -F '\t' '{print $10}' OFS=';' refseq_hits.txt > taxids_refseq_inter.txt
awk -F ';' '{print $1}' OFS='\t' taxids_refseq_inter.txt > taxids_refseq.txt
rm *_refseq_inter.txt
for i in $(ls taxids_refseq.txt);do
	(
	awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' $tool_dir/rankedlineage_edited.dmp OFS='\t' $i> ${i%taxids_refseq*}refseq_inter.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi 
done
paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' OFS='\t' refseq_hits.txt ) <(awk -F '\t' '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' refseq_inter.txt) > refseq_hits_genus.txt
rm refseq_inter.txt && rm taxids_refseq.txt && rm refseq_hits.txt
cd $proj_dir/metagenome/sighits/sighits_species/
for i in $(ls *_refseq_sighits.txt);do
	mv $i $proj_dir/metagenome/sighits/sighits_species/species_gene_annotation
done
cd $proj_dir/metagenome/sighits/sighits_species/species_gene_annotation
find . -type f -name '*_refseq_sighits.txt' -exec cat {} + > refseq_prehits_temp.txt
awk -F '\t' '!/abundance/' refseq_prehits_temp.txt > refseq_prehits.txt && rm refseq_prehits_temp.txt
awk -F '\t' '{print $2, $14}' OFS='\t' refseq_prehits.txt > prehits.txt
awk -F '\t' '{print $1, $14}' OFS='\t' refseq_hits_genus.txt > hits.txt
awk -F '\t' 'NR==FNR{c[$1$2]++;next};c[$1$2] > 0' prehits.txt hits.txt > consistent_hits_taxids.txt
awk -F '\t' '{print $1}' consistent_hits_taxids.txt > consistent_hits.txt && rm consistent_hits_taxids.txt
for i in $(ls *_refseq_sighits.txt);do
	awk -F '\t' 'NR==FNR{c[$1]++;next};c[$2] > 0' consistent_hits.txt $i > ${i%*_refseq_sighits.txt}_refseq_filtered.txt
done
rm consistent_hits.txt && rm hits.txt && rm prehits.txt && rm refseq_prehits.txt && rm *_refseq_sighits.txt && rm *refseq_hits_genus
for i in $(ls *_refseq_filtered.txt);do
	mv $i $proj_dir/metagenome/sighits/sighits_species/
done 
cd $proj_dir/metagenome/sighits/sighits_species/
for i in $(ls *refseq_filtered.txt);do
	echo $'abundance\tqseqid\tsseqid\tlength\tmismatch\tevalue\tpident\tqcovs\tqseq\tsseq\tstaxids\tstitle\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tdomain' | \
	cat - $i > ${i%*_refseq_filtered.txt}_sighits.txt
done
rm *_refseq_filtered.txt
}
if [ "$refseq_filter" == "TRUE" ] && [ "$species_level" == "TRUE" ]; then
	time species_refseq_filter 2>> $proj_dir/log.out
fi

species_stats(){
cd $proj_dir/metagenome/sighits/sighits_species
find . -type f -name '*_sighits.txt' -exec cat {} + > sighits.txt
awk -F '\t' '{print $13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20}' sighits.txt > rankedlineage_subhits.txt && rm sighits.txt
cd $proj_dir/metagenome/sighits/sighits_species/
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -u | awk '!/species/' > species_taxa_mean_temp1.txt
echo -e 'species' | cat - species_taxa_mean_temp1.txt > species_taxa_mean_temp.txt && rm species_taxa_mean_temp1.txt
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -u | awk '!/species/' > species_taxa_unique_sequences_temp1.txt
echo -e 'species' | cat - species_taxa_unique_sequences_temp1.txt > species_taxa_unique_sequences_temp.txt && rm species_taxa_unique_sequences_temp1.txt
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -u | awk '!/species/' > species_taxa_quantification_accuracy_temp1.txt
echo -e 'species' | cat - species_taxa_quantification_accuracy_temp1.txt > species_taxa_quantification_accuracy_temp.txt && rm species_taxa_quantification_accuracy_temp1.txt
species_level=species
for i in $(ls *_sighits.txt);do
	Rscript $tool_dir/Rscripts/stats_summary.R $i $min_uniq $species_level
	echo $'species\tmean\tuniq_reads\tstderr' | cat - stats1.txt > stats2.txt 
	id=${i%_sighits*}_mean && awk -v id=$id '{gsub(/mean/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$2}' > holdmean.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdmean.txt species_taxa_mean_temp.txt > holdmean2.txt && cat holdmean2.txt > species_taxa_mean_temp.txt
	id=${i%_sighits*}_uniq_reads && awk -v id=$id '{gsub(/uniq_reads/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$3}' > holduniq_reads.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holduniq_reads.txt species_taxa_unique_sequences_temp.txt > holduniq_reads2.txt && cat holduniq_reads2.txt > species_taxa_unique_sequences_temp.txt
	id=${i%_sighits*}_stderr && awk -v id=$id '{gsub(/stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$4}' > holdstderr.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdstderr.txt species_taxa_quantification_accuracy_temp.txt > holdstderr2.txt && cat holdstderr2.txt > species_taxa_quantification_accuracy_temp.txt
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $1,"\t",$2,"\t",$3,"\t",$4,"\t", a[$1]}'  rankedlineage_subhits.txt stats2.txt > stats3.txt
	awk '{print $1,$2,$3,$4,$6,$7,$8,$9,$10,$11 }' stats3.txt | awk '{gsub(/ /,"\t"); print }' > ${i%_sighits*}_taxastats.txt
	rm *stats1* *stats2* *stats3* *hold*
done
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' species_taxa_mean_temp.txt > species_taxa_mean_temp2.txt && rm species_taxa_mean_temp.txt
awk '{gsub(/ /,"\t"); print}' species_taxa_mean_temp2.txt > ../../results/species_level/species_taxa_mean.txt && rm species_taxa_mean_temp2.txt
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' species_taxa_unique_sequences_temp.txt > species_taxa_unique_sequences_temp2.txt && rm species_taxa_unique_sequences_temp.txt
awk '{gsub(/ /,"\t"); print}' species_taxa_unique_sequences_temp2.txt > ../../results/species_level/species_taxa_unique_sequences.txt && rm species_taxa_unique_sequences_temp2.txt
awk '{gsub(/ /,"\t"); print}' species_taxa_quantification_accuracy_temp.txt > species_taxa_quantification_accuracy_temp2.txt && rm species_taxa_quantification_accuracy_temp.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/species_level/species_taxa_mean.txt species_taxa_quantification_accuracy_temp2.txt > ../../results/species_level/species_taxa_quantification_accuracy.txt && rm species_taxa_quantification_accuracy_temp2.txt

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
cd $proj_dir/metagenome/results/species_level
for i in {mean,unique_sequences,quantification_accuracy}; do
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_species/rankedlineage_subhits.txt species_taxa_${i}.txt | \
	awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' > species_taxainfo_${i}.txt
done
rm *_taxa_*

echo -e "${YELLOW}- creating species-level visualizations"
species_taxa_mean=species_taxainfo_mean.txt
species_taxa_uniq=species_taxainfo_unique_sequences.txt
species_taxa_quant=species_taxainfo_quantification_accuracy.txt
percent_thresh=5
Rscript $tool_dir/Rscripts/species_level_corr.R $species_taxa_mean $percent_thresh &>/dev/null
}
if [ "$species_level" == "TRUE" ]; then
	time species_stats 2>> $proj_dir/log.out
fi
################################################################################################################
genus_filter(){
cd $proj_dir/metagenome/sighits
mkdir sighits_genus
cd $proj_dir/metagenome/results
mkdir genus_level
cp seqcov.txt ./genus_level

echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is Preforming Genus-Level Classification \n\e[97m########################################################\n"
cd $proj_dir/metagenome/alignment
echo -e "${YELLOW}- preforming exact-matching algorithm"
for i in $(ls *_haplotig_nd.megablast);do
	awk '$6>=97' $i | awk 'gsub(" ","_",$0)' > ../sighits/sighits_genus/${i%_haplotig*}_filter.txt
	awk 'NR==FNR {h[$1] = $2; next} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,h[$1]}' ../haplotig/${i%_haplotig*}_normalized.txt ../sighits/sighits_genus/${i%_haplotig*}_filter.txt | 
	awk '{print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' | awk 'gsub(" ","\t",$0)' > ../sighits/sighits_genus/${i%_haplotig*}_sighits.txt
done
rm ../sighits/sighits_genus/*_filter.txt 
echo -e "${YELLOW}- compiling genus-level OTUs with multi-alignment algorithm"
cd $proj_dir/metagenome/sighits/sighits_genus/
for i in $(ls *_sighits.txt);do
	awk -F '\t' '{a[$2]++;b[$2]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i > ${i%_sighits*}_genus_unique_reads.txt
done
cd $proj_dir/metagenome/sighits/sighits_genus
for i in $(ls *_sighits.txt);do
	awk -F'|' 'FNR==NR{a[$1,$2]=1; next}  !a[$1,$2]' ${i%_sighits.txt}_genus_unique_reads.txt $i > ${i%_sighits*}_dup.txt
done
cd $proj_dir/metagenome/sighits/sighits_genus
for i in $(ls *_dup.txt);do
	awk -F '\t' '{print $11}' OFS=';' $i > ${i%_dup*}_taxids_dup_inter.txt
done
for i in $(ls *_dup_inter.txt);do
	awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_dup_inter*}_taxids_dup.txt
done
rm *_taxids_dup_inter.txt
cd $proj_dir/metagenome/sighits/sighits_genus
for i in $(ls *_taxids_dup.txt);do
	(
	awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' $tool_dir/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_dup*}_dup_inter.txt 
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
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
	paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' OFS='\t' ${i%*_dup.txt}_species_taxa.txt) > ${i%_dup*}_genus_duplicates_uncultured.txt
done

rm *_species_taxid.txt && rm *_dup_inter.txt && rm *_dup.txt && rm *_species_column.txt && rm *_species_taxa.txt
for i in $(ls *_genus_duplicates_uncultured.txt);do
	awk -F '\t' '!/Uncultured/' $i > ${i%*_genus_duplicates_uncultured*}_genus_duplicates.txt
done
rm *_genus_duplicates_uncultured.txt
for i in $(ls *_genus_duplicates.txt);do
	(
	awk -F '\t' '{print $1, $2"~"$15, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_genus_duplicates*}_genus_inter.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
for i in $(ls *_genus_inter.txt);do
	(
	awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_genus_inter*}_genus_inter2.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
for i in $(ls *_genus_inter2.txt);do
	(
	awk -F '\t' '{dups[$1]++} END {for (num in dups) {print num}}' $i | sort -k1,1  > ${i%_genus_inter2*}_duplicate_count.txt 
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
for i in $(ls *_duplicate_count.txt);do
	(
	awk -F '~' '{a[$1]++;b[$1]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i | sort -k1,1 > ${i%_duplicate_count*}_multialign_genus_reads.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
for i in $(ls *_multialign_genus_reads.txt);do
	(
	awk -F '\t'  'FNR==NR {a[$1]; next}; $1 in a' $i ${i%_multialign_genus_reads*}_genus_inter2.txt | sort -u -k1,1 | awk 'gsub("~","\t",$0)'| awk -F '\t' '{print $3, $1, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' > ${i%_multialign_genus_reads*}_genus_OTU.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
for i in $(ls *_genus_unique_reads.txt);do
	awk -F '\t' '{print $11}' OFS=';' $i > ${i%_genus_unique_reads*}_taxids_uniq_inter.txt
done

for i in $(ls *_uniq_inter.txt);do
	awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_uniq_inter*}_taxids_uniq.txt
done

rm *_taxids_uniq_inter.txt
cd $proj_dir/metagenome/sighits/sighits_genus
for i in $(ls *_taxids_uniq.txt);do
	(
	awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' $tool_dir/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_uniq*}_uniq_inter.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi 
done
wait
for i in $(ls *_uniq_inter.txt);do
	(
	awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_uniq_inter*}_species_taxid.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
rm *_taxids_uniq.txt
for i in $(ls *_species_taxid.txt);do
	awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_species_taxid*}_species_column.txt
done
for i in $(ls *_species_column.txt);do
	paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column*}_species_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_species_taxa.txt
done
for i in $(ls *_unique_reads.txt);do
	paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%*_genus_unique_reads*}_species_taxa.txt) > ${i%_uniq*}_genus_unique_uncultured.txt
done
rm *_species_taxid.txt && rm *_uniq_inter.txt && rm *_species_column.txt && rm *_species_taxa.txt
for i in $(ls *_genus_unique_uncultured.txt);do
	awk -F '\t' '!/Uncultured/' $i > ${i%*_genus_unique_uncultured*}_genus_unique_sequences.txt
done
rm *_genus_unique_uncultured.txt && rm *_genus_inter.txt && rm *_genus_inter2.txt && rm *_duplicate_count.txt && rm *_multialign_genus_reads.txt && rm *_genus_duplicates.txt
for i in $(ls *_genus_OTU.txt);do
	cat $i ${i%_genus_OTU*}_genus_genus_unique_sequences.txt > ${i%_genus_OTU*}_complete_genus_reads.txt
done 
echo -e "${YELLOW}- compiling taxonomic information"
for i in $(ls *_complete_genus_reads.txt);do
	awk -F '\t' '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' $i > ${i%_complete_genus_reads*}_sighits_temp.txt
	echo $'abundance\tqseqid\tsseqid\tlength\tmismatch\tevalue\tpident\tqcovs\tqseq\tsseq\tstaxids\tstitle\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tdomain' | \
	cat - ${i%_complete_genus_reads*}_sighits_temp.txt > ${i%_complete_genus_reads*}_sighits.txt
done
rm *_complete_genus_reads.txt && rm *_sighits_temp.txt && rm *_genus_OTU.txt && rm *_unique_reads.txt && rm *_genus_genus_unique_sequences.txt
cd $proj_dir/metagenome/sighits/sighits_genus
find . -type f -name '*_sighits.txt' -exec cat {} + > sighits.txt
awk -F '\t' '{print $14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20}' sighits.txt > rankedlineage_subhits_temp.txt && rm sighits.txt
awk -F '\t' '!/genus/' rankedlineage_subhits_temp.txt > rankedlineage_subhits_temp2.txt && rankedlineage_subhits_temp.txt
echo $'genus\tfamily\torder\tclass\tphylum\tkingdom\tdomain' | cat - rankedlineage_subhits_temp2.txt > rankedlineage_subhits.txt && rankedlineage_subhits_temp2.txt
}
if [ "$genus_level" == "TRUE" ]; then
	time genus_filter 2>> $proj_dir/log.out
fi
genus_refseq_filter(){
cd $proj_dir/metagenome/sighits/sighits/genus
mkdir genus_gene_annotation
for i in $(ls *_sighits.txt);do
	cp $i $proj_dir/metagenome/sighits/sighits_genus/genus_gene_annotation
done
cd $proj_dir/metagenome/sighits/sighis_genus/genus_gene_annotation
for i in $(ls *_sighits.txt);do
	awk -F '\t' -v var="${i%*_sighits.txt}" '{print $1"\t"var"_"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20}' $i > ${i%*_sighits.txt}_refseq_sighits.txt
done
rm *_sighits.txt
for i in $(ls *_refseq_sighits.txt);do
	awk 'FNR==NR{ h=$0 }NR>1{ print (!a[$14]++? h ORS $0 : $0) >> $14".txt" }' $i
done

for i in $(ls *.txt);do
	sed -i '/^$/d' $i 
done
for i in $(ls *_refseq_sighits.txt);do
	mv $i $proj_dir/metagenome/sighits/sighits_genus
done
for i in $(ls *.txt);do
	awk -F '\t' '{print ">"$2"\n"$9}' $i > ${i%.txt}.fasta
done
rm *.txt
for i in $(ls *.fasta);do
	$tool_dir/ncbi-blast-2.8.1+/bin/blastn -task megablast -query $i -db $refseq_db_dir -num_threads $threads -evalue 1e-10 -max_target_seqs 1 -outfmt "6 qseqid sseqid length mismatch evalue pident qcovs qseq sseq staxids stitle" -out ${i%*.fasta}_refseq.megablast
done
for i in $(ls *_refseq.megablast);do
	sort -u $i > ${i%*_refseq.megablast}_refseq_nr.txt
done
rm *_refseq.megablast && rm *.fasta
find . -type f -name '*_refseq_nr.txt' -exec cat {} + > refseq_hits.txt
awk -F '\t' '{print $10}' OFS=';' refseq_hits.txt > taxids_refseq_inter.txt
awk -F ';' '{print $1}' OFS='\t' taxids_refseq_inter.txt > taxids_refseq.txt
rm *_refseq_inter.txt
for i in $(ls taxids_refseq.txt);do
	(
	awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' $tool_dir/rankedlineage_edited.dmp OFS='\t' $i> ${i%taxids_refseq*}refseq_inter.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi 
done
paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' OFS='\t' refseq_hits.txt ) <(awk -F '\t' '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' refseq_inter.txt) > refseq_hits_genus.txt
rm refseq_inter.txt && rm taxids_refseq.txt && rm refseq_hits.txt
cd $proj_dir/metagenome/sighits/sighits_genus/
for i in $(ls *_refseq_sighits.txt);do
	mv $i $proj_dir/metagenome/sighits/sighits_genus/genus_gene_annotation
done
cd $proj_dir/metagenome/sighits/sighits_genus/genus_gene_annotation
find . -type f -name '*_refseq_sighits.txt' -exec cat {} + > refseq_prehits_temp.txt
awk -F '\t' '!/abundance/' refseq_prehits_temp.txt > refseq_prehits.txt && rm refseq_prehits_temp.txt
awk -F '\t' '{print $2, $14}' OFS='\t' refseq_prehits.txt > prehits.txt
awk -F '\t' '{print $1, $14}' OFS='\t' refseq_hits_genus.txt > hits.txt
awk -F '\t' 'NR==FNR{c[$1$2]++;next};c[$1$2] > 0' prehits.txt hits.txt > consistent_hits_taxids.txt
awk -F '\t' '{print $1}' consistent_hits_taxids.txt > consistent_hits.txt && rm consistent_hits_taxids.txt
for i in $(ls *_refseq_sighits.txt);do
	awk -F '\t' 'NR==FNR{c[$1]++;next};c[$2] > 0' consistent_hits.txt $i > ${i%*_refseq_sighits.txt}_refseq_filtered.txt
done
rm consistent_hits.txt && rm hits.txt && rm prehits.txt && rm refseq_prehits.txt && rm *_refseq_sighits.txt && rm *refseq_hits_genus
for i in $(ls *_refseq_filtered.txt);do
	mv $i $proj_dir/metagenome/sighits/sighits_genus/
done 
cd $proj_dir/metagenome/sighits/sighits_genus/
for i in $(ls *refseq_filtered.txt);do
	echo $'abundance\tqseqid\tsseqid\tlength\tmismatch\tevalue\tpident\tqcovs\tqseq\tsseq\tstaxids\tstitle\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tdomain' | \
	cat - $i > ${i%*_refseq_filtered.txt}_sighits.txt
done
rm *_refseq_filtered.txt
}
if [ "$refseq_filter" == "TRUE" ] && [ "$genus_level" == "TRUE" ]; then
	time genus_refseq_filter 2>> $proj_dir/log.out
fi

genus_stats(){
cd $proj_dir/metagenome/sighits/sighits_genus
echo -e "${YELLOW}- quantifying the genus-level taxonomy"
awk '{print $1}' rankedlineage_subhits.txt | sort -u > genus_taxa_mean_temp1.txt
echo -e 'genus' | cat - genus_taxa_mean_temp1.txt > genus_taxa_mean_temp.txt && rm genus_taxa_mean_temp1.txt
awk '{print $1}' rankedlineage_subhits.txt | sort -u > genus_taxa_unique_sequences_temp1.txt
echo -e 'genus' | cat - genus_taxa_unique_sequences_temp1.txt > genus_taxa_unique_sequences_temp.txt && rm genus_taxa_unique_sequences_temp1.txt
awk '{print $1}' rankedlineage_subhits.txt | sort -u > genus_taxa_quantification_accuracy_temp1.txt
echo -e 'genus' | cat - genus_taxa_quantification_accuracy_temp1.txt > genus_taxa_quantification_accuracy_temp.txt && rm genus_taxa_quantification_accuracy_temp1.txt 
genus_level=genus
for i in $(ls *_sighits.txt);do
	Rscript $tool_dir/Rscripts/stats_summary.R $i $min_uniq $genus_level
	echo $'genus\tmean\tuniq_reads\tstderr' | cat - stats1.txt > stats2.txt 
	id=${i%_sighits*}_mean && awk -v id=$id '{gsub(/mean/,id); print }' stats2.txt | awk '{print $1,"\t",$2}' > holdmean.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdmean.txt genus_taxa_mean_temp.txt > holdmean2.txt && cat holdmean2.txt > genus_taxa_mean_temp.txt
	id=${i%_sighits*}_uniq_reads && awk -v id=$id '{gsub(/uniq_reads/,id); print }' stats2.txt | awk '{print $1,"\t",$3}' > holduniq_reads.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holduniq_reads.txt genus_taxa_unique_sequences_temp.txt > holduniq_reads2.txt && cat holduniq_reads2.txt > genus_taxa_unique_sequences_temp.txt
	id=${i%_sighits*}_stderr && awk -v id=$id '{gsub(/stderr/,id); print }' stats2.txt | awk '{print $1,"\t",$4}' > holdstderr.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdstderr.txt genus_taxa_quantification_accuracy_temp.txt > holdstderr2.txt && cat holdstderr2.txt > genus_taxa_quantification_accuracy_temp.txt
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $1,"\t",$2,"\t",$3,"\t",$4,"\t", a[$1]}'  rankedlineage_subhits.txt stats2.txt > stats3.txt
	awk '{print $1,$2,$3,$4,$6,$7,$8,$9,$10 }' stats3.txt | awk '{gsub(/ /,"\t"); print }' > ${i%_sighits*}_taxastats.txt
	rm *stats1* *stats2* *stats3* *hold*
done
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' genus_taxa_mean_temp.txt > genus_taxa_mean_temp2.txt && rm genus_taxa_mean_temp.txt
awk '{gsub(/ /,"\t"); print}' genus_taxa_mean_temp2.txt > ../../results/genus_level/genus_taxa_mean.txt && rm genus_taxa_mean_temp2.txt
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' genus_taxa_unique_sequences_temp.txt > genus_taxa_unique_sequences_temp2.txt && rm genus_taxa_unique_sequences_temp.txt
awk '{gsub(/ /,"\t"); print}' genus_taxa_unique_sequences_temp2.txt > ../../results/genus_level/genus_taxa_unique_sequences.txt && rm genus_taxa_unique_sequences_temp2.txt
awk '{gsub(/ /,"\t"); print}' genus_taxa_quantification_accuracy_temp.txt > genus_taxa_quantification_accuracy_temp2.txt && rm genus_taxa_quantification_accuracy_temp.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/genus_level/genus_taxa_mean.txt genus_taxa_quantification_accuracy_temp2.txt > genus_taxa_quantification_accuracy_temp3.txt && rm genus_taxa_quantification_accuracy_temp2.txt
sed '2,${/genus/d;}' genus_taxa_quantification_accuracy_temp3.txt > ../../results/genus_level/genus_taxa_quantification_accuracy.txt && rm genus_taxa_quantification_accuracy_temp3.txt

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
cd $proj_dir/metagenome/results/genus_level
for i in {mean,unique_sequences,quantification_accuracy}; do
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_genus/rankedlineage_subhits.txt genus_taxa_${i}.txt | \
	awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' > genus_taxainfo_${i}.txt
done
rm *_taxa_*
echo -e "${YELLOW}- creating genus-level visualizations"
genus_taxa_mean=genus_taxainfo_mean.txt
genus_taxa_uniq=genus_taxainfo_unique_sequences.txt
genus_taxa_quant=genus_taxainfo_quantification_accuracy.txt
percent_thresh=5
genus_mean=genus_mean.txt
genus_uniq=genus_unique_sequences.txt
genus_quant=genus_quantification_accuracy.txt
Rscript $tool_dir/Rscripts/genus_level_corr.R $genus_mean $percent_thresh &>/dev/null
}
if [ "$genus_level" == "TRUE" ]; then
	time genus_stats 2>> $proj_dir/log.out
fi
################################################################################################################
family_filter(){
cd $proj_dir/metagenome/sighits
mkdir sighits_family
cd $proj_dir/metagenome/results
mkdir family_level
cp seqcov.txt ./family_level
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is Preforming Family-Level Classification \n\e[97m########################################################\n"
cd $proj_dir/metagenome/alignment
echo -e "${YELLOW}- preforming exact-matching algorithm"
for i in $(ls *_haplotig_nd.megablast);do
	awk '$6>=97' $i | awk 'gsub(" ","_",$0)' > ../sighits/sighits_family/${i%_haplotig*}_filter.txt
	awk 'NR==FNR {h[$1] = $2; next} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,h[$1]}' ../haplotig/${i%_haplotig*}_normalized.txt ../sighits/sighits_family/${i%_haplotig*}_filter.txt | 
	awk '{print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' | awk 'gsub(" ","\t",$0)' > ../sighits/sighits_family/${i%_haplotig*}_sighits.txt
done
rm ../sighits/sighits_family/*_filter.txt 
echo -e "${YELLOW}- compiling family-level OTUs with multi-alignment algorithm"
cd $proj_dir/metagenome/sighits/sighits_family/
for i in $(ls *_sighits.txt);do
	awk -F '\t' '{a[$2]++;b[$2]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i > ${i%_sighits*}_family_unique_reads.txt
done
cd $proj_dir/metagenome/sighits/sighits_family
for i in $(ls *_sighits.txt);do
	awk -F'|' 'FNR==NR{a[$1,$2]=1; next}  !a[$1,$2]' ${i%_sighits.txt}_family_unique_reads.txt $i > ${i%_sighits*}_dup.txt
done
cd $proj_dir/metagenome/sighits/sighits_family
for i in $(ls *_dup.txt);do
	awk -F '\t' '{print $11}' OFS=';' $i > ${i%_dup*}_taxids_dup_inter.txt
done

for i in $(ls *_dup_inter.txt);do
	awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_dup_inter*}_taxids_dup.txt
done

rm *_taxids_dup_inter.txt
cd $proj_dir/metagenome/sighits/sighits_family
for i in $(ls *_taxids_dup.txt);do
	(
	awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' $tool_dir/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_dup*}_dup_inter.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi 
done
wait
for i in $(ls *_dup_inter.txt);do
	(
	awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_dup_inter*}_species_taxid.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
rm *_taxids_dup.txt
for i in $(ls *_species_taxid.txt);do
	awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_species_taxid*}_species_column.txt
done
for i in $(ls *_species_column.txt);do
	paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column*}_species_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_species_taxa.txt
done
for i in $(ls *_dup.txt);do
	paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' OFS='\t' ${i%*_dup.txt}_species_taxa.txt) > ${i%_dup*}_family_duplicates_uncultured.txt
done
rm *_species_taxid.txt && rm *_dup_inter.txt && rm *_dup.txt && rm *_species_column.txt && rm *_species_taxa.txt
for i in $(ls *_family_duplicates_uncultured.txt);do
	awk -F '\t' '!/Uncultured/' $i > ${i%*_family_duplicates_uncultured*}_family_duplicates.txt
done
rm *_family_duplicates_uncultured.txt
for i in $(ls *_family_duplicates.txt);do
	(
	awk -F '\t' '{print $1, $2"~"$16, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_family_duplicates*}_family_inter.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
for i in $(ls *_family_inter.txt);do
	(
	awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_family_inter*}_family_inter2.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
for i in $(ls *_family_inter2.txt);do
	(
	awk -F '\t' '{dups[$1]++} END {for (num in dups) {print num}}' $i | sort -k1,1  > ${i%_family_inter2*}_duplicate_count.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
for i in $(ls *_duplicate_count.txt);do
	(
	awk -F '~' '{a[$1]++;b[$1]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i | sort -k1,1 > ${i%_duplicate_count*}_multialign_family_reads.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
for i in $(ls *_multialign_family_reads.txt);do
	(
	awk -F '\t'  'FNR==NR {a[$1]; next}; $1 in a' $i ${i%_multialign_family_reads*}_family_inter2.txt | sort -u -k1,1 | awk 'gsub("~","\t",$0)'| awk -F '\t' '{print $3, $1, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' > ${i%_multialign_family_reads*}_family_OTU.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait

for i in $(ls *_family_unique_reads.txt);do
	awk -F '\t' '{print $11}' OFS=';' $i > ${i%_family_unique_reads*}_taxids_uniq_inter.txt
done

for i in $(ls *_uniq_inter.txt);do
	awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_uniq_inter*}_taxids_uniq.txt
done

rm *_taxids_uniq_inter.txt
cd $proj_dir/metagenome/sighits/sighits_family
for i in $(ls *_taxids_uniq.txt);do
	(
	awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' $tool_dir/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_uniq*}_uniq_inter.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi 
done
wait
for i in $(ls *_uniq_inter.txt);do
	(
	awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_uniq_inter*}_species_taxid.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
rm *_taxids_uniq.txt
for i in $(ls *_species_taxid.txt);do
	awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_species_taxid*}_species_column.txt
done
for i in $(ls *_species_column.txt);do
	paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column*}_species_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_species_taxa.txt
done
for i in $(ls *_unique_reads.txt);do
	paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%*_family_unique_reads*}_species_taxa.txt) > ${i%_uniq*}_family_unique_uncultured.txt
done
rm *_species_taxid.txt && rm *_uniq_inter.txt && rm *_uniq.txt && rm *_species_column.txt && rm *_species_taxa.txt
for i in $(ls *_family_unique_uncultured.txt);do
	awk -F '\t' '!/Uncultured/' $i > ${i%*_family_unique_uncultured*}_family_unique_sequences.txt
done
rm *_family_unique_uncultured.txt && rm *_family_inter.txt && rm *_family_inter2.txt && rm *_duplicate_count.txt && rm *_multialign_family_reads.txt && rm *_family_duplicates.txt
for i in $(ls *_family_OTU.txt);do
	(
	cat $i ${i%_family_OTU*}_family_family_unique_sequences.txt > ${i%_family_OTU*}_complete_family_reads.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
rm *_family_family_unique_sequences.txt && rm *_family_OTU.txt && rm *_family_unique_reads.txt
for i in $(ls *_complete_family_reads.txt);do
	awk -F '\t' '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' $i > ${i%_complete_family_reads*}_sighits_temp.txt
	echo $'abundance\tqseqid\tsseqid\tlength\tmismatch\tevalue\tpident\tqcovs\tqseq\tsseq\tstaxids\tstitle\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tdomain' | \
	cat - ${i%_complete_family_reads*}_sighits_temp.txt > ${i%_complete_family_reads*}_sighits.txt
done
}
if [ "$family_level" == "TRUE" ]; then
	time family_filter 2>> $proj_dir/log.out
fi
family_refseq_filter(){
cd $proj_dir/metagenome/sighits/sighits/family
mkdir family_gene_annotation
for i in $(ls *_sighits.txt);do
	cp $i $proj_dir/metagenome/sighits/sighits_family/family_gene_annotation
done
cd $proj_dir/metagenome/sighits/sighis_family/family_gene_annotation
for i in $(ls *_sighits.txt);do
	awk -F '\t' -v var="${i%*_sighits.txt}" '{print $1"\t"var"_"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20}' $i > ${i%*_sighits.txt}_refseq_sighits.txt
done
rm *_sighits.txt
for i in $(ls *_refseq_sighits.txt);do
	awk 'FNR==NR{ h=$0 }NR>1{ print (!a[$15]++? h ORS $0 : $0) >> $15".txt" }' $i
done

for i in $(ls *.txt);do
	sed -i '/^$/d' $i 
done
for i in $(ls *_refseq_sighits.txt);do
	mv $i $proj_dir/metagenome/sighits/sighits_family
done
for i in $(ls *.txt);do
	awk -F '\t' '{print ">"$2"\n"$9}' $i > ${i%.txt}.fasta
done
rm *.txt
for i in $(ls *.fasta);do
	$tool_dir/ncbi-blast-2.8.1+/bin/blastn -task megablast -query $i -db $refseq_db_dir -num_threads $threads -evalue 1e-10 -max_target_seqs 1 -outfmt "6 qseqid sseqid length mismatch evalue pident qcovs qseq sseq staxids stitle" -out ${i%*.fasta}_refseq.megablast
done
for i in $(ls *_refseq.megablast);do
	sort -u $i > ${i%*_refseq.megablast}_refseq_nr.txt
done
rm *_refseq.megablast && rm *.fasta
find . -type f -name '*_refseq_nr.txt' -exec cat {} + > refseq_hits.txt
awk -F '\t' '{print $10}' OFS=';' refseq_hits.txt > taxids_refseq_inter.txt
awk -F ';' '{print $1}' OFS='\t' taxids_refseq_inter.txt > taxids_refseq.txt
rm *_refseq_inter.txt
for i in $(ls taxids_refseq.txt);do
	(
	awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' $tool_dir/rankedlineage_edited.dmp OFS='\t' $i> ${i%taxids_refseq*}refseq_inter.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi 
done
paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' OFS='\t' refseq_hits.txt ) <(awk -F '\t' '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' refseq_inter.txt) > refseq_hits_genus.txt
rm refseq_inter.txt && rm taxids_refseq.txt && rm refseq_hits.txt
cd $proj_dir/metagenome/sighits/sighits_family/
for i in $(ls *_refseq_sighits.txt);do
	mv $i $proj_dir/metagenome/sighits/sighits_family/family_gene_annotation
done
cd $proj_dir/metagenome/sighits/sighits_family/family_gene_annotation
find . -type f -name '*_refseq_sighits.txt' -exec cat {} + > refseq_prehits_temp.txt
awk -F '\t' '!/abundance/' refseq_prehits_temp.txt > refseq_prehits.txt && rm refseq_prehits_temp.txt
awk -F '\t' '{print $2, $14}' OFS='\t' refseq_prehits.txt > prehits.txt
awk -F '\t' '{print $1, $14}' OFS='\t' refseq_hits_genus.txt > hits.txt
awk -F '\t' 'NR==FNR{c[$1$2]++;next};c[$1$2] > 0' prehits.txt hits.txt > consistent_hits_taxids.txt
awk -F '\t' '{print $1}' consistent_hits_taxids.txt > consistent_hits.txt && rm consistent_hits_taxids.txt
for i in $(ls *_refseq_sighits.txt);do
	awk -F '\t' 'NR==FNR{c[$1]++;next};c[$2] > 0' consistent_hits.txt $i > ${i%*_refseq_sighits.txt}_refseq_filtered.txt
done
rm consistent_hits.txt && rm hits.txt && rm prehits.txt && rm refseq_prehits.txt && rm *_refseq_sighits.txt && rm *refseq_hits_family
for i in $(ls *_refseq_filtered.txt);do
	mv $i $proj_dir/metagenome/sighits/sighits_family/
done 
cd $proj_dir/metagenome/sighits/sighits_family/
for i in $(ls *refseq_filtered.txt);do
	echo $'abundance\tqseqid\tsseqid\tlength\tmismatch\tevalue\tpident\tqcovs\tqseq\tsseq\tstaxids\tstitle\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tdomain' | \
	cat - $i > ${i%*_refseq_filtered.txt}_sighits.txt
done
rm *_refseq_filtered.txt
}
if [ "$refseq_filter" == "TRUE" ] && [ "$family_level" == "TRUE" ]; then
	time family_refseq_filter 2>> $proj_dir/log.out
fi

family_stats(){
cd $proj_dir/metagenome/sighits/sighits_family
echo -e "${YELLOW}- compiling taxonomic information"
find . -type f -name '*_sighits.txt' -exec cat {} + > sighits.txt
awk -F '\t' '{print $15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20}' sighits.txt > rankedlineage_subhits_temp.txt && rm sighits.txt
awk -F '\t' '!/family/' rankedlineage_subhits_temp.txt > rankedlineage_subhits_temp2.txt && rankedlineage_subhits_temp.txt
echo $'family\torder\tclass\tphylum\tkingdom\tdomain' | cat - rankedlineage_subhits_temp2.txt > rankedlineage_subhits.txt && rankedlineage_subhits_temp2.txt


cd $proj_dir/metagenome/sighits/sighits_family
echo -e "${YELLOW}- quantifying the family-level taxonomy"
awk '{print $1}' rankedlineage_subhits.txt | sort -u > family_taxa_mean_temp1.txt
echo -e 'family' | cat - family_taxa_mean_temp1.txt > family_taxa_mean_temp.txt && rm family_taxa_mean_temp1.txt
awk '{print $1}' rankedlineage_subhits.txt | sort -u > family_taxa_unique_sequences_temp1.txt
echo -e 'family' | cat - family_taxa_unique_sequences_temp1.txt > family_taxa_unique_sequences_temp.txt && rm family_taxa_unique_sequences_temp1.txt
awk '{print $1}' rankedlineage_subhits.txt | sort -u > family_taxa_quantification_accuracy_temp1.txt
echo -e 'family' | cat - family_taxa_quantification_accuracy_temp1.txt > family_taxa_quantification_accuracy_temp.txt && rm family_taxa_quantification_accuracy_temp1.txt 
family_level=family
for i in $(ls *_sighits.txt);do
	Rscript $tool_dir/Rscripts/stats_summary.R $i $min_uniq $family_level
	echo $'family\tmean\tuniq_reads\tstderr' | cat - stats1.txt > stats2.txt 
	id=${i%_sighits*}_mean && awk -v id=$id '{gsub(/mean/,id); print }' stats2.txt | awk '{print $1,"\t",$2}' > holdmean.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdmean.txt family_taxa_mean_temp.txt > holdmean2.txt && cat holdmean2.txt > family_taxa_mean_temp.txt
	id=${i%_sighits*}_uniq_reads && awk -v id=$id '{gsub(/uniq_reads/,id); print }' stats2.txt | awk '{print $1,"\t",$3}' > holduniq_reads.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holduniq_reads.txt family_taxa_unique_sequences_temp.txt > holduniq_reads2.txt && cat holduniq_reads2.txt > family_taxa_unique_sequences_temp.txt
	id=${i%_sighits*}_stderr && awk -v id=$id '{gsub(/stderr/,id); print }' stats2.txt | awk '{print $1,"\t",$4}' > holdstderr.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdstderr.txt family_taxa_quantification_accuracy_temp.txt > holdstderr2.txt && cat holdstderr2.txt > family_taxa_quantification_accuracy_temp.txt
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $1,"\t",$2,"\t",$3,"\t",$4,"\t", a[$1]}'  rankedlineage_subhits.txt stats2.txt > stats3.txt
	awk '{print $1,$2,$3,$4,$6,$7,$8,$9 }' stats3.txt | awk '{gsub(/ /,"\t"); print }' > ${i%_sighits*}_taxastats.txt
	rm *stats1* *stats2* *stats3* *hold*
done
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' family_taxa_mean_temp.txt > family_taxa_mean_temp2.txt && rm family_taxa_mean_temp.txt
awk '{gsub(/ /,"\t"); print}' family_taxa_mean_temp2.txt > ../../results/family_level/family_taxa_mean.txt && rm family_taxa_mean_temp2.txt
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' family_taxa_unique_sequences_temp.txt > family_taxa_unique_sequences_temp2.txt && rm family_taxa_unique_sequences_temp.txt
awk '{gsub(/ /,"\t"); print}' family_taxa_unique_sequences_temp2.txt > ../../results/family_level/family_taxa_unique_sequences.txt && rm family_taxa_unique_sequences_temp2.txt
awk '{gsub(/ /,"\t"); print}' family_taxa_quantification_accuracy_temp.txt > family_taxa_quantification_accuracy_temp2.txt && rm family_taxa_quantification_accuracy_temp.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/family_level/family_taxa_mean.txt family_taxa_quantification_accuracy_temp2.txt > family_taxa_quantification_accuracy_temp3.txt && rm family_taxa_quantification_accuracy_temp2.txt
sed '2,${/family/d;}' family_taxa_quantification_accuracy_temp3.txt > ../../results/family_level/family_taxa_quantification_accuracy.txt && rm family_taxa_quantification_accuracy_temp3.txt

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
cd $proj_dir/metagenome/results/family_level
for i in {mean,unique_sequences,quantification_accuracy}; do
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_family/rankedlineage_subhits.txt family_taxa_${i}.txt | \
	awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' > family_taxainfo_${i}.txt
done
rm *_taxa_*
echo -e "${YELLOW}- creating family-level visualizations"
family_taxa_mean=family_taxainfo_mean.txt
family_taxa_uniq=family_taxainfo_unique_sequences.txt
family_taxa_quant=family_taxainfo_quantification_accuracy.txt
family_mean=family_mean.txt
family_uniq=family_unique_sequences.txt
family_quant=family_quantification_accuracy.txt
percent_thresh=5
Rscript $tool_dir/Rscripts/family_level_corr.R $family_mean $percent_thresh &>/dev/null
}
if [ "$family_level" == "TRUE" ]; then
	time family_stats 2>> $proj_dir/log.out
fi

order_filter(){
cd $proj_dir/metagenome/sighits
mkdir sighits_order
cd $proj_dir/metagenome/results
mkdir order_level
cp seqcov.txt ./order_level
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is Preforming Order-Level Classification \n\e[97m########################################################\n"
cd $proj_dir/metagenome/alignment
echo -e "${YELLOW}- preforming exact-matching algorithm"
for i in $(ls *_haplotig_nd.megablast);do
	awk '$6>=97' $i | awk 'gsub(" ","_",$0)' > ../sighits/sighits_order/${i%_haplotig*}_filter.txt
	awk 'NR==FNR {h[$1] = $2; next} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,h[$1]}' ../haplotig/${i%_haplotig*}_normalized.txt ../sighits/sighits_order/${i%_haplotig*}_filter.txt | 
	awk '{print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' | awk 'gsub(" ","\t",$0)' > ../sighits/sighits_order/${i%_haplotig*}_sighits.txt
done
rm ../sighits/sighits_order/*_filter.txt 
echo -e "${YELLOW}- compiling order-level OTUs with multi-alignment algorithm"
cd $proj_dir/metagenome/sighits/sighits_order/
for i in $(ls *_sighits.txt);do
	awk -F '\t' '{a[$2]++;b[$2]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i > ${i%_sighits*}_order_unique_reads.txt
done
cd $proj_dir/metagenome/sighits/sighits_order
for i in $(ls *_sighits.txt);do
		awk -F'|' 'FNR==NR{a[$1,$2]=1; next}  !a[$1,$2]' ${i%_sighits.txt}_order_unique_reads.txt $i > ${i%_sighits*}_dup.txt
done
cd $proj_dir/metagenome/sighits/sighits_order
for i in $(ls *_dup.txt);do
	awk -F '\t' '{print $11}' OFS=';' $i > ${i%_dup*}_taxids_dup_inter.txt
done

for i in $(ls *_dup_inter.txt);do
	awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_dup_inter*}_taxids_dup.txt
done

rm *_taxids_dup_inter.txt
cd $proj_dir/metagenome/sighits/sighits_order
for i in $(ls *_taxids_dup.txt);do
	(
	awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' /home/brandon/Desktop/Qmatey/tools/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_dup*}_dup_inter.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi 
done
wait
for i in $(ls *_dup_inter.txt);do
	(
	awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_dup_inter*}_species_taxid.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
		wait
	fi
done
rm *_taxids_dup.txt
for i in $(ls *_species_taxid.txt);do
	awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_species_taxid*}_species_column.txt
done
for i in $(ls *_species_column.txt);do
	paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column*}_species_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_species_taxa.txt
done
for i in $(ls *_dup.txt);do
	paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' OFS='\t' ${i%*_dup.txt}_species_taxa.txt) > ${i%_dup*}_order_duplicates_uncultured.txt
done
rm *_species_taxid.txt && rm *_dup_inter.txt && rm *_dup.txt && rm *_species_column.txt && rm *_species_taxa.txt
for i in $(ls *_order_duplicates_uncultured.txt);do
	awk -F '\t' '!/Uncultured/' $i > ${i%*_order_duplicates_uncultured*}_order_duplicates.txt
done
rm *_order_duplicates_uncultured.txt
for i in $(ls *_order_duplicates.txt);do
	(
	awk -F '\t' '{print $1, $2"~"$17, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_order_duplicates*}_order_inter.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
for i in $(ls *_order_inter.txt);do
	(
	awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_order_inter*}_order_inter2.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
for i in $(ls *_order_inter2.txt);do
	(
	awk -F '\t' '{dups[$1]++} END {for (num in dups) {print num}}' $i | sort -k1,1  > ${i%_order_inter2*}_duplicate_count.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
for i in $(ls *_duplicate_count.txt);do
	(
	awk -F '~' '{a[$1]++;b[$1]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i | sort -k1,1 > ${i%_duplicate_count*}_multialign_order_reads.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
for i in $(ls *_multialign_order_reads.txt);do
	(
	awk -F '\t'  'FNR==NR {a[$1]; next}; $1 in a' $i ${i%_multialign_order_reads*}_order_inter2.txt | sort -u -k1,1 | awk 'gsub("~","\t",$0)'| awk -F '\t' '{print $3, $1, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' > ${i%_multialign_order_reads*}_order_OTU.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
rm *_order_inter.txt && rm *_order_inter2.txt && rm *_duplicate_count.txt && rm *_multialign_order_reads.txt && rm *_order_duplicates.txt
for i in $(ls *_order_unique_reads.txt);do
	awk -F '\t' '{print $11}' OFS=';' $i > ${i%_order_unique_reads*}_taxids_uniq_inter.txt
done

for i in $(ls *_uniq_inter.txt);do
	awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_uniq_inter*}_taxids_uniq.txt
done

rm *_taxids_uniq_inter.txt
cd $proj_dir/metagenome/sighits/sighits_order
for i in $(ls *_taxids_uniq.txt);do
	(
	awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' $tool_dir/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_uniq*}_uniq_inter.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi 
done
wait
for i in $(ls *_uniq_inter.txt);do
	(
	awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_uniq_inter*}_species_taxid.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
rm *_taxids_uniq.txt
for i in $(ls *_species_taxid.txt);do
	awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_species_taxid*}_species_column.txt
done
for i in $(ls *_species_column.txt);do
	paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column*}_species_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_species_taxa.txt
done
for i in $(ls *_unique_reads.txt);do
	paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%*_order_unique_reads*}_species_taxa.txt) > ${i%_uniq*}_order_unique_uncultured.txt
done
rm *_species_taxid.txt && rm *_uniq_inter.txt && rm *_uniq.txt && rm *_species_column.txt && rm *_species_taxa.txt
for i in $(ls *_order_unique_uncultured.txt);do
	awk -F '\t' '!/Uncultured/' $i > ${i%*_order_unique_uncultured*}_order_unique_sequences.txt
done
rm *_order_unique_uncultured.txt &&
for i in $(ls *_order_OTU.txt);do
	(
	cat $i ${i%_order_OTU*}_order_order_unique_sequences.txt > ${i%_order_OTU*}_complete_order_reads.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done 
wait
for i in $(ls *_complete_order_reads.txt);do
	awk -F '\t' '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' $i > ${i%_complete_order_reads*}_sighits_temp.txt
	echo $'abundance\tqseqid\tsseqid\tlength\tmismatch\tevalue\tpident\tqcovs\tqseq\tsseq\tstaxids\tstitle\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tdomain' | \
	cat - ${i%_complete_order_reads*}_sighits_temp.txt > ${i%_complete_order_reads*}_sighits.txt
done
rm *_complete_order_reads.txt && rm *_sighits_temp.txt && rm *_order_OTU.txt && rm *_unique_reads.txt
}
if [ "$order_level" == "TRUE" ]; then
	time order_filter 2>> $proj_dir/log.out
fi
order_refseq_filter(){
cd $proj_dir/metagenome/sighits/sighits/order
mkdir order_gene_annotation
for i in $(ls *_sighits.txt);do
	cp $i $proj_dir/metagenome/sighits/sighits_order/order_gene_annotation
done
cd $proj_dir/metagenome/sighits/sighis_order/order_gene_annotation
for i in $(ls *_sighits.txt);do
	awk -F '\t' -v var="${i%*_sighits.txt}" '{print $1"\t"var"_"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20}' $i > ${i%*_sighits.txt}_refseq_sighits.txt
done
rm *_sighits.txt
for i in $(ls *_refseq_sighits.txt);do
	awk 'FNR==NR{ h=$0 }NR>1{ print (!a[$16]++? h ORS $0 : $0) >> $16".txt" }' $i
done

for i in $(ls *.txt);do
	sed -i '/^$/d' $i 
done
for i in $(ls *_refseq_sighits.txt);do
	mv $i $proj_dir/metagenome/sighits/sighits_order
done
for i in $(ls *.txt);do
	awk -F '\t' '{print ">"$2"\n"$9}' $i > ${i%.txt}.fasta
done
rm *.txt
for i in $(ls *.fasta);do
	$tool_dir/ncbi-blast-2.8.1+/bin/blastn -task megablast -query $i -db $refseq_db_dir -num_threads $threads -evalue 1e-10 -max_target_seqs 1 -outfmt "6 qseqid sseqid length mismatch evalue pident qcovs qseq sseq staxids stitle" -out ${i%*.fasta}_refseq.megablast
done
for i in $(ls *_refseq.megablast);do
	sort -u $i > ${i%*_refseq.megablast}_refseq_nr.txt
done
rm *_refseq.megablast && rm *.fasta
find . -type f -name '*_refseq_nr.txt' -exec cat {} + > refseq_hits.txt
awk -F '\t' '{print $10}' OFS=';' refseq_hits.txt > taxids_refseq_inter.txt
awk -F ';' '{print $1}' OFS='\t' taxids_refseq_inter.txt > taxids_refseq.txt
rm *_refseq_inter.txt
for i in $(ls taxids_refseq.txt);do
	(
	awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' $tool_dir/rankedlineage_edited.dmp OFS='\t' $i> ${i%taxids_refseq*}refseq_inter.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi 
done
paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' OFS='\t' refseq_hits.txt ) <(awk -F '\t' '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' refseq_inter.txt) > refseq_hits_genus.txt
rm refseq_inter.txt && rm taxids_refseq.txt && rm refseq_hits.txt
cd $proj_dir/metagenome/sighits/sighits_order/
for i in $(ls *_refseq_sighits.txt);do
	mv $i $proj_dir/metagenome/sighits/sighits_order/order_gene_annotation
done
cd $proj_dir/metagenome/sighits/sighits_order/order_gene_annotation
find . -type f -name '*_refseq_sighits.txt' -exec cat {} + > refseq_prehits_temp.txt
awk -F '\t' '!/abundance/' refseq_prehits_temp.txt > refseq_prehits.txt && rm refseq_prehits_temp.txt
awk -F '\t' '{print $2, $14}' OFS='\t' refseq_prehits.txt > prehits.txt
awk -F '\t' '{print $1, $14}' OFS='\t' refseq_hits_genus.txt > hits.txt
awk -F '\t' 'NR==FNR{c[$1$2]++;next};c[$1$2] > 0' prehits.txt hits.txt > consistent_hits_taxids.txt
awk -F '\t' '{print $1}' consistent_hits_taxids.txt > consistent_hits.txt && rm consistent_hits_taxids.txt
for i in $(ls *_refseq_sighits.txt);do
	awk -F '\t' 'NR==FNR{c[$1]++;next};c[$2] > 0' consistent_hits.txt $i > ${i%*_refseq_sighits.txt}_refseq_filtered.txt
done
rm consistent_hits.txt && rm hits.txt && rm prehits.txt && rm refseq_prehits.txt && rm *_refseq_sighits.txt && rm *refseq_hits_order
for i in $(ls *_refseq_filtered.txt);do
	mv $i $proj_dir/metagenome/sighits/sighits_order/
done 
cd $proj_dir/metagenome/sighits/sighits_order/
for i in $(ls *refseq_filtered.txt);do
	echo $'abundance\tqseqid\tsseqid\tlength\tmismatch\tevalue\tpident\tqcovs\tqseq\tsseq\tstaxids\tstitle\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tdomain' | \
	cat - $i > ${i%*_refseq_filtered.txt}_sighits.txt
done
rm *_refseq_filtered.txt
}
if [ "$refseq_filter" == "TRUE" ] && [ "$order_level" == "TRUE" ]; then
	time order_refseq_filter 2>> $proj_dir/log.out
fi


order_stats(){
cd $proj_dir/metagenome/sighits/sighits_order
echo -e "${YELLOW}- compiling taxonomic information"
find . -type f -name '*_sighits.txt' -exec cat {} + > sighits.txt
awk -F '\t' '{print $16"\t"$17"\t"$18"\t"$19"\t"$20}' sighits.txt > rankedlineage_subhits_temp.txt && rm sighits.txt
awk -F '\t' '!/order/' rankedlineage_subhits_temp.txt > rankedlineage_subhits_temp2.txt && rankedlineage_subhits_temp.txt
echo $'order\tclass\tphylum\tkingdom\tdomain' | cat - rankedlineage_subhits_temp2.txt > rankedlineage_subhits.txt && rankedlineage_subhits_temp2.txt


cd $proj_dir/metagenome/sighits/sighits_order
echo -e "${YELLOW}- quantifying the order-level taxonomy"
awk '{print $1}' rankedlineage_subhits.txt | sort -u > order_taxa_mean_temp1.txt
echo -e 'order' | cat - order_taxa_mean_temp1.txt > order_taxa_mean_temp.txt && rm order_taxa_mean_temp1.txt
awk '{print $1}' rankedlineage_subhits.txt | sort -u > order_taxa_unique_sequences_temp1.txt
echo -e 'order' | cat - order_taxa_unique_sequences_temp1.txt > order_taxa_unique_sequences_temp.txt && rm order_taxa_unique_sequences_temp1.txt
awk '{print $1}' rankedlineage_subhits.txt | sort -u > order_taxa_quantification_accuracy_temp1.txt
echo -e 'order' | cat - order_taxa_quantification_accuracy_temp1.txt > order_taxa_quantification_accuracy_temp.txt && rm order_taxa_quantification_accuracy_temp1.txt
order_level=order
for i in $(ls *_sighits.txt);do
	Rscript $tool_dir/Rscripts/stats_summary.R $i $min_uniq $order_level
	echo $'order\tmean\tuniq_reads\tstderr' | cat - stats1.txt > stats2.txt 
	id=${i%_sighits*}_mean && awk -v id=$id '{gsub(/mean/,id); print }' stats2.txt | awk '{print $1,"\t",$2}' > holdmean.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdmean.txt order_taxa_mean_temp.txt > holdmean2.txt && cat holdmean2.txt > order_taxa_mean_temp.txt
	id=${i%_sighits*}_uniq_reads && awk -v id=$id '{gsub(/uniq_reads/,id); print }' stats2.txt | awk '{print $1,"\t",$3}' > holduniq_reads.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holduniq_reads.txt order_taxa_unique_sequences_temp.txt > holduniq_reads2.txt && cat holduniq_reads2.txt > order_taxa_unique_sequences_temp.txt
	id=${i%_sighits*}_stderr && awk -v id=$id '{gsub(/stderr/,id); print }' stats2.txt | awk '{print $1,"\t",$4}' > holdstderr.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdstderr.txt order_taxa_quantification_accuracy_temp.txt > holdstderr2.txt && cat holdstderr2.txt > order_taxa_quantification_accuracy_temp.txt
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $1,"\t",$2,"\t",$3,"\t",$4,"\t", a[$1]}'  rankedlineage_subhits.txt stats2.txt > stats3.txt
	awk '{print $1,$2,$3,$4,$6,$7,$8}' stats3.txt | awk '{gsub(/ /,"\t"); print }' > ${i%_sighits*}_taxastats.txt
	rm *stats1* *stats2* *stats3* *hold*
done
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' order_taxa_mean_temp.txt > order_taxa_mean_temp2.txt && rm order_taxa_mean_temp.txt
awk '{gsub(/ /,"\t"); print}' order_taxa_mean_temp2.txt > ../../results/order_level/order_taxa_mean.txt && rm order_taxa_mean_temp2.txt
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' order_taxa_unique_sequences_temp.txt > order_taxa_unique_sequences_temp2.txt && rm order_taxa_unique_sequences_temp.txt
awk '{gsub(/ /,"\t"); print}' order_taxa_unique_sequences_temp2.txt > ../../results/order_level/order_taxa_unique_sequences.txt && rm order_taxa_unique_sequences_temp2.txt
awk '{gsub(/ /,"\t"); print}' order_taxa_quantification_accuracy_temp.txt > order_taxa_quantification_accuracy_temp2.txt && rm order_taxa_quantification_accuracy_temp.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/order_level/order_taxa_mean.txt order_taxa_quantification_accuracy_temp2.txt > order_taxa_quantification_accuracy_temp3.txt && rm order_taxa_quantification_accuracy_temp2.txt
sed '2,${/order/d;}' order_taxa_quantification_accuracy_temp3.txt > ../../results/order_level/order_taxa_quantification_accuracy.txt && rm *order_taxa_quantification_accuracy_temp3.txt

cd $proj_dir/metagenome/results/order_level
i="_mean$"
awk -vp="$i" 'NR==1{for(i=1; i<=NF; i++) if ($i~p) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' order_taxa_mean.txt > temp_mean.txt
i="uniq_reads$"
awk -vp="$i" 'NR==1{for(i=1; i<=NF; i++) if ($i~p) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' order_taxa_unique_sequences.txt > temp_uniq.txt
paste temp_mean.txt temp_uniq.txt | awk '/^[0-9]/ {for(i=1; i<=NF/2; i++) {s=s OFS $i*$(NF/2+i); }sub(/^ /,x,s);$0=s; s=""} !/[0-9]/{$0=$1;}1' > temp_uniq_mean.txt
tail -n +2 temp_uniq_mean.txt > temp_uniq_mean_2.txt
awk '{for (i=1; i<=NF; i++) sum[i]+=$i;}; END{for (i in sum) print sum[i];}' temp_uniq_mean_2.txt > temp_uniq_mean_3.txt
awk '{sum+=$1}END{print sum}' temp_uniq_mean_3.txt > mean_uniq.txt
paste mean_uniq.txt seqcov.txt | awk '{print(($1/$2)* 100)}' > order_percent_coverage.txt
rm *temp* *mean_uniq*
cd $proj_dir/metagenome/results/order_level
for i in {mean,unique_sequences,quantification_accuracy}; do
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_order/rankedlineage_subhits.txt order_taxa_${i}.txt | \
	awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' > order_taxainfo_${i}.txt
done
rm *_taxa_*
echo -e "${YELLOW}- creating order-level visualizations"
order_taxa_mean=order_taxainfo_mean.txt
order_taxa_uniq=order_taxainfo_unique_sequences.txt
order_taxa_quant=order_taxainfo_quantification_accuracy.txt
order_mean=order_mean.txt
order_uniq=order_unique_sequences.txt
order_quant=order_quantification_accuracy.txt
percent_thresh=5
Rscript $tool_dir/Rscripts/order_level_corr.R $order_mean $percent_thresh &>/dev/null
}
if [ "$order_level" == "TRUE" ]; then
	time order_stats 2>> $proj_dir/log.out
fi
################################################################################################################
class_filter(){
cd $proj_dir/metagenome/sighits
mkdir sighits_class
cd $proj_dir/metagenome/results
mkdir class_level
cp seqcov.txt ./class_level
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is Preforming Class-Level Classification \n\e[97m########################################################\n"
cd $proj_dir/metagenome/alignment
echo -e "${YELLOW}- preforming exact-matching algorithm"
for i in $(ls *_haplotig_nd.megablast);do
	awk '$6>=97' $i | awk 'gsub(" ","_",$0)' > ../sighits/sighits_class/${i%_haplotig*}_filter.txt
	awk 'NR==FNR {h[$1] = $2; next} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,h[$1]}' ../haplotig/${i%_haplotig*}_normalized.txt ../sighits/sighits_class/${i%_haplotig*}_filter.txt | 
	awk '{print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' | awk 'gsub(" ","\t",$0)' > ../sighits/sighits_class/${i%_haplotig*}_sighits.txt
done
rm ../sighits/sighits_class/*_filter.txt 
echo -e "${YELLOW}- compiling class-level OTUs with multi-alignment algorithm"
cd $proj_dir/metagenome/sighits/sighits_class/
for i in $(ls *_sighits.txt);do
	awk -F '\t' '{a[$2]++;b[$2]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i > ${i%_sighits*}_class_unique_reads.txt
done
cd $proj_dir/metagenome/sighits/sighits_class
for i in $(ls *_sighits.txt);do
		awk -F'|' 'FNR==NR{a[$1,$2]=1; next}  !a[$1,$2]' ${i%_sighits.txt}_class_unique_reads.txt $i > ${i%_sighits*}_dup.txt
done
cd $proj_dir/metagenome/sighits/sighits_class
for i in $(ls *_dup.txt);do
	awk -F '\t' '{print $11}' OFS=';' $i > ${i%_dup*}_taxids_dup_inter.txt
done

for i in $(ls *_dup_inter.txt);do
	awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_dup_inter*}_taxids_dup.txt
done

rm *_taxids_dup_inter.txt
cd $proj_dir/metagenome/sighits/sighits_class
for i in $(ls *_taxids_dup.txt);do
	(
	awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' /home/brandon/Desktop/Qmatey/tools/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_dup*}_dup_inter.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi 
done
wait
for i in $(ls *_dup_inter.txt);do
	(
	awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_dup_inter*}_species_taxid.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
rm *_taxids_dup.txt
for i in $(ls *_species_taxid.txt);do
	awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_species_taxid*}_species_column.txt
done
for i in $(ls *_species_column.txt);do
	paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column*}_species_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_species_taxa.txt
done
for i in $(ls *_dup.txt);do
	paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' OFS='\t' ${i%*_dup.txt}_species_taxa.txt) > ${i%_dup*}_class_duplicates_uncultured.txt
done
rm *_species_taxid.txt && rm *_dup_inter.txt && rm *_dup.txt && rm *_species_column.txt && rm *_species_taxa.txt
for i in $(ls *_class_duplicates_uncultured.txt);do
	awk -F '\t' '!/Uncultured/' $i > ${i%*_class_duplicates_uncultured*}_class_duplicates.txt
done
rm *_class_duplicates_uncultured.txt
for i in $(ls *_class_duplicates.txt);do
	(
	awk -F '\t' '{print $1, $2"~"$18, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_class_duplicates*}_class_inter.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
for i in $(ls *_class_inter.txt);do
	(
	awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_class_inter*}_class_inter2.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
for i in $(ls *_class_inter2.txt);do
	(
	awk -F '\t' '{dups[$1]++} END {for (num in dups) {print num}}' $i | sort -k1,1  > ${i%_class_inter2*}_duplicate_count.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
for i in $(ls *_duplicate_count.txt);do
	(
	awk -F '~' '{a[$1]++;b[$1]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i | sort -k1,1 > ${i%_duplicate_count*}_multialign_class_reads.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
for i in $(ls *_multialign_class_reads.txt);do
	(
	awk -F '\t'  'FNR==NR {a[$1]; next}; $1 in a' $i ${i%_multialign_class_reads*}_class_inter2.txt | sort -u -k1,1 | awk 'gsub("~","\t",$0)'| awk -F '\t' '{print $3, $1, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' > ${i%_multialign_class_reads*}_class_OTU.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
rm *_class_inter.txt && rm *_class_inter2.txt && rm *_duplicate_count.txt && rm *_multialign_class_reads.txt && rm *_class_duplicates.txt
for i in $(ls *_class_unique_reads.txt);do
	awk -F '\t' '{print $11}' OFS=';' $i > ${i%_class_unique_reads*}_taxids_uniq_inter.txt
done

for i in $(ls *_uniq_inter.txt);do
	awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_uniq_inter*}_taxids_uniq.txt
done

rm *_taxids_uniq_inter.txt
cd $proj_dir/metagenome/sighits/sighits_class
for i in $(ls *_taxids_uniq.txt);do
	(
	awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' $tool_dir/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_uniq*}_uniq_inter.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi 
done
wait
for i in $(ls *_uniq_inter.txt);do
	(
	awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_uniq_inter*}_species_taxid.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
rm *_taxids_uniq.txt
for i in $(ls *_species_taxid.txt);do
	awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_species_taxid*}_species_column.txt
done
for i in $(ls *_species_column.txt);do
	paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column*}_species_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_species_taxa.txt
done
for i in $(ls *_unique_reads.txt);do
	paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%*_class_unique_reads*}_species_taxa.txt) > ${i%_uniq*}_class_unique_uncultured.txt
done
rm *_species_taxid.txt && rm *_uniq_inter.txt && rm *_uniq.txt && rm *_species_column.txt && rm *_species_taxa.txt
for i in $(ls *_class_unique_uncultured.txt);do
	awk -F '\t' '!/Uncultured/' $i > ${i%*_class_unique_uncultured*}_class_unique_sequences.txt
done
for i in $(ls *_class_OTU.txt);do
	(
	cat $i ${i%_class_OTU*}_class_class_unique_sequences.txt > ${i%_class_OTU*}_complete_class_reads.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done 
wait
for i in $(ls *_complete_class_reads.txt);do
	awk -F '\t' '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' $i > ${i%_complete_class_reads*}_sighits_temp.txt
	echo $'abundance\tqseqid\tsseqid\tlength\tmismatch\tevalue\tpident\tqcovs\tqseq\tsseq\tstaxids\tstitle\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tdomain' | \
	cat - ${i%_complete_class_reads*}_sighits_temp.txt > ${i%_complete_class_reads*}_sighits.txt
done
rm *_complete_class_reads.txt && rm *_sighits_temp.txt && rm *_class_OTU.txt && rm *_unique_reads.txt
}
if [ "$class_level" == "TRUE" ]; then
	time class_filter 2>> $proj_dir/log.out
fi 
class_refseq_filter(){
cd $proj_dir/metagenome/sighits/sighits/class
mkdir class_gene_annotation
for i in $(ls *_sighits.txt);do
	cp $i $proj_dir/metagenome/sighits/sighits_class/class_gene_annotation
done
cd $proj_dir/metagenome/sighits/sighis_class/class_gene_annotation
for i in $(ls *_sighits.txt);do
	awk -F '\t' -v var="${i%*_sighits.txt}" '{print $1"\t"var"_"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20}' $i > ${i%*_sighits.txt}_refseq_sighits.txt
done
rm *_sighits.txt
for i in $(ls *_refseq_sighits.txt);do
	awk 'FNR==NR{ h=$0 }NR>1{ print (!a[$17]++? h ORS $0 : $0) >> $17".txt" }' $i
done

for i in $(ls *.txt);do
	sed -i '/^$/d' $i 
done
for i in $(ls *_refseq_sighits.txt);do
	mv $i $proj_dir/metagenome/sighits/sighits_class
done
for i in $(ls *.txt);do
	awk -F '\t' '{print ">"$2"\n"$9}' $i > ${i%.txt}.fasta
done
rm *.txt
for i in $(ls *.fasta);do
	$tool_dir/ncbi-blast-2.8.1+/bin/blastn -task megablast -query $i -db $refseq_db_dir -num_threads $threads -evalue 1e-10 -max_target_seqs 1 -outfmt "6 qseqid sseqid length mismatch evalue pident qcovs qseq sseq staxids stitle" -out ${i%*.fasta}_refseq.megablast
done
for i in $(ls *_refseq.megablast);do
	sort -u $i > ${i%*_refseq.megablast}_refseq_nr.txt
done
rm *_refseq.megablast && rm *.fasta
find . -type f -name '*_refseq_nr.txt' -exec cat {} + > refseq_hits.txt
awk -F '\t' '{print $10}' OFS=';' refseq_hits.txt > taxids_refseq_inter.txt
awk -F ';' '{print $1}' OFS='\t' taxids_refseq_inter.txt > taxids_refseq.txt
rm *_refseq_inter.txt
for i in $(ls taxids_refseq.txt);do
	(
	awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' $tool_dir/rankedlineage_edited.dmp OFS='\t' $i> ${i%taxids_refseq*}refseq_inter.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi 
done
paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' OFS='\t' refseq_hits.txt ) <(awk -F '\t' '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' refseq_inter.txt) > refseq_hits_genus.txt
rm refseq_inter.txt && rm taxids_refseq.txt && rm refseq_hits.txt
cd $proj_dir/metagenome/sighits/sighits_class/
for i in $(ls *_refseq_sighits.txt);do
	mv $i $proj_dir/metagenome/sighits/sighits_class/class_gene_annotation
done
cd $proj_dir/metagenome/sighits/sighits_class/class_gene_annotation
find . -type f -name '*_refseq_sighits.txt' -exec cat {} + > refseq_prehits_temp.txt
awk -F '\t' '!/abundance/' refseq_prehits_temp.txt > refseq_prehits.txt && rm refseq_prehits_temp.txt
awk -F '\t' '{print $2, $14}' OFS='\t' refseq_prehits.txt > prehits.txt
awk -F '\t' '{print $1, $14}' OFS='\t' refseq_hits_genus.txt > hits.txt
awk -F '\t' 'NR==FNR{c[$1$2]++;next};c[$1$2] > 0' prehits.txt hits.txt > consistent_hits_taxids.txt
awk -F '\t' '{print $1}' consistent_hits_taxids.txt > consistent_hits.txt && rm consistent_hits_taxids.txt
for i in $(ls *_refseq_sighits.txt);do
	awk -F '\t' 'NR==FNR{c[$1]++;next};c[$2] > 0' consistent_hits.txt $i > ${i%*_refseq_sighits.txt}_refseq_filtered.txt
done
rm consistent_hits.txt && rm hits.txt && rm prehits.txt && rm refseq_prehits.txt && rm *_refseq_sighits.txt && rm *refseq_hits_class
for i in $(ls *_refseq_filtered.txt);do
	mv $i $proj_dir/metagenome/sighits/sighits_class/
done 
cd $proj_dir/metagenome/sighits/sighits_class/
for i in $(ls *refseq_filtered.txt);do
	echo $'abundance\tqseqid\tsseqid\tlength\tmismatch\tevalue\tpident\tqcovs\tqseq\tsseq\tstaxids\tstitle\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tdomain' | \
	cat - $i > ${i%*_refseq_filtered.txt}_sighits.txt
done
rm *_refseq_filtered.txt
}
if [ "$refseq_filter" == "TRUE" ] && [ "$class_level" == "TRUE" ]; then
	time class_refseq_filter 2>> $proj_dir/log.out
fi
class_stats(){
cd $proj_dir/metagenome/sighits/sighits_class
echo -e "${YELLOW}- compiling taxonomic information"
find . -type f -name '*_sighits.txt' -exec cat {} + > sighits.txt
awk -F '\t' '{print $17"\t"$18"\t"$19"\t"$20}' sighits.txt > rankedlineage_subhits_temp.txt && rm sighits.txt
awk -F '\t' '!/class/' rankedlineage_subhits_temp.txt > rankedlineage_subhits_temp2.txt && rankedlineage_subhits_temp.txt
echo $'class\tphylum\tkingdom\tdomain' | cat - rankedlineage_subhits_temp2.txt > rankedlineage_subhits.txt && rankedlineage_subhits_temp2.txt


cd $proj_dir/metagenome/sighits/sighits_class
echo -e "${YELLOW}- quantifying the class-level taxonomy"
awk '{print $1}' rankedlineage_subhits.txt | sort -u > class_taxa_mean_temp1.txt
echo -e 'class' | cat - class_taxa_mean_temp1.txt > class_taxa_mean_temp.txt && rm class_taxa_mean_temp1.txt
awk '{print $1}' rankedlineage_subhits.txt | sort -u > class_taxa_unique_sequences_temp1.txt
echo -e 'class' | cat - class_taxa_unique_sequences_temp1.txt > class_taxa_unique_sequences_temp.txt && rm class_taxa_unique_sequences_temp1.txt
awk '{print $1}' rankedlineage_subhits.txt | sort -u > class_taxa_quantification_accuracy_temp1.txt
echo -e 'class' | cat - class_taxa_quantification_accuracy_temp1.txt > class_taxa_quantification_accuracy_temp.txt && rm class_taxa_quantification_accuracy_temp1.txt
class_level=class
for i in $(ls *_sighits.txt);do
	Rscript $tool_dir/Rscripts/stats_summary.R $i $min_uniq $class_level
	echo $'class\tmean\tuniq_reads\tstderr' | cat - stats1.txt > stats2.txt 
	id=${i%_sighits*}_mean && awk -v id=$id '{gsub(/mean/,id); print }' stats2.txt | awk '{print $1,"\t",$2}' > holdmean.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdmean.txt class_taxa_mean_temp.txt > holdmean2.txt && cat holdmean2.txt > class_taxa_mean_temp.txt
	id=${i%_sighits*}_uniq_reads && awk -v id=$id '{gsub(/uniq_reads/,id); print }' stats2.txt | awk '{print $1,"\t",$3}' > holduniq_reads.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holduniq_reads.txt class_taxa_unique_sequences_temp.txt > holduniq_reads2.txt && cat holduniq_reads2.txt > class_taxa_unique_sequences_temp.txt
	id=${i%_sighits*}_stderr && awk -v id=$id '{gsub(/stderr/,id); print }' stats2.txt | awk '{print $1,"\t",$4}' > holdstderr.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdstderr.txt class_taxa_quantification_accuracy_temp.txt > holdstderr2.txt && cat holdstderr2.txt > class_taxa_quantification_accuracy_temp.txt
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $1,"\t",$2,"\t",$3,"\t",$4,"\t", a[$1]}'  rankedlineage_subhits.txt stats2.txt > stats3.txt
	awk '{print $1,$2,$3,$4,$6,$7}' stats3.txt | awk '{gsub(/ /,"\t"); print }' > ${i%_sighits*}_taxastats.txt
	rm *stats1* *stats2* *stats3* *hold*
done
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' class_taxa_mean_temp.txt > class_taxa_mean_temp2.txt && rm class_taxa_mean_temp.txt
awk '{gsub(/ /,"\t"); print}' class_taxa_mean_temp2.txt > ../../results/class_level/class_taxa_mean.txt && rm class_taxa_mean_temp2.txt
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' class_taxa_unique_sequences_temp.txt > class_taxa_unique_sequences_temp2.txt && rm class_taxa_unique_sequences_temp.txt
awk '{gsub(/ /,"\t"); print}' class_taxa_unique_sequences_temp2.txt > ../../results/class_level/class_taxa_unique_sequences.txt && rm class_taxa_unique_sequences_temp2.txt
awk '{gsub(/ /,"\t"); print}' class_taxa_quantification_accuracy_temp.txt > class_taxa_quantification_accuracy_temp2.txt && rm class_taxa_quantification_accuracy_temp.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/class_level/class_taxa_mean.txt class_taxa_quantification_accuracy_temp2.txt > class_taxa_quantification_accuracy_temp3.txt && rm class_taxa_quantification_accuracy_temp2.txt
sed '2,${/class/d;}' class_taxa_quantification_accuracy_temp3.txt > ../../results/class_level/class_taxa_quantification_accuracy.txt && rm class_taxa_quantification_accuracy_temp3.txt

cd $proj_dir/metagenome/results/class_level
i="_mean$"
awk -vp="$i" 'NR==1{for(i=1; i<=NF; i++) if ($i~p) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' class_taxa_mean.txt > temp_mean.txt
i="uniq_reads$"
awk -vp="$i" 'NR==1{for(i=1; i<=NF; i++) if ($i~p) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' class_taxa_unique_sequences.txt > temp_uniq.txt
paste temp_mean.txt temp_uniq.txt | awk '/^[0-9]/ {for(i=1; i<=NF/2; i++) {s=s OFS $i*$(NF/2+i); }sub(/^ /,x,s);$0=s; s=""} !/[0-9]/{$0=$1;}1' > temp_uniq_mean.txt
tail -n +2 temp_uniq_mean.txt > temp_uniq_mean_2.txt
awk '{for (i=1; i<=NF; i++) sum[i]+=$i;}; END{for (i in sum) print sum[i];}' temp_uniq_mean_2.txt > temp_uniq_mean_3.txt
awk '{sum+=$1}END{print sum}' temp_uniq_mean_3.txt > mean_uniq.txt
paste mean_uniq.txt seqcov.txt | awk '{print(($1/$2)* 100)}' > class_percent_coverage.txt
rm *temp* *mean_uniq*
cd $proj_dir/metagenome/results/class_level
for i in {mean,unique_sequences,quantification_accuracy}; do
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_class/rankedlineage_subhits.txt class_taxa_${i}.txt | \
	awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' > class_taxainfo_${i}.txt
done
rm *_taxa_*
echo -e "${YELLOW}- creating class-level visualizations"
class_taxa_mean=class_taxainfo_mean.txt
class_taxa_uniq=class_taxainfo_unique_sequences.txt
class_taxa_quant=class_taxainfo_quantification_accuracy.txt
class_mean=class_mean.txt
class_uniq=class_unique_sequences.txt
class_quant=class_quantification_accuracy.txt
percent_thresh=5
Rscript $tool_dir/Rscripts/class_level_corr.R $class_mean $percent_thresh &>/dev/null
}
if [ "$class_level" == "TRUE" ]; then
	time class_stats 2>> $proj_dir/log.out
fi

phylum_filter(){
cd $proj_dir/metagenome/sighits
mkdir sighits_phylum
cd $proj_dir/metagenome/results
mkdir phylum_level
cp seqcov.txt ./phylum_level
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is Preforming Phylum-Level Classification \n\e[97m########################################################\n"
cd $proj_dir/metagenome/alignment
echo -e "${YELLOW}- preforming exact-matching algorithm"
for i in $(ls *_haplotig_nd.megablast);do
	awk '$6>=97' $i | awk 'gsub(" ","_",$0)' > ../sighits/sighits_phylum/${i%_haplotig*}_filter.txt
	awk 'NR==FNR {h[$1] = $2; next} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,h[$1]}' ../haplotig/${i%_haplotig*}_normalized.txt ../sighits/sighits_phylum/${i%_haplotig*}_filter.txt | 
	awk '{print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' | awk 'gsub(" ","\t",$0)' > ../sighits/sighits_phylum/${i%_haplotig*}_sighits.txt
done
rm ../sighits/sighits_phylum/*_filter.txt 
echo -e "${YELLOW}- compiling phylum-level OTUs with multi-alignment algorithm"
cd $proj_dir/metagenome/sighits/sighits_phylum/
for i in $(ls *_sighits.txt);do
	awk -F '\t' '{a[$2]++;b[$2]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i > ${i%_sighits*}_phylum_unique_reads.txt
done
cd $proj_dir/metagenome/sighits/sighits_phylum
for i in $(ls *_sighits.txt);do
		awk -F'|' 'FNR==NR{a[$1,$2]=1; next}  !a[$1,$2]' ${i%_sighits.txt}_phylum_unique_reads.txt $i > ${i%_sighits*}_dup.txt
done
cd $proj_dir/metagenome/sighits/sighits_phylum
for i in $(ls *_dup.txt);do
	awk -F '\t' '{print $11}' OFS=';' $i > ${i%_dup*}_taxids_dup_inter.txt
done

for i in $(ls *_dup_inter.txt);do
	awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_dup_inter*}_taxids_dup.txt
done

rm *_taxids_dup_inter.txt
cd $proj_dir/metagenome/sighits/sighits_phylum
for i in $(ls *_taxids_dup.txt);do
	(
	awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' /home/brandon/Desktop/Qmatey/tools/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_dup*}_dup_inter.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi 
done
wait
for i in $(ls *_dup_inter.txt);do
	(
	awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_dup_inter*}_species_taxid.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
rm *_taxids_dup.txt
for i in $(ls *_species_taxid.txt);do
	awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_species_taxid*}_species_column.txt
done
for i in $(ls *_species_column.txt);do
	paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column*}_species_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_species_taxa.txt
done
for i in $(ls *_dup.txt);do
	paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' OFS='\t' ${i%*_dup.txt}_species_taxa.txt) > ${i%_dup*}_phylum_duplicates_uncultured.txt
done
rm *_species_taxid.txt && rm *_dup_inter.txt && rm *_dup.txt && rm *_species_column.txt && rm *_species_taxa.txt
for i in $(ls *_phylum_duplicates_uncultured.txt);do
	awk -F '\t' '!/Uncultured/' $i > ${i%*_phylum_duplicates_uncultured*}_phylum_duplicates.txt
done
rm *_phylum_duplicates_uncultured.txt
for i in $(ls *_phylum_duplicates.txt);do
	(
	awk -F '\t' '{print $1, $2"~"$19, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_phylum_duplicates*}_phylum_inter.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
for i in $(ls *_phylum_inter.txt);do
	(
	awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_phylum_inter*}_phylum_inter2.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
for i in $(ls *_phylum_inter2.txt);do
	(
	awk -F '\t' '{dups[$1]++} END {for (num in dups) {print num}}' $i | sort -k1,1  > ${i%_phylum_inter2*}_duplicate_count.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
for i in $(ls *_duplicate_count.txt);do
	(
	awk -F '~' '{a[$1]++;b[$1]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i | sort -k1,1 > ${i%_duplicate_count*}_multialign_phylum_reads.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
for i in $(ls *_multialign_phylum_reads.txt);do
	(
	awk -F '\t'  'FNR==NR {a[$1]; next}; $1 in a' $i ${i%_multialign_phylum_reads*}_phylum_inter2.txt | sort -u -k1,1 | awk 'gsub("~","\t",$0)'| awk -F '\t' '{print $3, $1, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' > ${i%_multialign_phylum_reads*}_phylum_OTU.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
rm *_phylum_inter.txt && rm *_phylum_inter2.txt && rm *_duplicate_count.txt && rm *_multialign_phylum_reads.txt && rm *_phylum_duplicates.txt
for i in $(ls *_phylum_unique_reads.txt);do
	awk -F '\t' '{print $11}' OFS=';' $i > ${i%_phylum_unique_reads*}_taxids_uniq_inter.txt
done

for i in $(ls *_uniq_inter.txt);do
	awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_uniq_inter*}_taxids_uniq.txt
done

rm *_taxids_uniq_inter.txt
cd $proj_dir/metagenome/sighits/sighits_phylum
for i in $(ls *_taxids_uniq.txt);do
	(
	awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' $tool_dir/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_uniq*}_uniq_inter.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi 
done
wait
for i in $(ls *_uniq_inter.txt);do
	(
	awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_uniq_inter*}_species_taxid.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done
wait
rm *_taxids_uniq.txt
for i in $(ls *_species_taxid.txt);do
	awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_species_taxid*}_species_column.txt
done
for i in $(ls *_species_column.txt);do
	paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column*}_species_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_species_taxa.txt
done
for i in $(ls *_unique_reads.txt);do
	paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%*_phylum_unique_reads*}_species_taxa.txt) > ${i%_uniq*}_phylum_unique_uncultured.txt
done
rm *_species_taxid.txt && rm *_uniq_inter.txt && rm *_species_column.txt && rm *_species_taxa.txt
for i in $(ls *_phylum_unique_uncultured.txt);do
	awk -F '\t' '!/Uncultured/' $i > ${i%*_phylum_unique_uncultured*}_phylum_unique_sequences.txt
done

for i in $(ls *_phylum_OTU.txt);do
	(
	cat $i ${i%_phylum_OTU*}_phylum_phylum_unique_sequences.txt > ${i%_phylum_OTU*}_complete_phylum_reads.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi
done 
wait
for i in $(ls *_complete_phylum_reads.txt);do
	awk -F '\t' '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' $i > ${i%_complete_phylum_reads*}_sighits_temp.txt
	echo $'abundance\tqseqid\tsseqid\tlength\tmismatch\tevalue\tpident\tqcovs\tqseq\tsseq\tstaxids\tstitle\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tdomain' | \
	cat - ${i%_complete_phylum_reads*}_sighits_temp.txt > ${i%_complete_phylum_reads*}_sighits.txt
done
rm *_complete_phylum_reads.txt && rm *_sighits_temp.txt && rm *_phylum_OTU.txt && rm *_unique_reads.txt
cd $proj_dir/metagenome/sighits/sighits_phylum
echo -e "${YELLOW}- compiling taxonomic information"
find . -type f -name '*_sighits.txt' -exec cat {} + > sighits.txt
awk -F '\t' '{print $18"\t"$19"\t"$20}' sighits.txt > rankedlineage_subhits_temp.txt && rm sighits.txt
awk -F '\t' '!/phylum/' rankedlineage_subhits_temp.txt > rankedlineage_subhits_temp2.txt && rankedlineage_subhits_temp.txt
echo $'phylum\tkingdom\tdomain' | cat - rankedlineage_subhits_temp2.txt > rankedlineage_subhits.txt && rankedlineage_subhits_temp2.txt


cd $proj_dir/metagenome/sighits/sighits_phylum
echo -e "${YELLOW}- quantifying the phylum-level taxonomy"
awk '{print $1}' rankedlineage_subhits.txt | sort -u > phylum_taxa_mean_temp1.txt
echo -e 'phylum' | cat - phylum_taxa_mean_temp1.txt > phylum_taxa_mean_temp.txt && rm phylum_taxa_mean_temp1.txt
awk '{print $1}' rankedlineage_subhits.txt | sort -u > phylum_taxa_unique_sequences_temp1.txt
echo -e 'phylum' | cat - phylum_taxa_unique_sequences_temp1.txt > phylum_taxa_unique_sequences_temp.txt && rm phylum_taxa_unique_sequences_temp1.txt
awk '{print $1}' rankedlineage_subhits.txt | sort -u > phylum_taxa_quantification_accuracy_temp1.txt
echo -e 'phylum' | cat - phylum_taxa_quantification_accuracy_temp1.txt > phylum_taxa_quantification_accuracy_temp.txt && rm phylum_taxa_quantification_accuracy_temp1.txt
}
if [ "$phylum_level" == "TRUE" ]; then
	time phylum_filter 2>> $proj_dir/log.out
fi
phylum_refseq_filter(){
cd $proj_dir/metagenome/sighits/sighits/phylum
mkdir phylum_gene_annotation
for i in $(ls *_sighits.txt);do
	cp $i $proj_dir/metagenome/sighits/sighits_phylum/phylum_gene_annotation
done
cd $proj_dir/metagenome/sighits/sighis_phylum/phylum_gene_annotation
for i in $(ls *_sighits.txt);do
	awk -F '\t' -v var="${i%*_sighits.txt}" '{print $1"\t"var"_"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20}' $i > ${i%*_sighits.txt}_refseq_sighits.txt
done
rm *_sighits.txt
for i in $(ls *_refseq_sighits.txt);do
	awk 'FNR==NR{ h=$0 }NR>1{ print (!a[$18]++? h ORS $0 : $0) >> $18".txt" }' $i
done

for i in $(ls *.txt);do
	sed -i '/^$/d' $i 
done
for i in $(ls *_refseq_sighits.txt);do
	mv $i $proj_dir/metagenome/sighits/sighits_phylum
done
for i in $(ls *.txt);do
	awk -F '\t' '{print ">"$2"\n"$9}' $i > ${i%.txt}.fasta
done
rm *.txt
for i in $(ls *.fasta);do
	$tool_dir/ncbi-blast-2.8.1+/bin/blastn -task megablast -query $i -db $refseq_db_dir -num_threads $threads -evalue 1e-10 -max_target_seqs 1 -outfmt "6 qseqid sseqid length mismatch evalue pident qcovs qseq sseq staxids stitle" -out ${i%*.fasta}_refseq.megablast
done
for i in $(ls *_refseq.megablast);do
	sort -u $i > ${i%*_refseq.megablast}_refseq_nr.txt
done
rm *_refseq.megablast && rm *.fasta
find . -type f -name '*_refseq_nr.txt' -exec cat {} + > refseq_hits.txt
awk -F '\t' '{print $10}' OFS=';' refseq_hits.txt > taxids_refseq_inter.txt
awk -F ';' '{print $1}' OFS='\t' taxids_refseq_inter.txt > taxids_refseq.txt
rm *_refseq_inter.txt
for i in $(ls taxids_refseq.txt);do
	(
	awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' $tool_dir/rankedlineage_edited.dmp OFS='\t' $i> ${i%taxids_refseq*}refseq_inter.txt
	) &
	if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
	wait
	fi 
done
paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}' OFS='\t' refseq_hits.txt ) <(awk -F '\t' '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' refseq_inter.txt) > refseq_hits_genus.txt
rm refseq_inter.txt && rm taxids_refseq.txt && rm refseq_hits.txt
cd $proj_dir/metagenome/sighits/sighits_phylum/
for i in $(ls *_refseq_sighits.txt);do
	mv $i $proj_dir/metagenome/sighits/sighits_phylum/phylum_gene_annotation
done
cd $proj_dir/metagenome/sighits/sighits_phylum/phylum_gene_annotation
find . -type f -name '*_refseq_sighits.txt' -exec cat {} + > refseq_prehits_temp.txt
awk -F '\t' '!/abundance/' refseq_prehits_temp.txt > refseq_prehits.txt && rm refseq_prehits_temp.txt
awk -F '\t' '{print $2, $14}' OFS='\t' refseq_prehits.txt > prehits.txt
awk -F '\t' '{print $1, $14}' OFS='\t' refseq_hits_genus.txt > hits.txt
awk -F '\t' 'NR==FNR{c[$1$2]++;next};c[$1$2] > 0' prehits.txt hits.txt > consistent_hits_taxids.txt
awk -F '\t' '{print $1}' consistent_hits_taxids.txt > consistent_hits.txt && rm consistent_hits_taxids.txt
for i in $(ls *_refseq_sighits.txt);do
	awk -F '\t' 'NR==FNR{c[$1]++;next};c[$2] > 0' consistent_hits.txt $i > ${i%*_refseq_sighits.txt}_refseq_filtered.txt
done
rm consistent_hits.txt && rm hits.txt && rm prehits.txt && rm refseq_prehits.txt && rm *_refseq_sighits.txt && rm *refseq_hits_phylum
for i in $(ls *_refseq_filtered.txt);do
	mv $i $proj_dir/metagenome/sighits/sighits_phylum/
done 
cd $proj_dir/metagenome/sighits/sighits_phylum/
for i in $(ls *refseq_filtered.txt);do
	echo $'abundance\tqseqid\tsseqid\tlength\tmismatch\tevalue\tpident\tqcovs\tqseq\tsseq\tstaxids\tstitle\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tdomain' | \
	cat - $i > ${i%*_refseq_filtered.txt}_sighits.txt
done
rm *_refseq_filtered.txt
}
if [ "$refseq_filter" == "TRUE" ] && [ "$phylum_level" == "TRUE" ]; then
	time phylum_refseq_filter 2>> $proj_dir/log.out
fi
phylum_stats(){
phylum_level=phylum
for i in $(ls *_sighits.txt);do
	Rscript $tool_dir/Rscripts/stats_summary.R $i $min_uniq $phylum_level
	echo $'phylum\tmean\tuniq_reads\tstderr' | cat - stats1.txt > stats2.txt 
	id=${i%_sighits*}_mean && awk -v id=$id '{gsub(/mean/,id); print }' stats2.txt | awk '{print $1,"\t",$2}' > holdmean.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdmean.txt phylum_taxa_mean_temp.txt > holdmean2.txt && cat holdmean2.txt > phylum_taxa_mean_temp.txt
	id=${i%_sighits*}_uniq_reads && awk -v id=$id '{gsub(/uniq_reads/,id); print }' stats2.txt | awk '{print $1,"\t",$3}' > holduniq_reads.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holduniq_reads.txt phylum_taxa_unique_sequences_temp.txt > holduniq_reads2.txt && cat holduniq_reads2.txt > phylum_taxa_unique_sequences_temp.txt
	id=${i%_sighits*}_stderr && awk -v id=$id '{gsub(/stderr/,id); print }' stats2.txt | awk '{print $1,"\t",$4}' > holdstderr.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdstderr.txt phylum_taxa_quantification_accuracy_temp.txt > holdstderr2.txt && cat holdstderr2.txt > phylum_taxa_quantification_accuracy_temp.txt
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $1,"\t",$2,"\t",$3,"\t",$4,"\t", a[$1]}'  rankedlineage_subhits.txt stats2.txt > stats3.txt
	awk '{print $1,$2,$3,$4,$6}' stats3.txt | awk '{gsub(/ /,"\t"); print }' > ${i%_sighits*}_taxastats.txt
	rm *stats1* *stats2* *stats3* *hold*
done
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' phylum_taxa_mean_temp.txt > phylum_taxa_mean_temp2.txt && rm phylum_taxa_mean_temp.txt
awk '{gsub(/ /,"\t"); print}' phylum_taxa_mean_temp2.txt > ../../results/phylum_level/phylum_taxa_mean.txt && rm phylum_taxa_mean_temp2.txt
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' phylum_taxa_unique_sequences_temp.txt > phylum_taxa_unique_sequences_temp2.txt && rm phylum_taxa_unique_sequences_temp.txt
awk '{gsub(/ /,"\t"); print}' phylum_taxa_unique_sequences_temp2.txt > ../../results/phylum_level/phylum_taxa_unique_sequences.txt && rm phylum_taxa_unique_sequences_temp2.txt
awk '{gsub(/ /,"\t"); print}' phylum_taxa_quantification_accuracy_temp.txt > phylum_taxa_quantification_accuracy_temp2.txt && rm phylum_taxa_quantification_accuracy_temp.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/phylum_level/phylum_taxa_mean.txt phylum_taxa_quantification_accuracy_temp2.txt > phylum_taxa_quantification_accuracy_temp3.txt && rm phylum_taxa_quantification_accuracy_temp2.txt
sed '2,${/phylum/d;}' phylum_taxa_quantification_accuracy_temp3.txt > ../../results/phylum_level/phylum_taxa_quantification_accuracy.txt && rm phylum_taxa_quantification_accuracy_temp3.txt

cd $proj_dir/metagenome/results/phylum_level
i="_mean$"
awk -vp="$i" 'NR==1{for(i=1; i<=NF; i++) if ($i~p) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' phylum_taxa_mean.txt > temp_mean.txt
i="uniq_reads$"
awk -vp="$i" 'NR==1{for(i=1; i<=NF; i++) if ($i~p) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' phylum_taxa_unique_sequences.txt > temp_uniq.txt
paste temp_mean.txt temp_uniq.txt | awk '/^[0-9]/ {for(i=1; i<=NF/2; i++) {s=s OFS $i*$(NF/2+i); }sub(/^ /,x,s);$0=s; s=""} !/[0-9]/{$0=$1;}1' > temp_uniq_mean.txt
tail -n +2 temp_uniq_mean.txt > temp_uniq_mean_2.txt
awk '{for (i=1; i<=NF; i++) sum[i]+=$i;}; END{for (i in sum) print sum[i];}' temp_uniq_mean_2.txt > temp_uniq_mean_3.txt
awk '{sum+=$1}END{print sum}' temp_uniq_mean_3.txt > mean_uniq.txt
paste mean_uniq.txt seqcov.txt | awk '{print(($1/$2)* 100)}' > phylum_percent_coverage.txt
rm *temp* *mean_uniq*
cd $proj_dir/metagenome/results/phylum_level
for i in {mean,unique_sequences,quantification_accuracy}; do
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_phylum/rankedlineage_subhits.txt phylum_taxa_${i}.txt | \
	awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' > phylum_taxainfo_${i}.txt
done
rm *_taxa_*
echo -e "${YELLOW}- creating phylum-level visualizations"
phylum_taxa_mean=phylum_taxainfo_mean.txt
phylum_taxa_uniq=phylum_taxainfo_unique_sequences.txt
phylum_taxa_quant=phylum_taxainfo_quantification_accuracy.txt
phylum_mean=phylum_mean.txt
phylum_uniq=phylum_unique_sequences.txt
phylum_quant=phylum_quantification_accuracy.txt
percent_thresh=5
Rscript $tool_dir/Rscripts/phylum_level_corr.R $phylum_mean $percent_thresh &>/dev/null
}
if [ "$phylum_level" == "TRUE" ]; then
	time phylum_stats 2>> $proj_dir/log.out
fi
cd $tool_dir
rm rankedlineage.dmp && rm rankedlineage_edited.dmp



















