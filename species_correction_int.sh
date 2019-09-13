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

################################################################################################################
echo -e "${YELLOW}------------------------------------------------------------------------------ \n \n Qmatey is performing species-level clustering \n \n------------------------------------------------------------------------------"
#Species-level sighit identification
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
	awk -F '\t' '{print $11}' OFS='\t' $i > ${i%_dup*}_taxids_dup.txt
done

################################################################################################################
cd $proj_dir/metagenome/sighits/sighits_species
#N=$threads
for i in $(ls *_taxids_dup.txt);do
	#((t=t%N)); ((t++==0)) && wait
	awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' /home/brandon/Desktop/Qmatey/tools/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_dup*}_dup_inter.txt #&
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
	paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' OFS='\t' ${i%*_dup.txt}_species_taxa.txt) > ${i%_dup*}_species_duplicates.txt
done



rm *_species_taxid.txt && rm *_dup_inter.txt && rm *_dup.txt && rm *_species_column.txt && rm *_species_taxa.txt
#################################################################################################################
#species-level clustering
#N=$threads
cd /home/brandon/Desktop/Qmatey/examples/project1/metagenome/sighits/sighits_species/
for i in $(ls *_species_duplicates.txt);do
	#((t=t%N)); ((t++==0)) && wait
	awk -F '\t' '{a[$2,$14]++;b[$2,$14]=$0}END{for(x in a)if(a[x]==1)print b[x]}' OFS='\t' $i > ${i%_species_duplicates*}_species_errors.txt
done

k=$(ls *_species_duplicates.txt)
for i in $(ls *_species_errors.txt);do
	awk -F '\t' 'FNR==NR{a[$2]=1; next}  !a[$2]' OFS='\t' $i ${i%_species_errors.txt}_species_duplicates.txt > ${i%_species_errors*}_filtered_species_sequences.txt
done

for i in $(ls *_filtered_species_sequences.txt);do
	awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i > ${i%_filtered_species_sequences*}_multi_species_reads.txt
done

rm *_filtered_species_sequences.txt

for i in $(ls *_multi_species_reads.txt);do
	cat $i ${i%_multi_species_reads.txt}_species_unique_reads.txt > ${i%_multi_species_reads*}_complete_species_reads.txt
done

rm *_multi_species_reads.txt && rm *_species_duplicates.txt && rm *_species_errors.txt && rm *_species_unique_reads.txt 


for i in $(ls *_complete_species_reads.txt);do
	awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' $i > ${i%_complete_species_reads*}_sighits_temp.txt
	echo $'abundance\tqseqid\tsseqid\tlength\tmismatch\tevalue\tpident\tqcovs\tqseq\tsseq\tstaxids\tstitle' | \
	cat - ${i%_complete_species_reads*}_sighits_temp.txt > ${i%_complete_species_reads*}_sighits.txt
done

rm *_complete_species_reads.txt && rm *_sighits_temp.txt
#################################################################################################################
#Combine all taxids for all files/individuals and perform single search against new_taxdump.
cd $proj_dir/metagenome/sighits/sighits_species
find . -type f -name '*_sighits.txt' -exec cat {} + > sighits.txt
awk '{print $11}' sighits.txt | awk '{gsub(";","\n"); print}' | sort -u -n | sed -e '1s/staxids/tax_id/' > taxids_sighits.txt && rm sighits.txt
awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}'  /home/brandon/Desktop/Qmatey/tools/rankedlineage_edited.dmp taxids_sighits.txt | \
awk '{gsub(/ /,"_"); print }' > rankedlineage_subhits.txt 

rm taxids_sighits.txt
#################################################################################################################
#Install datamash
#Generate file with mean, number of unique reads per taxID, and standard error
#Now, perform merge subsetted new_taxdump with each file while retaining only taxids in each file
cd $proj_dir/metagenome/sighits/sighits_species/
awk '{print $1}' rankedlineage_subhits.txt > ../../results/species_level/proj_taxa_mean.txt
awk '{print $1}' rankedlineage_subhits.txt > ../../results/species_level/proj_taxa_uniq_reads.txt
awk '{print $1}' rankedlineage_subhits.txt > ../../results/species_level/proj_taxa_stderr.txt 
for i in $(ls *_sighits.txt); do
	cut -f 1,11 $i | awk '{print $2,"\t",$1}' | datamash --header-in --sort --group 1 mean 2 sstdev 2 count 2 | \
	awk '{ print $1,"\t",$2,"\t",$4,"\t",$3 / sqrt($4) }' > stats1.txt
	echo $'tax_id\tmean\tuniq_reads\tstderr' | cat - stats1.txt > stats2.txt
	id=${i%_sighits*}_mean && awk -v id=$id '{gsub(/mean/,id); print }' stats2.txt | awk '{print $1,"\t",$2}' > holdmean.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdmean.txt ../../results/species_level/proj_taxa_mean.txt > holdmean2.txt && cat holdmean2.txt > ../../results/species_level/proj_taxa_mean.txt
	id=${i%_sighits*}_uniq_reads && awk -v id=$id '{gsub(/uniq_reads/,id); print }' stats2.txt | awk '{print $1,"\t",$3}' > holduniq_reads.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holduniq_reads.txt ../../results/species_level/proj_taxa_uniq_reads.txt > holduniq_reads2.txt && cat holduniq_reads2.txt > ../../results/species_level/proj_taxa_uniq_reads.txt
	id=${i%_sighits*}_stderr && awk -v id=$id '{gsub(/stderr/,id); print }' stats2.txt | awk '{print $1,"\t",$4}' > holdstderr.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdstderr.txt ../../results/species_level/proj_taxa_stderr.txt > holdstderr2.txt && cat holdstderr2.txt > ../../results/species_level/proj_taxa_stderr.txt
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $1,"\t",$2,"\t",$3,"\t",$4,"\t", a[$1]}'  rankedlineage_subhits.txt stats2.txt > stats3.txt
	awk '{print $1,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14 }' stats3.txt | awk '{gsub(/ /,"\t"); print }' > ${i%_sighits*}_taxastats.txt
	rm *stats1* *stats2* *stats3* *hold*
done
#################################################################################################################
#Calculates percent coverage of sighits from unmatched microbiome input files
cd $proj_dir/metagenome/results/species_level/
i="_mean$"
awk -vp="$i" 'NR==1{for(i=1; i<=NF; i++) if ($i~p) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' proj_taxa_mean.txt > temp_mean.txt
i="uniq_reads$"
awk -vp="$i" 'NR==1{for(i=1; i<=NF; i++) if ($i~p) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' proj_taxa_uniq_reads.txt > temp_uniq.txt
paste temp_mean.txt temp_uniq.txt | awk '/^[0-9]/ {for(i=1; i<=NF/2; i++) {s=s OFS $i*$(NF/2+i); }sub(/^ /,x,s);$0=s; s=""} !/[0-9]/{$0=$1;}1' > temp_uniq_mean.txt
tail -n +2 temp_uniq_mean.txt > temp_uniq_mean_2.txt
awk '{for (i=1; i<=NF; i++) sum[i]+=$i;}; END{for (i in sum) print sum[i];}' temp_uniq_mean_2.txt > temp_uniq_mean_3.txt
awk '{sum+=$1}END{print sum}' temp_uniq_mean_3.txt > mean_uniq.txt
paste mean_uniq.txt seqcov.txt | awk '{print(($1/$2)* 100)}' > proj_percent_coverage.txt
rm *temp* *mean_uniq*
################################################################################################################
cd $proj_dir/metagenome/results/species_level
for i in {mean,uniq_reads,stderr}; do
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_species/rankedlineage_subhits.txt proj_taxa_${i}.txt | \
	awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' > proj_taxainfo_${i}.txt
done
rm *_taxa_*
################################################################################################################
#Genus-level sighit identification
echo -e "${YELLOW}------------------------------------------------------------------------------ \n \n Qmatey is performing genus-level clustering \n \n------------------------------------------------------------------------------"
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
	awk '{print $11}' $i > ${i%_dup*}_taxids_dup.txt
done
################################################################################################################
cd /home/brandon/Desktop/Qmatey/examples/project1/metagenome/sighits/sighits_genus/
#N=$threads
for i in $(ls *_taxids_dup.txt);do
	#((t=t%N)); ((t++==0)) && wait
	awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' /home/brandon/Desktop/Qmatey/tools/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_dup*}_dup_inter.txt #&
done


for i in $(ls *_dup_inter.txt);do
	awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_dup_inter*}_genus_taxid.txt
done

rm *_taxids_dup.txt

for i in $(ls *_genus_taxid.txt);do
	awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_genus_taxid*}_species_column.txt
done

for i in $(ls *_species_column.txt);do
	paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column.txt}_genus_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_genus_taxa.txt
done

for i in $(ls *_dup.txt);do
	paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' OFS='\t' ${i%_dup.txt}_genus_taxa.txt) > ${i%_dup*}_genus_duplicates.txt
done

rm *_genus_taxid.txt && rm *_dup_inter.txt && rm *_dup.txt && rm *_species_column.txt && rm *_genus_taxa.txt 
################################################################################################################
#Genus-level clustering
cd /home/brandon/Desktop/Qmatey/examples/project1/metagenome/sighits/sighits_genus/
for i in $(ls *_genus_duplicates.txt);do
	#((t=t%N)); ((t++==0)) && wait
	awk -F '\t' '{a[$2,$15]++;b[$2,$15]=$0}END{for(x in a)if(a[x]==1)print b[x]}' OFS='\t' $i > ${i%_genus_duplicates*}_genus_errors.txt
done

for i in $(ls *_genus_errors.txt);do
	awk -F '\t' 'FNR==NR{a[$2]=1; next}  !a[$2]' OFS='\t' $i ${i%_genus_errors.txt}_genus_duplicates.txt > ${i%_genus_errors*}_filtered_genus_sequences.txt
done

for i in $(ls *_filtered_genus_sequences.txt);do
	awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i > ${i%_filtered_genus_sequences*}_multi_genus_reads.txt
done

rm *_filtered_genus_sequences.txt

for i in $(ls *_multi_genus_reads.txt);do
	cat $i ${i%_multi_genus_reads.txt}_genus_unique_reads.txt > ${i%_multi_genus_reads*}_complete_genus_reads.txt
done
rm *_multi_genus_reads.txt && rm *_genus_duplicates.txt && rm *_genus_unique_reads.txt && rm *_genus_errors.txt

for i in $(ls *_complete_genus_reads.txt);do
	awk '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' $i > ${i%_complete_genus_reads*}_sighits_temp.txt
	echo $'abundance\tqseqid\tsseqid\tlength\tmismatch\tevalue\tpident\tqcovs\tqseq\tsseq\tstaxids\tstitle' | \
	cat - ${i%_complete_genus_reads*}_sighits_temp.txt > ${i%_complete_genus_reads*}_sighits.txt
done

rm *_complete_genus_reads.txt && rm *_sighits_temp.txt
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
awk '{print $1}' rankedlineage_subhits.txt > ../../results/genus_level/proj_taxa_mean.txt
awk '{print $1}' rankedlineage_subhits.txt > ../../results/genus_level/proj_taxa_uniq_reads.txt
awk '{print $1}' rankedlineage_subhits.txt > ../../results/genus_level/proj_taxa_stderr.txt 
for i in $(ls *_sighits.txt); do
	cut -f 1,11 $i | awk '{print $2,"\t",$1}' | datamash --header-in --sort --group 1 mean 2 sstdev 2 count 2 | \
	awk '{ print $1,"\t",$2,"\t",$4,"\t",$3 / sqrt($4) }' > stats1.txt
	echo $'tax_id\tmean\tuniq_reads\tstderr' | cat - stats1.txt > stats2.txt
	id=${i%_sighits*}_mean && awk -v id=$id '{gsub(/mean/,id); print }' stats2.txt | awk '{print $1,"\t",$2}' > holdmean.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdmean.txt ../../results/genus_level/proj_taxa_mean.txt > holdmean2.txt && cat holdmean2.txt > ../../results/genus_level/proj_taxa_mean.txt
	id=${i%_sighits*}_uniq_reads && awk -v id=$id '{gsub(/uniq_reads/,id); print }' stats2.txt | awk '{print $1,"\t",$3}' > holduniq_reads.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holduniq_reads.txt ../../results/genus_level/proj_taxa_uniq_reads.txt > holduniq_reads2.txt && cat holduniq_reads2.txt > ../../results/genus_level/proj_taxa_uniq_reads.txt
	id=${i%_sighits*}_stderr && awk -v id=$id '{gsub(/stderr/,id); print }' stats2.txt | awk '{print $1,"\t",$4}' > holdstderr.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdstderr.txt ../../results/genus_level/proj_taxa_stderr.txt > holdstderr2.txt && cat holdstderr2.txt > ../../results/genus_level/proj_taxa_stderr.txt
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $1,"\t",$2,"\t",$3,"\t",$4,"\t", a[$1]}'  rankedlineage_subhits.txt stats2.txt > stats3.txt
	awk '{print $1,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14 }' stats3.txt | awk '{gsub(/ /,"\t"); print }' > ${i%_sighits*}_taxastats.txt
	rm *stats1* *stats2* *stats3* *hold*
done
#################################################################################################################
#Calculates percent coverage of sighits from unmatched microbiome input files
cd $proj_dir/metagenome/results/genus_level
i="_mean$"
awk -vp="$i" 'NR==1{for(i=1; i<=NF; i++) if ($i~p) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' proj_taxa_mean.txt > temp_mean.txt
i="uniq_reads$"
awk -vp="$i" 'NR==1{for(i=1; i<=NF; i++) if ($i~p) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' proj_taxa_uniq_reads.txt > temp_uniq.txt
paste temp_mean.txt temp_uniq.txt | awk '/^[0-9]/ {for(i=1; i<=NF/2; i++) {s=s OFS $i*$(NF/2+i); }sub(/^ /,x,s);$0=s; s=""} !/[0-9]/{$0=$1;}1' > temp_uniq_mean.txt
tail -n +2 temp_uniq_mean.txt > temp_uniq_mean_2.txt
awk '{for (i=1; i<=NF; i++) sum[i]+=$i;}; END{for (i in sum) print sum[i];}' temp_uniq_mean_2.txt > temp_uniq_mean_3.txt
awk '{sum+=$1}END{print sum}' temp_uniq_mean_3.txt > mean_uniq.txt
paste mean_uniq.txt seqcov.txt | awk '{print(($1/$2)* 100)}' > proj_percent_coverage.txt
rm *temp* *mean_uniq*
################################################################################################################
cd $proj_dir/metagenome/results/genus_level
for i in {mean,uniq_reads,stderr}; do
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_genus/rankedlineage_subhits.txt proj_taxa_${i}.txt | \
	awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' > proj_taxainfo_${i}.txt
done
rm *_taxa_*
################################################################################################################


