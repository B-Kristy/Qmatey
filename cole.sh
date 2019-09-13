cd /home/brandon/Desktop/Qmatey/examples/project1/metagenome/sighits/sighits_species
for i in $(ls *_filtered_species_sequences.txt);do
	grep 'Metazoa' $i > ${i%_filtered_species_sequences*}_incorrect_reads.txt 
done

for i in $(ls *_incorrect_reads.txt);do
	awk -F '\t' '{print $3}' $i | awk -F "|" '{print $2}' > ${i%_incorrect_reads*}_masked_gi_list.gi
done
$k= $(ls *_masked_gi_list.gi)

for i in $(ls *_incorrect_reads.txt);do
	awk '{print " " $2 "\t" $10}' test_input_incorrect_reads.txt | awk 'gsub(" ", ">")' > test_input_incorrect_sequence_id.txt
done


for i in $(ls *_incorrect_sequence_id.txt);do
	awk '{print $1"\n"$2}' test_input_incorrect_sequence_id.txt > test_input_incorrect_sequence.fasta
done

for i in $(ls *_incorrect_sequence.fasta);do
	/home/brandon/Desktop/Qmatey/tools/ncbi-blast-2.8.1+/bin/blastn -task megablast -query $i -db /home/brandon/Desktop/db_dir/nt.00 -num_threads 1 -evalue 1e-10 -negative_gilist test_input_masked_gi_list.gi -max_target_seqs 5 -outfmt \
	"6 qseqid sseqid length mismatch evalue pident qcovs qseq sseq staxids stitle" \
	-out ${i%_incorrect_fasta*}_incorrect.megablast
done
