#Before species-level clustering
awk 'NR>1 gsub($3, $2)' metagenome1.txt > metagenome2.txt
awk '{print $3}' metagenome2.txt > metagenome3.txt
sed 's/^.*_//' metagenome3.txt > metagenome4.txt
paste <(awk '{print $1, $2}' metagenome2.txt) <(awk '{print $1}' metagenome4.txt ) > metagenome5.txt



