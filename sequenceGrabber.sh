# Fetch genomes from NCBI for addition to GFKB.
# Very conservative (hard coded) wait time, not suitable for very large numbers of accessions.

while read p; do
  sleep 10s & efetch -db nuccore -id $p -format fasta; wait
done <accessions.txt
