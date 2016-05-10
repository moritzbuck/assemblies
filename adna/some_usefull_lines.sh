# a loop to run kraken in a loop, first copy your DB 


module load Kraken

MY_KDB=/scratch/kraken_db
mkdir $MY_KDB
cp -av $KRAKEN_DB/* $MY_KDB


for d in `find . -type d | grep clean`
do
echo "kraken $d"
kraken --preload --threads 16 --db $MY_KDB $d/*.fastq  > ${d%%clean}.kraken
done

# running mapDamage in parallel over many samples

parallel -j 16 --xapply mapDamage --merge-reference-sequences -d {1}/mapDamage -i {1}/map_hits.sam  -r ../final.contigs.fa ::: Batch*

# rnuning map damage on a single sample :

mapDamage --merge-reference-sequences -d output_folder -i map_frile_from_bbmap_forexample.sam  -r reference_it_was_mapped_too.fasta

