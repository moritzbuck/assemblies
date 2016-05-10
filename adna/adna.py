from pandas import DataFrame, Index
from os.path import join as pjoin
import subprocess
from Bio import SeqIO
from tqdm import tqdm
import gzip
from generic import Assembly
from generic import Sample
import os
import json
from generic.scripts import *

project = "b2015426"
root = "/proj/b2015426/private/20151202_NJ7P"
outfolder = pjoin(root, "processed/")
ass_path = "/proj/b2015426/private/20151202_NJ7P/processed/assemblies/megahit/"

# find folders containing each sample
folders = subprocess.Popen(["find", root, "-type", "d" ], stdout = subprocess.PIPE)
folders = [f[:-1] for f in folders.stdout if "Sample_" in f and "/PRI_" not in f ]
folders = {"-".join(f.split("_")[-2:]) : f for f in folders}

# loadind and cleaning metadata
samples_metadata = DataFrame.from_csv("adna/samples.csv")
samples_metadata.index = Index([i.replace("(repeat)","").replace("(Repeat)","").replace(" ","") for i in samples_metadata.index])


# do the simple samples
samples = {b : {} for  b in folders.keys()}
for b in folders.keys(): 
    samples[b]["reads"] = [{"U" : pjoin(folders[b],f) } for f in os.listdir(folders[b]) if f[-8:] == "fastq.gz"]


# take care of the samples that bugged (e.g. the sequencing facility did not demultiplex for some reason

buggy_lane = pjoin(root, "20160113_data/ftp.dna.ku.dk/20160105_hiseq4a/lane2/")
buggy_libs = [pjoin(buggy_lane,b) for b in os.listdir(buggy_lane)  if b[-8:] == "fastq.gz"]
buggy_samples = [k for k in list(samples_metadata.index) if k not in folders.keys()]
index2sample = {v :k for k,v in samples_metadata.loc[buggy_samples]['Index sequence'].to_dict().iteritems()}
run = False
if run:
    deplexed = {b : [] for b in buggy_samples}
    for f in  buggy_libs:
        print "processing ", f
        with gzip.open(f) as handle:
            for s in tqdm(SeqIO.parse(handle,"fastq")):
                tag = s.description.split(":")[-1]
                if index2sample.has_key(tag):
                    deplexed[index2sample[tag]] += [s]
            

for b in buggy_samples:
    sample_folder = pjoin(outfolder, "samples", b)
    if run:
        print "writting samples ", b 
        if not os.path.exists(sample_folder):
            os.makedirs(sample_folder)
        with gzip.open(pjoin(sample_folder,b + "_reads.fastq.gz"),"w") as handle:
            SeqIO.write(deplexed[b],handle,"fastq")
#    samples[b] = {}
#    samples[b]['reads'] = [{"U" : pjoin(sample_folder,b + "_reads.fastq.gz")}]
    
assembly_data = {"name": "adna", "out_path" : outfolder, "samples" : samples}
with open(pjoin(outfolder, "assembly_data.json"),"w") as handle:
    json.dump(assembly_data, handle)

assembly = Assembly(pjoin(outfolder, "assembly_data.json"), "b2015426")
all_shufed =  [{ "U" : pjoin(outfolder, "samples", f) } for f in os.listdir(pjoin(outfolder, "samples/")) if "all_trimmed" in f]
combined_sample = Sample(name = "shuffled", reads = all_shufed , path = pjoin(outfolder, "samples/"), no_cleans = True)




if run : 
    assembly.clean_all(submit = True)
    """
for d in `find . -type d | grep clean`
do
echo "kraken $d"
kraken --preload --threads 16 --db $MY_KDB $d/*.fastq  > ${d%%clean}.kraken
done
"""
    assembly.megahit_assembly(name = "megahit", reads = ['/proj/b2015426/private/20151202_NJ7P/processed/samples/shuffled/normalized/shuffled_unpaired_normalised.fastq'], submit = True)

    # mapping trimmed reads:
    reference = pjoin(ass_path, "final.contigs.fa")
    for samp in assembly.samples:
        reads = [pjoin(samp.path,p) for p in os.listdir(samp.path) if ".fastq.gz" in p][0]
        path = pjoin(os.path.dirname(reference), "mapping_"+ os.path.basename(reference)  , samp.name )
        script = make_slurm_header(path, "mapping_to_" + samp.name, project,  time = "1-00:00:00")
        script += make_single_bbmapping_script(reference, path, reads, only_bowtie = False)

        with open(pjoin(path,  "mapping_script.sh"),"w") as handle:
            handle.writelines(script)
    
        sh.sbatch(pjoin(path,  "mapping_script.sh"))

    for samp in assembly.samples:
        reads = ".gz ".join(sum(samp.get_singles().values(),[]))
        path = pjoin(os.path.dirname(reference), "untrimmed_map_"+ os.path.basename(reference)  , samp.name )
        script = make_slurm_header(path, "untrimmed_map_" + samp.name, project,  time = "1-00:00:00")
        script += make_single_bbmapping_script(reference, path, reads, only_bowtie = True)

        with open(pjoin(path,  "mapping_script.sh"),"w") as handle:
            handle.writelines(script)
    
        sh.sbatch(pjoin(path,  "mapping_script.sh"))

    
    for samp in assembly.samples:
        reads = ".gz ".join(sum(samp.get_singles().values(),[]))
        path = pjoin(os.path.dirname(reference), "untrimmed_map_"+ os.path.basename(reference)  , samp.name )
        script = make_slurm_header(path, "untrimmed_map_" + samp.name, project,  time = "1-00:00:00")
        script += make_single_bbmapping_script(reference, path, reads, only_bowtie = True)

        with open(pjoin(path,  "mapping_script.sh"),"w") as handle:
            handle.writelines(script)

    reference = "/home/moritz/glob/data/genomes/Ptrichocarpa_210_v3.0.hardmasked.fa"
    
    for samp in assembly.samples:
        reads = ".gz ".join(sum(samp.get_singles().values(),[]))
        path = pjoin(os.path.dirname(reference), "untrimmed_map_"+ os.path.basename(reference)  , samp.name )
        script = make_slurm_header(path, "untrimmed_map_" + samp.name, project,  time = "1-00:00:00")
        script += make_single_bbmapping_script(reference, path, reads, only_bowtie = True)
        with open(pjoin(path,  "mapping_script.sh"),"w") as handle:
            handle.writelines(script)

