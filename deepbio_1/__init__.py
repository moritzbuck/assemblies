from __future__ import division

import sh
import os
import json
from tqdm import tqdm
from pandas import DataFrame

#from plumbing.autopaths import AutoPaths

raw_slurm = """#!/bin/bash
#SBATCH -D %s
#SBATCH -J ray_deepbio_1_%s
#SBATCH -o %sass.out
#SBATCH -e %sass.err
#SBATCH -A b2011032
#SBATCH -t 2-00:00:00
#SBATCH -n 16
#SBATCH -p node
#SBATCH -C mem256GB
#SBATCH --mail-user murumbii@gmail.com
#SBATCH --mail-type=ALL
"""

raw_map_slurm = """#!/bin/bash
#SBATCH -D %s
#SBATCH -J map_%s_%s
#SBATCH -o %smap.out
#SBATCH -e %smap.err
#SBATCH -A b2011032
#SBATCH -t 16:00:00
#SBATCH -n 16
#SBATCH -p node
#SBATCH --mail-user murumbii@gmail.com
#SBATCH --mail-type=ALL

"""

raw_newbler_slurm = """#!/bin/bash
#SBATCH -D %s
#SBATCH -J newbler_%s
#SBATCH -o %snewbler.out
#SBATCH -e %snewbler.err
#SBATCH -A b2011105
#SBATCH -t 1-00:00:00
#SBATCH -n 16
#SBATCH -p node
#SBATCH --mail-user murumbii@gmail.com
#SBATCH --mail-type=ALL

"""

raw_concoct_slurm = """#!/bin/bash
#SBATCH -D %s
#SBATCH -J concoct_%s
#SBATCH -o %sconcoct.out
#SBATCH -e %sconcoct.err
#SBATCH -A b2014318
#SBATCH -t 19:00:00
#SBATCH -n 16
#SBATCH -p node
#SBATCH --mail-user murumbii@gmail.com
#SBATCH --mail-type=ALL

"""

raw_concoct = """
concoct --coverage_file  %s --composition_file %s -l %s -b %s --no_total_coverage
"""

raw_checkm = """
checkm lineage_wf -t 16 -x fasta . %s/checkm > %scheckm.txt
"""            


raw_hmm_pipe = """#!/bin/bash
#SBATCH -D %s
#SBATCH -J hmmer_pipe
#SBATCH -o %shmmer.out
#SBATCH -e %shmmer.err
#SBATCH -A b2014318
#SBATCH -t 18:00:00
#SBATCH -n 16
#SBATCH -p node
#SBATCH --mail-user murumbii@gmail.com
#SBATCH --mail-type=ALL

prokka --metagenome --prefix %sannotation/ --locustag metagenome  %s454AllContigs.fna

parallel --xapply -j 16 hmmersearch --max {1}  %s/annotation/metagenome.faa ::: `ls ~/glob/data/pfam/pfam-a/*.hmm` 
""" 


parallel_raw_map = """
#bowtie2-build %s
parallel --xapply -j 16 /home/moritz/repos/metassemble/scripts/map/map-bowtie2-markduplicates.sh -ct1 {1}_R1_001.clean.fastq {1}_R2_001.clean.fastq {2} %s map map_{2}  '>' map_{2}.out  '2>' map_{2}.err ::: %s ::: %s
"""

raw_map = "module load BEDTools\n/home/moritz/repos/metassemble/scripts/map/map-bowtie2-markduplicates_half.sh -ct 16 %s %s %s %s %s"
raw_unmap = "bowtie2 -ct 16 -p 16 -x %s -1 %s -2 %s  --un-conc %s > /dev/null" 


raw_ray = """
k=`cat %scurrent_k.txt`
mpiexec -n 16 Ray231 -k $k -o %s_$k %s
k=$((k + 10))
echo $k > %scurrent_k.txt
if [ `cat current_k.txt` -le 91 ]; then sbatch %s ; else echo "thats all folks"; fi
"""


raw_newbler = """
/home/moritz/repos/metassemble/scripts/assembly/merge-asm-newbler.sh %s %s
"""

raw_path =  "/home/moritz/GEFES/views/pools/run005-pool01/clean/"
nr_threads = 16
executable = 'Ray231'
out_root = "/home/moritz/people/valerie/deepbio/deepbio_1/assemblies/"
stats_out = out_root + "stats/"
raw_ass = out_root + "raw_assemblies/"
merged_ass = out_root + "merged_assemblies/"
newbler = sh.Command("/home/moritz/repos/metassemble/scripts/assembly/merge-asm-newbler.sh")
minimus = sh.Command("/home/moritz/repos/metassemble/scripts/assembly/merge-asm-minimus2.sh")

class Assembly(object):

    def __repr__(self): return '<%s object with k-mer length %i>' % \
        (self.__class__.__name__, -1 if len(self.success_ks) == 0 else max(self.success_ks) )
 
    
    def __init__(self, k=41):
        self.dir = raw_path
        self.coass = False
        self.pair = sum([[ self.dir + d for d in os.listdir(self.dir) if ".fastq" in d and "fwd" in d], [self.dir + d for d in os.listdir(self.dir) if ".fastq" in d and "rev" in d]],[])

        self.k = k
        self.out_dir = raw_ass 
        self.merge_dir = raw_ass +  "merged/" 

        if not os.path.exists(raw_ass):
            os.makedirs(raw_ass)
        if not os.path.exists(stats_out):
            os.makedirs(stats_out)
        self.success_ks = [int(f.split("_")[1]) for f in os.listdir(raw_ass) if "ray" in f]
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        if os.path.exists(self.out_dir + "_stats.json"):
            print "Loading stats" 
            with open(self.out_dir + "_stats.json") as handle:
                self.stats = json.load(handle)
         

    def run(self):
        self.arg =  "-p " + " ".join(self.pair)
        self.script = raw_slurm % ( self.out_dir, str(self.k), self.out_dir, self.out_dir)
            
        self.script = self.script + raw_ray % (self.out_dir, self.out_dir + "ray", self.arg, self.out_dir, self.out_dir + "slurm_script.sh")
        with open(self.out_dir + "current_k.txt","w") as handle:
            handle.writelines([str(self.k)])                    
        with open(self.out_dir + "slurm_script.sh","w") as handle:
            handle.writelines(self.script)
        #sh.sbatch(self.out_dir + "slurm_script.sh")

    def make_stats(self, rerun = False):
        self.stats = {}
        if os.path.exists(self.out_dir + "mapper/"):
            if os.path.exists(self.out_dir + self.name + "_stats.json") and not rerun:
                print "Loading stats for", self.name
                with open(self.out_dir + self.name + "_stats.json") as handle:
                    self.stats = json.load(handle)
            else:
                print "Computing stats for", self.name
                self.stats['raw_reads'] = self.pairs
                self.stats['output'] = self.out_dir
                print "Computing raw number of paired-reads"
                self.stats['read_count'] = int(sh.zgrep("-Ec", "$", self.reads))/4 
                print "Computing clean number of paired-reads"
                self.stats['clean_read_count'] = int(sh.grep("-Ec", "$", self.reads.replace("fastq.gz","clean.fastq")))/4
                self.stats['ratios'] = float(self.stats['clean_read_count'])/float(self.stats['read_count'])
                for k in tqdm(self.success_ks):
                    self.stats[k] = {}
                    self.stats[k]['Success'] =  os.path.exists(self.out_dir + "ray_" + str(k) + "/Contigs.fasta")
                    if self.stats[k]['Success']:
                        t=sh.assemstats(0,self.out_dir + "ray_" + str(k) + "/Contigs.fasta" )
                        self.stats[k]['Assembly'] = {a:b for a,b in  zip(*[[w.strip() for w in  l.split("\t")][1:] for l in str(t).split("\n")[0:2]])}
                        self.stats[k]['mapped_frac'] = float(sh.grep("overall",  self.out_dir + "map_" + str(k) + ".err").split(" ")[0][0:-1])/100
                        self.stats[k]['dupli_frac'] = float(sh.grep("Unknown", self.out_dir + "mapper/ray_" +str(k)+ "_" + self.name+ "-smd.metrics").split("\t")[7])
            
        with open(self.out_dir + self.name + "_stats.json", 'w') as handle:
            json.dump(self.stats, handle)
        
    def map_all(self):
        for k in self.success_ks:
            self.map(k)
                    
    def map_merged(self):
        if os.path.exists(self.merge_dir + "454AllContigs.fna"):
            if not os.path.exists(self.out_dir + "mapper_" + "merged" + "/"):
                os.makedirs(self.out_dir + "mapper_" + "merged" + "/")

            slurm_header = raw_map_slurm % ( self.out_dir, "deepbio_1", "merged", self.out_dir + "mapper_" + "merged" + "/", self.out_dir + "mapper_" + "merged" + "/")
            mapper_cmd = raw_map % (" ".join(self.pair), "deepbio_1", self.merge_dir + "454AllContigs.fna" , "merged" , self.out_dir + "mapper_" + "merged" + "/")
            self.map_script = slurm_header + mapper_cmd
        
            with open(self.out_dir + "mapper_" + "merged" + "/" + "map_slurm_script_" + "merged" + ".sh","w") as handle:
                handle.writelines(self.map_script)
            sh.sbatch(self.out_dir + "mapper_" + "merged" + "/" + "map_slurm_script_" + "merged" + ".sh")
            print "Launched mapper for", "merged deepbio_1"

    def map(self,k):
        if os.path.exists(self.out_dir + "ray_" + str(k) + "/Contigs.fasta"):
            if not os.path.exists(self.out_dir + "mapper_" + str(k) + "/"):
                os.makedirs(self.out_dir + "mapper_" + str(k) + "/")

            slurm_header = raw_map_slurm % ( self.out_dir, "deepbio_1", str(k), self.out_dir + "mapper_" + str(k) + "/", self.out_dir + "mapper_" + str(k) + "/")
            mapper_cmd = raw_map % (" ".join(self.pair), "deepbio_1", self.out_dir + "ray_" + str(k) + "/Contigs.fasta", "ray_" + str(self.k), self.out_dir + "mapper_" + str(k) + "/")
            self.map_script = slurm_header + mapper_cmd
        
            with open(self.out_dir + "mapper_" + str(k) + "/" + "map_slurm_script_" + str(k) + ".sh","w") as handle:
                handle.writelines(self.map_script)
            sh.sbatch(self.out_dir + "mapper_" + str(k) + "/" + "map_slurm_script_" + str(k) + ".sh")
            print "Launched mapper for", "deepbio_1", "with k-mer size", k

    def map2ref(self,ref):
        slurm_header = raw_map_slurm % ( self.out_dir, self.name, "2ref", self.out_dir + "mapper_" + "2ref" + "/", self.out_dir + "mapper_" + "2ref" + "/")
        mapper_cmd = raw_map % (self.pair, "deepbio_1", ref, "ray_" + "2ref", self.out_dir + "mapper_" + "2ref" + "/")
        self.map_script = slurm_header + mapper_cmd
        
        with open(self.out_dir + "mapper_" + "2ref" + "/" + "map_slurm_script_" + "2ref" + ".sh","w") as handle:
            handle.writelines(self.map_script)
        sh.sbatch(self.out_dir + "mapper_" + "2ref" + "/" + "map_slurm_script_" + "2ref" + ".sh")
        print "Launched mapper for", self.name, "with reference", ref

            
    def format_lib_stats(self):
        ff = lambda x :  str(round(x/1000000,2))+"M"
        pp = lambda x :  str(round(x*100,2))+"%"
        return  "\t".join([self.name, ff(self.stats['pairs_count']), ff(self.stats['clean_pairs_count'])],pp(self.stats['ratios']))

    def format_assem_stats(self):
        ff = lambda x :  str(round(float(x)/1000000,2))+"M"
        pp = lambda x :  str(round(x*100,2))+"%"
        kk = lambda x :  str(round(float(x)/1000,2))+"k"
        lines = []
        for k in tqdm(self.success_ks):
            if self.stats[k]['Success']:
                line = "\t".join([self.name, str(k),ff(self.stats[k]['Assembly']['sum']), ff(self.stats[k]['Assembly']['n']), kk(self.stats[k]['Assembly']['max']),kk(self.stats[k]['Assembly']['n50_len']), pp(self.stats[k]['map_stats']['dupli_frac']), pp(self.stats[k]['map_stats']['mapped_frac']), self.stats[k]['output'] + "ray_" + k + "/Contigs.fasta" ])
            else:
                line = "\t".join([self.name, str(k) ,"NA" ,"NA","NA","NA","NA","NA","NA"])
        lines.append(line)
        return lines

    def merge(self):
        fastas = [self.out_dir + "ray_" + str(k) + "/Contigs.fasta" for k in self.success_ks if os.path.exists(self.out_dir + "ray_" + str(k) + "/Contigs.fasta")]
        if not os.path.exists(self.merge_dir) : os.makedirs(self.merge_dir)
        in_size = sum([os.stat(f).st_size for f in fastas])
        newbler = raw_newbler % (self.merge_dir, " ".join(fastas))
        header = raw_newbler_slurm % ( self.merge_dir, "deepbio_1", self.merge_dir, self.merge_dir)
        with open(self.merge_dir + "merge_slurm_script.sh","w") as handle:
            handle.writelines(header + newbler)
    #        sh.sbatch(self.merge_dir + "merge_slurm_script.sh")

    def concoct_the_hell(self):
        if not os.path.exists(self.merge_dir + "/concoct/" ) : os.makedirs(self.merge_dir + "/concoct/")
        script = raw_concoct_slurm % (self.merge_dir + "/concoct/" , self.name ,self.merge_dir + "/concoct/", self.merge_dir + "/concoct/")
        script += raw_concoct % (stats_out + self.name +"_contig_coverages_for_concoct.csv",self.merge_dir + "/454AllContigs.fna", 1000, self.merge_dir  + "/concoct/")
        with open(self.merge_dir + "concoct/" + "concoct_slurm_script.sh","w") as handle:
            handle.writelines(script)
    #        sh.sbatch(merge_dir  + "/concoct/" + "concoct_slurm_script.sh")

    def compile_contig_matrices(self):
        fs = sh.find(self.merge_dir, "-name", "*-smds.coverage.percontig").stdout.split("\n")[:-1]

        df = DataFrame.from_csv(fs[0],sep="\t")
        glob = df.ix[:,0:2]
        for f in fs:
            id = [c.replace("map_","") for c in f.split("/") if "map_" in c][0]
            values =  DataFrame.from_csv(f,sep="\t")["cov_mean_sample_0"]
            assert sum([a!=b for a,b in zip(values.index, glob.index)]) == 0
            if sum([a!=b for a,b in zip(values.index, glob.index)]) == 0:
                glob[id] = values
            else:
                print f, "is weird"

        samples = list(set([c.split("_")[2] for c in mat.columns if "EUSM" in c]))
        processed = {}
        processed['length'] = glob['length']
        processed['GC'] = glob['GC']
        for s in tqdm(samples):
            processed[s]=  glob[[c for c in glob.columns if s in c]].apply(sum, axis=1)
        glob = DataFrame.from_dict(processed)
        glob.to_csv(stats_out  + self.name + "_contig_coverages.csv",sep="\t")
        glob[samples].to_csv(stats_out  + self.name + "_contig_coverages_for_concoct.csv",sep="\t")
        return glob

    def hmmer_the_hell(self):
#        fs = sh.find(self.merge_dir + "bins/annotations/", "-name", "*.faa").stdout.split("\n")[:-1]
        script = raw_hmm_pipe  % (self.merge_dir, self.merge_dir,self.merge_dir, self.merge_dir, self.merge_dir, self.merge_dir)
        with open(self.merge_dir + "hmmer_slurm_script.sh","w") as handle:
            handle.writelines(script)
    #        sh.sbatch(merge_dir  + "/concoct/" + "concoct_slurm_script.sh")

    def unmapped_pulling(self):
        if os.path.exists(self.out_dir + "merged/454AllContigs.fna" ):
            if not os.path.exists(self.merge_dir + "unmapped_reads/"):
                os.makedirs(self.merge_dir + "unmapped_reads/")

            slurm_header = raw_map_slurm % ( self.out_dir, "deepbio_1", "unmap", self.merge_dir + "unmapped_reads/" , self.merge_dir + "unmapped_reads/")
            mapper_cmd = raw_unmap % ( self.merge_dir + "454AllContigs.fna", self.pair[0], self.pair[1], self.merge_dir + "unmapped_reads/" )
            self.map_script = slurm_header + mapper_cmd
        
            with open(self.merge_dir + "unmapped_reads/" + "unmap_slurm_script" + ".sh","w") as handle:
                handle.writelines(self.map_script)
#            sh.sbatch(self.merge_dir + "unmapped_reads/" + "unmap_slurm_script" + ".sh")
            print "Launched unmapper for", "deepbio_1"

