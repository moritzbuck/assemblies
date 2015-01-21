from __future__ import division

import sh
import os
import json
from tqdm import tqdm
from pandas import DataFrame

#from plumbing.autopaths import AutoPaths

raw_slurm = """#!/bin/bash
#SBATCH -D %s
#SBATCH -J ray_%s_%s
#SBATCH -o %sass.out
#SBATCH -e %sass.err
#SBATCH -A b2011035
#SBATCH -t 7-00:00:00
#SBATCH -n 16
#SBATCH -p node
#SBATCH -C mem512GB
#SBATCH --mail-user murumbii@gmail.com
#SBATCH --mail-type=ALL
"""

raw_map_slurm = """#!/bin/bash
#SBATCH -D %s
#SBATCH -J map_%s_%s
#SBATCH -o %smap.out
#SBATCH -e %smap.err
#SBATCH -A b2011035
#SBATCH -t 08:00:00
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
#SBATCH -A b2011035
#SBATCH -t 3-00:00:00
#SBATCH -n 16
#SBATCH -p node
#SBATCH --mail-user murumbii@gmail.com
#SBATCH --mail-type=ALL

"""

#### remove the paired reads make it to single reads
raw_map = """
/home/moritz/repos/metassemble/scripts/map/map-bowtie2-markduplicates_singles.sh -ct 16 %s %s %s %s %s
"""


raw_ray = """
k=`cat %scurrent_k.txt`
mpiexec -n 16 Ray231 -k $k -o %s_$k  %s
k=$((k + 10))
echo $k > %scurrent_k.txt
sbatch %s
"""

raw_newbler = """
/home/moritz/repos/metassemble/scripts/assembly/merge-asm-newbler.sh %s %s
"""

raw_path =  "/pica/v3/b2011138/TroutBogHypo/"
nr_threads = 16
executable = 'Ray231'
out_root = raw_path
stats_out = out_root + "stats/"
raw_ass = out_root + "raw_assemblies/"
merged_ass = out_root + "merged_assemblies/"
newbler = sh.Command("/home/moritz/repos/metassemble/scripts/assembly/merge-asm-newbler.sh")
minimus = sh.Command("/home/moritz/repos/metassemble/scripts/assembly/merge-asm-minimus2.sh")

class Assembly(object):

    def __repr__(self): return '<%s object %s with k-mer length %i>' % \
        (self.__class__.__name__, self.name, max(self.success_ks))
 
    
    def __init__(self, name, k=31, name_coass = ""):
        if isinstance(name, list):
            self.name =  "coass"+name_coass
            self.reads = [raw_path + n + ".fastq.gz" for n in name]
            self.coass = True
        else :
            self.coass = False
            self.name = name
            self.reads = raw_path + name + ".fastq.gz"
        self.k = k
        if not os.path.exists(raw_ass + self.name):
            os.makedirs(raw_ass +self.name)
        self.success_ks = [int(f.split("_")[1]) for f in os.listdir(raw_ass + self.name) if "ray" in f]
        self.out_dir = raw_ass + self.name + "/"
        self.merge_dir = raw_ass + self.name + "/merged/" 
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        if os.path.exists(self.out_dir + self.name + "_stats.json"):
            print "Loading stats for", self.name
            with open(self.out_dir + self.name + "_stats.json") as handle:
                self.stats = json.load(handle)
         

    def run(self):
        self.script = raw_slurm % ( self.out_dir, self.name, str(self.k), self.out_dir, self.out_dir)
        if not self.coass:
            cleans = "sickle se -f " + self.reads + " -t illumina -n -o " + self.reads.replace("fastq.gz","clean.fastq") +"\n" if not os.path.exists(self.reads.replace(".fastq.gz",".clean.fastq")) else ""
            self.arg = self.reads.replace("fastq.gz","clean.fastq")
            self.script = raw_slurm % ( self.out_dir, self.name, str(self.k), self.out_dir, self.out_dir)
            for c in cleans:
                self.script = self.script + c
        else :
            self.arg = " ".join([" -s " + r.replace("fastq.gz","clean.fastq") for r in self.reads])

        self.script = self.script + raw_ray % (self.out_dir, self.out_dir + "ray", self.arg, self.out_dir, self.out_dir + "slurm_script.sh")
        with open(self.out_dir + "current_k.txt","w") as handle:
            handle.writelines([str(self.k)])                    
        with open(self.out_dir + "slurm_script.sh","w") as handle:
            handle.writelines(self.script)
        sh.sbatch(self.out_dir + "slurm_script.sh")

    def make_stats(self, rerun = False):
        self.stats = {}
        if True: #os.path.exists(self.out_dir + "mapper/"):
            if os.path.exists(self.out_dir + self.name + "_stats.json") and not rerun:
                print "Loading stats for", self.name
                with open(self.out_dir + self.name + "_stats.json") as handle:
                    self.stats = json.load(handle)
            else:
                print "Computing stats for", self.name
                self.stats['raw_reads'] = self.reads
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
                        self.stats[k]['mapped_frac'] = float(sh.grep("overall",  self.out_dir + "mapper_" + str(k) + "/map.err").split(" ")[0][0:-1])/100
                        self.stats[k]['dupli_frac'] = float(sh.grep("Unknown", self.out_dir + "mapper/ray_" +str(k)+ "_" + self.name+ "-smd.metrics").split("\t")[7])
            
        with open(self.out_dir + self.name + "_stats.json", 'w') as handle:
            json.dump(self.stats, handle)
        
    def map_all(self):
        for k in self.success_ks:
            self.map(k)
                    
    def map_merged(self):
        if os.path.exists(self.merge_dir + "454AllContigs.fna"):
            slurm_header = raw_map_slurm % ( self.out_dir, self.name, "merged", self.merge_dir, self.merge_dir)
            mapper_cmd = raw_map % (self.reads, self.name, self.merge_dir + "454AllContigs.fna" , "merged", self.merge_dir + "mapper/")
            self.map_script = slurm_header + mapper_cmd
        
            with open(self.merge_dir +  "map_slurm_script.sh","w") as handle:
                handle.writelines(self.map_script)
            sh.sbatch(self.merge_dir + "map_slurm_script.sh")
            print "Launched mapper for", self.name, "with merged assembly"

    def map(self,k):
        if os.path.exists(self.out_dir + "ray_" + str(k) + "/Contigs.fasta"):
            if not os.path.exists(self.out_dir + "mapper_" + str(k) + "/"):
                os.makedirs(self.out_dir + "mapper_" + str(k) + "/")

            slurm_header = raw_map_slurm % ( self.out_dir, self.name, str(k), self.out_dir + "mapper_" + str(k) + "/", self.out_dir + "mapper_" + str(k) + "/")
            mapper_cmd = raw_map % (self.reads, self.name, self.out_dir + "ray_" + str(k) + "/Contigs.fasta", "ray_" + str(self.k), self.out_dir + "mapper_" + str(k) + "/")
            self.map_script = slurm_header + mapper_cmd
        
            with open(self.out_dir + "mapper_" + str(k) + "/" + "map_slurm_script_" + str(k) + ".sh","w") as handle:
                handle.writelines(self.map_script)
            sh.sbatch(self.out_dir + "mapper_" + str(k) + "/" + "map_slurm_script_" + str(k) + ".sh")
            print "Launched mapper for", self.name, "with k-mer size", k

    def map2ref(self,ref):
        if not os.path.exists( self.out_dir + "mapper_" + "2ref" + "/"):
            os.makedirs( self.out_dir + "mapper_" + "2ref" + "/")
            
        slurm_header = raw_map_slurm % ( self.out_dir, self.name, "2ref", self.out_dir + "mapper_" + "2ref" + "/", self.out_dir + "mapper_" + "2ref" + "/")
        mapper_cmd = raw_map % (self.reads, self.name, ref, "ray_" + "2ref", self.out_dir + "mapper_" + "2ref" + "/")
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
        header = raw_newbler_slurm % ( self.merge_dir, self.name, self.merge_dir, self.merge_dir)
        with open(self.merge_dir + "merge_slurm_script.sh","w") as handle:
            handle.writelines(header + newbler)
        sh.sbatch(self.merge_dir + "merge_slurm_script.sh")

def run_all():
    asses = [Assembly(p.replace(".fastq.gz","")) for p in os.listdir(raw_path) if ".fastq.gz" in p and not "clean" in p]
    for a in asses:
        a.run()

def merge_all():
    asses = [Assembly(p.replace(".fastq.gz","")) for p in os.listdir(raw_path) if ".fastq.gz" in p and not "clean" in p]
    for a in asses:
        a.merge()

                
def map_all_singles():
    asses = [Assembly(p.replace(".fastq.gz","")) for p in os.listdir(raw_path) if ".fastq.gz" in p and not "clean" in p]
    for a in asses:
        a.map_all()

def map_all_merged():
    asses = [Assembly(p.replace(".fastq.gz","")) for p in os.listdir(raw_path) if ".fastq.gz" in p and not "clean" in p]
    for a in asses:
        a.map_merged()

def merge_final():
    asses = [Assembly(p.replace(".fastq.gz","")) for p in os.listdir(raw_path) if ".fastq.gz" in p and not "clean" in p]
    fastas = [ass.merge_dir + "454AllContigs.fna" for ass in asses]
    if not os.path.exists(merged_ass) : os.makedirs(merged_ass)
    in_size = sum([os.stat(f).st_size for f in fastas])
    newbler = raw_newbler % (merged_ass, " ".join(fastas))
    header = raw_newbler_slurm % ( merged_ass,"merged_all", merged_ass,merged_ass)
    with open(merged_ass + "merge_slurm_script.sh","w") as handle:
        handle.writelines(header + newbler)
    sh.sbatch(merged_ass + "merge_slurm_script.sh")

def map_final():
    asses = [Assembly(p.replace(".fastq.gz","")) for p in os.listdir(raw_path) if ".fastq.gz" in p and not "clean" in p]
    for a in asses:
        a.map2ref(merged_ass +"454AllContigs.fna")

def compile_contig_matrices():
    fs = sh.find(raw_ass, "-name", "*-smds.coverage.percontig").stdout.split("\n")[:-1]
    fs = [f for f in fs if "_2ref" in f]
    df = DataFrame.from_csv(fs[0],sep="\t")
    glob = df.ix[:,0:2]
    for f in fs:
        id = [c for c in f.split("/") if "IH" in c][0]
        values =  DataFrame.from_csv(f,sep="\t")["cov_mean_sample_0"]
        assert sum([a!=b for a,b in zip(values.index, glob.index)]) == 0
        if sum([a!=b for a,b in zip(values.index, glob.index)]) == 0:
            glob[id] = values
        else:
            print f, "is weird"
    glob.to_csv(stats_out  + "all_contig_coverages.csv",sep="\t")
    

def write_lib_stats(assemblies):
    header = "\traw number of read-pairs\tqc-ed number of read-pairs\tpercentage of clean reads"
    lines = [[header]]+[a.format_lib_stats() for a in assemblies]
    lines = sum(lines, [])
    lines = [l+"\n" for l in lines]
    with open(stats_out + "library_stats.csv","w") as handle:
        handle.writelines(lines)
    return lines

def write_assembly_stats(assemblies):
    header = "sample\tk-mer size\tassembly size\tnumber of contigs\tmax contig length\tn50 length\tproportion of duplicated reads\tproportion of successfuly mapped reads\tassembly path"
    lines = [header]+[a.format_assem_stats() for a in assemblies]
    lines = [l+"\n" for l in lines]
    with open(stats_out + "assemblies_stats.csv","w") as handle:
        handle.writelines(lines)
    return lines
    
def write_merged_stats(samples):
    header = "assembly size\tnumber of contigs\tmax contig length\tn50 length\tassembly path"
    lines = [header]
    ff = lambda x :  str(round(float(x)/1000000,2))+"M"
    pp = lambda x :  str(int(x*100))+"%"
    kk = lambda x :  str(round(float(x)/1000,2))+"k"

    for s in samples:
        t=sh.assemstats(0, merged_ass + s + "/454AllContigs.fna")
        stats = {a:b for a,b in  zip(*[[w.strip() for w in  l.split("\t")][1:] for l in str(t).split("\n")[0:2]])}
        lines.append("\t".join([s,ff(stats['sum']), ff(stats['n']), kk(stats['max']),kk(stats['n50_len']), merged_ass + s + "/454AllContigs.fna"]))
    lines = [l+"\n" for l in lines]

    with open(stats_out + "merged_assemblies_stats.csv","w") as handle:
        handle.writelines(lines)
    return lines

all_2007_asses = ["IHXI","IHWF","IHUS","IHWI","IHWN","IHWG","IHWH","IHUZ","IHSA","IHWA","IHWY","IHWB","IHWU","IHXG","IHXO","IHXN","IHXW"]
all_2008_asses = ["IHPY","IHSB","IHWS","IHWX","IHWW","IHXP","IHXT","IHXH","IHXA","IHXF","IHXU","IHXS","IHXX","IHPO","IHWP","IHPN"]
all_2009_asses = ["IHWT","IHTZ","IHXY","IHUN","IHUO","IHUP","IHUA","IHUB","IHUC","IHUH","IHUI","IHUF"]
#all_31_asses =   [Assembly(sam, samples[sam],k=31)  for sam in samples]
#all_41_asses =   [Assembly(sam, samples[sam],k=41)  for sam in samples]
#all_51_asses =   [Assembly(sam, samples[sam],k=51)  for sam in samples]
            
#frac_0p1 = [ass for ass in all_31_asses + all_41_asses + all_51_asses if ass.name in [k for k,v in samples.iteritems() if v['fraction'] == 0.1]]
#frac_0p8 = [ass for ass in all_31_asses + all_41_asses + all_51_asses if ass.name in [k for k,v in samples.iteritems() if v['fraction'] == 0.8]]
#frac_3p0 = [ass for ass in all_31_asses + all_41_asses + all_51_asses if ass.name in [k for k,v in samples.iteritems() if v['fraction'] == 3.0]]


#all_asses = all_31_asses + all_41_asses + all_51_asses
