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
#SBATCH -A b2013151
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
#SBATCH -A b2011032
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
#SBATCH -A b2013151
#SBATCH -t 3-00:00:00
#SBATCH -n 16
#SBATCH -p node
#SBATCH --mail-user murumbii@gmail.com
#SBATCH --mail-type=ALL

"""


raw_map = """
/home/moritz/repos/metassemble/scripts/map/map-bowtie2-markduplicates.sh -ct 16 %s %s %s %s %s %s
"""


raw_ray = """
mpiexec -n 16 Ray231 -amos -enable-neighbourhoods -k %s -o %s %s
"""

raw_newbler = """
/home/moritz/repos/metassemble/scripts/assembly/merge-asm-newbler.sh %s %s
"""

raw_path =  "/pica/v8/b2013151/INBOX/"
rub_path = raw_path + "B.Bergman_14_02/"
ill_path = raw_path + "B.Bergman_14_01/"
samples_json = "/home/moritz/people/mimebs/samples.json"
nr_threads = 16
executable = 'Ray231'
out_root = '/pica/v8/b2013151/private/Bergman_14/'
stats_out = out_root + "stats/"
raw_ass = out_root + "raw_assemblies/"
merged_ass = out_root + "merged_assemblies/"
newbler = sh.Command("/home/moritz/repos/metassemble/scripts/assembly/merge-asm-newbler.sh")
minimus = sh.Command("/home/moritz/repos/metassemble/scripts/assembly/merge-asm-minimus2.sh")


with open(samples_json) as handle:
    samples = json.load(handle)

class Assembly(object):

    def __repr__(self): return '<%s object %s with k-mer length %i>' % \
        (self.__class__.__name__, self.name, self.k)
 
    
    def __init__(self, name, sample_dic, k=31):
        self.dirs = [ rub_path + sample_dic['rubicon_path'] + "/", ill_path + sample_dic['illumina_path'] + "/"]
        self.name = name
        self.subdirs = sum([[d+t+"/"  for t in os.listdir(d)]  for d in self.dirs],[])
        self.pairs = sum([zip([t+f for f in os.listdir(t) if ".gz" in f and "md5" not in f and "_2." not in f],[t+f for f in os.listdir(t) if ".gz" in f and  "md5" not in f and "_1." not in f])  for t in self.subdirs],[])
        self.k = k
        self.out_dir = raw_ass + self.name + "/" + str(self.k) + "/"
        self.frac = sample_dic['fraction']
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        if os.path.exists(self.out_dir + self.name + "_stats.json"):
            print "Loading stats for", self.name
            with open(self.out_dir + self.name + "_stats.json") as handle:
                self.stats = json.load(handle)
         

    def run(self):
        if not os.path.exists( self.out_dir + "ray" ) and os.path.exists(raw_ass + self.name + "/41/ray/Contigs.fasta"):
            cleans = ["sickle pe -f " + p[0] + " -r " + p[1] + " -t sanger -n -o " + p[0].replace("fastq.gz","clean.fastq") + " -p " + p[1].replace("fastq.gz","clean.fastq") + " -s " +  p[0].replace("_1.fastq.gz","_single.clean.fastq") +"\n" for p in self.pairs if not os.path.exists(p[0].replace("_1.fastq.gz","_single.clean.fastq"))]     
            self.arg = sum([sum([["-p"],list(p)],[]) for p in self.pairs],[])
            self.arg = [p.replace("fastq.gz","clean.fastq") for p in self.arg]
            self.script = raw_slurm % ( self.out_dir, self.name, str(self.k), self.out_dir, self.out_dir)
            for c in cleans:
                self.script = self.script + c 
            self.script = self.script + raw_ray % (str(self.k), self.out_dir + "ray", " ".join(self.arg))
        
            with open(self.out_dir + "slurm_script.sh","w") as handle:
                handle.writelines(self.script)
            sh.sbatch(self.out_dir + "slurm_script.sh")

    def make_stats(self, rerun = False):
        self.stats = {}

        if os.path.exists(self.out_dir + "mapper/"):
            if os.path.exists(self.out_dir + self.name + "_stats.json") and not rerun:
                print "Loading stats for", self.name
                with open(self.out_dir + self.name + "_stats.json") as handle:
                    self.stats = json.load(handle)
            else:
                print "Computing stats for", self.name
                self.stats['k'] = self.k
                self.stats['raw_reads'] = self.pairs
                self.stats['output'] = self.out_dir
                print "Computing raw number of paired-reads"
                self.stats['pairs_count'] = {"_".join(p[0].split("/")[-1].split("_")[0:-1]) : int(sh.zgrep("-Ec", "$", p[0]))/4 for p in tqdm(self.pairs)}
                print "Computing clean number of paired-reads"
                self.stats['clean_pairs_count'] = {"_".join(p[0].split("/")[-1].split("_")[0:-1]) : int(sh.grep("-Ec", "$", p[0].replace("fastq.gz","clean.fastq")))/4 for p in tqdm(self.pairs)}
                self.stats['ratios'] = {k : self.stats['clean_pairs_count'][k]/tot for k,tot in self.stats['pairs_count'].iteritems()}
                self.stats['map_stats'] = {}
                self.stats['map_stats']['mapped_frac'] = float(sh.grep("overall",  self.out_dir + "map.err").split(" ")[0][0:-1])/100
                self.stats['map_stats']['dupli_frac'] = float(sh.grep("Unknown", self.out_dir + "mapper/ray_" +str(self.k)+ "_" + self.name+ "-smd.metrics").split("\t")[7])
                self.stats['Success'] =  os.path.exists(self.out_dir + "ray/Contigs.fasta")
                if self.stats['Success']:
                    t=sh.assemstats(0,self.out_dir + "ray/Contigs.fasta" )
                    self.stats['Assembly'] = {a:b for a,b in  zip(*[[w.strip() for w in  l.split("\t")][1:] for l in str(t).split("\n")[0:2]])}
        else:
            self.stats['Success'] =  os.path.exists(self.out_dir + "ray/Contigs.fasta")
            
        with open(self.out_dir + self.name + "_stats.json", 'w') as handle:
            json.dump(self.stats, handle)
        
        
    def map(self):
        if os.path.exists(self.out_dir + "ray/Contigs.fasta") and not os.path.exists(self.out_dir + "mapper/"):
            lefties = ",".join([p[0] for p in self.pairs])
            righties = ",".join([p[1] for p in self.pairs])

            slurm_header = raw_map_slurm % ( self.out_dir, self.name, str(self.k), self.out_dir, self.out_dir)
            mapper_cmd = raw_map % (lefties, righties, self.name, self.out_dir + "ray/Contigs.fasta", "ray_" + str(self.k), self.out_dir + "mapper/")
            self.map_script = slurm_header + mapper_cmd
        
            with open(self.out_dir + "map_slurm_script.sh","w") as handle:
                handle.writelines(self.map_script)
            sh.sbatch(self.out_dir + "map_slurm_script.sh")
            print "Launched mapper for", self.name, "with k-mer size", self.k

    def format_lib_stats(self):
        libs = self.stats['clean_pairs_count'].keys()
        lib_type = {l: "illumina" if "P1203" in l else "rubicon" for l in libs }
        libs = sorted(libs, key=lambda l : lib_type[l])
        ff = lambda x :  str(round(x/1000000,2))+"M"
        pp = lambda x :  str(round(x*100,2))+"%"
        lines = ["\t".join([l, ff(self.stats['pairs_count'][l]), ff(self.stats['clean_pairs_count'][l]), pp(self.stats['ratios'][l]), lib_type[l]]) for l in libs]
        lines = [(self.name if i==0 else "") + "\t" + l  for i,l in enumerate(lines)]
        return lines

    def format_assem_stats(self):
        ff = lambda x :  str(round(float(x)/1000000,2))+"M"
        pp = lambda x :  str(round(x*100,2))+"%"
        kk = lambda x :  str(round(float(x)/1000,2))+"k"
        if self.stats['Success']:
            line = "\t".join([self.name,str(self.frac), str(self.k),ff(self.stats['Assembly']['sum']), ff(self.stats['Assembly']['n']), kk(self.stats['Assembly']['max']),kk(self.stats['Assembly']['n50_len']), pp(self.stats['map_stats']['dupli_frac']), pp(self.stats['map_stats']['mapped_frac']), self.stats['output'] + "ray/Contigs.fasta" ])
        else:
            line = "\t".join([self.name,str(self.frac), str(self.k) ,"NA" ,"NA","NA","NA","NA","NA","NA"])
        return line

    
def merge(assemblies, name):
    fastas = [ass.out_dir + "ray/Contigs.fasta" for ass in assemblies if os.path.exists(ass.out_dir + "ray/Contigs.fasta")]
    if not os.path.exists(merged_ass + name) : os.makedirs(merged_ass + name)
    in_size = sum([os.stat(f).st_size for f in fastas])
    newbler = raw_newbler % (merged_ass + name, " ".join(fastas))
    header = raw_newbler_slurm % ( merged_ass + name +"/", name, merged_ass + name+"/", merged_ass + name+"/")
    with open(merged_ass + name + "/" + "map_slurm_script.sh","w") as handle:
        handle.writelines(header + newbler)

def map_all(assemblies, name):
    for a in assemblies:
        lefties = ",".join([p[0] for p in a.pairs])
        righties = ",".join([p[1] for p in a.pairs])
        out_dir = merged_ass + name + "/mappers/" + a.name + "/"
        if not os.path.exists(out_dir) : os.makedirs(out_dir)
        slurm_header = raw_map_slurm % ( out_dir, name, a.name, out_dir, out_dir)
        mapper_cmd = raw_map % (lefties, righties, a.name,  merged_ass + name + "/454AllContigs.fna" , name, out_dir)
        map_script = slurm_header + mapper_cmd
        
        with open(out_dir + "map_slurm_script.sh","w") as handle:
            handle.writelines(map_script)
            
        sh.sbatch(out_dir + "map_slurm_script.sh")
        print "Launched mapper for", name, "with sample", a.name

def compile_contig_matrices():
    for m in ["frac_0p1", "frac_0p8", "frac_3p0"]:
        fs = sh.find(merged_ass + m, "-name", "*-smds.coverage.percontig").stdout.split("\n")[:-1]
        df = DataFrame.from_csv(fs[0],sep="\t")
        glob = df.ix[:,0:2]
        for f in fs:
            id = [c for c in f.split("/") if "GS8" in c][0]
            values =  DataFrame.from_csv(f,sep="\t")["cov_mean_sample_0"]
            assert sum([a!=b for a,b in zip(values.index, glob.index)]) == 0
            glob[id] = values
        glob.to_csv(stats_out + m + "_contig_coverages.csv",sep="\t")
    
def merge_all():
    merge(frac_0p8, "frac_0p8")
    merge(frac_0p1, "frac_0p1")
    merge(frac_3p0, "frac_3p0")

def write_lib_stats(assemblies):
    header = "\tlibrary\traw number of read-pairs\tqc-ed number of read-pairs\tpercentage of clean reads\tlibrary type"
    lines = [[header]]+[a.format_lib_stats() for a in assemblies]
    lines = sum(lines, [])
    lines = [l+"\n" for l in lines]
    with open(stats_out + "library_stats.csv","w") as handle:
        handle.writelines(lines)
    return lines

def write_assembly_stats(assemblies):
    header = "sample\tfraction\tk-mer size\tassembly size\tnumber of contigs\tmax contig length\tn50 length\tproportion of duplicated reads\tproportion of successfuly mapped reads\tassembly path"
    lines = [header]+[a.format_assem_stats() for a in assemblies]
    lines = [l+"\n" for l in lines]
    with open(stats_out + "assemblies_stats.csv","w") as handle:
        handle.writelines(lines)
    return lines
    
def write_merged_stats():
    header = "fraction\tassembly size\tnumber of contigs\tmax contig length\tn50 length\tassembly path"
    sets = [(frac_0p1, "frac_0p1"),(frac_0p8, "frac_0p8"),(frac_3p0, "frac_3p0")]
    lines = [header]
    ff = lambda x :  str(round(float(x)/1000000,2))+"M"
    pp = lambda x :  str(int(x*100))+"%"
    kk = lambda x :  str(round(float(x)/1000,2))+"k"

    for s in sets:
        t=sh.assemstats(0, merged_ass + s[1] + "/454AllContigs.fna")
        stats = {a:b for a,b in  zip(*[[w.strip() for w in  l.split("\t")][1:] for l in str(t).split("\n")[0:2]])}
        lines.append("\t".join([s[1],ff(stats['sum']), ff(stats['n']), kk(stats['max']),kk(stats['n50_len']), merged_ass + s[1] + "/454AllContigs.fna"]))
    lines = [l+"\n" for l in lines]

    with open(stats_out + "merged_assemblies_stats.csv","w") as handle:
        handle.writelines(lines)
    return lines

        
    
all_31_asses =   [Assembly(sam, samples[sam],k=31)  for sam in samples]
all_41_asses =   [Assembly(sam, samples[sam],k=41)  for sam in samples]
all_51_asses =   [Assembly(sam, samples[sam],k=51)  for sam in samples]
            
frac_0p1 = [ass for ass in all_31_asses + all_41_asses + all_51_asses if ass.name in [k for k,v in samples.iteritems() if v['fraction'] == 0.1]]
frac_0p8 = [ass for ass in all_31_asses + all_41_asses + all_51_asses if ass.name in [k for k,v in samples.iteritems() if v['fraction'] == 0.8]]
frac_3p0 = [ass for ass in all_31_asses + all_41_asses + all_51_asses if ass.name in [k for k,v in samples.iteritems() if v['fraction'] == 3.0]]


all_asses = all_31_asses + all_41_asses + all_51_asses
