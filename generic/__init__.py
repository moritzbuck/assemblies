from __future__ import division

import sh
import os
import json
from  generic.scripts import *
from tqdm import tqdm
from pandas import DataFrame
import random


class Sample(object):
    def __repr__(self):
        return '<%s object %s with %i libraries %i of which are paired>' % (self.__class__.__name__, self.name, len(self.reads), sum([len(r) == 2 for r in self.reads]) )        
        
    def __init__(self, name, reads, path, metadata = None, no_cleans = False):
        self.no_cleans = no_cleans
        self.name = name
        self.reads = reads
        
        self.path = path + self.name + "/"
        self.metadata = metadata
        self.clean_path = self.path + "clean/"
        self.normalize_path = self.path + "normalized/"
        self.normalized_sample = self.normalize_path + self.name + "_unpaired_normalised.fastq"
        if not os.path.exists(self.clean_path):
            os.makedirs(self.clean_path)


    def clean(self, project = "b2011138", submit = False):
        script = make_slurm_header(self.path, "clean_" + self.name, proj = project, cores = 1)

        libraries = self.get_libraries()
        
        script += make_parallel_sickle_script(libraries, threads = 1)

        with open(self.path +  "clean_" + self.name + "_script.sh","w") as handle:
            handle.writelines(script)
        if submit :
            clean_job = sh.sbatch(self.path +  "clean_" + self.name + "_script.sh")
            print "Launched clean_" + self.name

    def normalize(self, project, submit = False, shuffle = True):
        cores = 4
        script = make_slurm_header(self.normalize_path, "khmer_norm_" + self.name, proj = project, cores = cores, time = "7-00:00:00")

        libraries = self.get_all_reads()
        if shuffle:
            random.shuffle(libraries)
        
        script += make_khmer_script(libraries, self.normalized_sample, threads = cores)

        with open(self.path +  "normed_" + self.name + "_script.sh","w") as handle:
            handle.writelines(script)
        if submit :
            clean_job = sh.sbatch(self.path +  "normed_" + self.name + "_script.sh")
            print "Launched clean_" + self.name
            
    def get_libraries(self):
        pairs = [p for p in self.reads if p.has_key('1')]
        interlaced = [p for p in self.reads if p.has_key('I')]
        singles = [p for p in self.reads if p.has_key('U')]
        all_libs = []
        all_libs += [{'1' : p['1'], '2' : p['2'], 'c1' : self.clean_path + "clean_fwd_" + os.path.basename(p['1']).replace(".gz",""), 'c2' : self.clean_path + "clean_rev_" + os.path.basename(p['2']).replace(".gz",""), 'cU' : self.clean_path + "clean_unp_" + os.path.basename(p['1']).replace(".gz","") } for p in pairs]
        all_libs += [{'I' : p['I'], 'c1' : self.clean_path + "clean_fwd_" + os.path.basename(p['I']).replace(".gz",""), 'c2' : self.clean_path + "clean_rev_" + os.path.basename(p['I']).replace(".gz",""), 'cU' : self.clean_path + "clean_unp_" + os.path.basename(p['I']).replace(".gz","") } for p in interlaced]
        all_libs += [{'U' :  p['U'],  'cU' : self.clean_path + "clean_unp_" + os.path.basename(p['U']).replace(".gz","")} for p in singles]
        return all_libs
        

    def get_pairs(self, clean = True):
        if self.no_cleans:
            clean = False
        pairs = [p for p in self.reads if p.has_key('1')]
        interlaced = [p for p in self.reads if p.has_key('I')]
        if clean:
            fwd = [self.clean_path + "clean_fwd_" + os.path.basename(p['1']).replace(".gz","") for p in pairs]
            rev = [self.clean_path + "clean_rev_" + os.path.basename(p['2']).replace(".gz","") for p in pairs]
            fwd += [self.clean_path + "clean_fwd_" + os.path.basename(p['I']).replace(".gz","") for p in interlaced]
            rev += [self.clean_path + "clean_rev_" + os.path.basename(p['I']).replace(".gz","") for p in interlaced]
            return {'1' : fwd, '2' : rev}
        else :
            fwd = [p['1'] for p in pairs]
            rev = [p['2'] for p in pairs]
            inl = [p['I'] for p in interlaced]
            return {'1' : fwd, '2' : rev, 'I': inl}

    def get_singles(self, clean = True):
        if self.no_cleans:
            clean = False

        pairs = [p for p in self.reads if p.has_key('1')]
        inls = [p for p in self.reads if p.has_key('I')]
        singles = [p for p in self.reads if p.has_key('U')]
        if clean:
            single = [self.clean_path + "clean_unp_" + os.path.basename(p['U']).replace(".gz","") for p in singles]
            single_pairs = [self.clean_path + "clean_unp_" + os.path.basename(p['1']).replace(".gz","") for p in pairs]
            single_inls = [self.clean_path + "clean_unp_" + os.path.basename(p['I']).replace(".gz","") for p in inls]
            return { 'U' : single, '1': single_pairs, 'I' : single_inls}
        else :
            single = [p['U'] for p in singles]
            return {'U': single}

    def get_all_reads(self, clean = True):
        if self.no_cleans:
            clean = False
        pairs = [p for p in self.reads if p.has_key('1')]
        interlaced = [p for p in self.reads if p.has_key('I')]
        singles = [p for p in self.reads if p.has_key('U')]
        all_reads = []
        if clean:
            all_reads += [self.clean_path + "clean_fwd_" + os.path.basename(p['1']).replace(".gz","") for p in pairs]
            all_reads += [self.clean_path + "clean_rev_" + os.path.basename(p['2']).replace(".gz","") for p in pairs]
            all_reads += [self.clean_path + "clean_fwd_" + os.path.basename(p['I']).replace(".gz","") for p in interlaced]
            all_reads += [self.clean_path + "clean_rev_" + os.path.basename(p['I']).replace(".gz","") for p in interlaced]
            all_reads += [self.clean_path + "clean_unp_" + os.path.basename(p['U']).replace(".gz","") for p in singles]
            all_reads += [self.clean_path + "clean_unp_" + os.path.basename(p['1']).replace(".gz","") for p in pairs]
            all_reads += [self.clean_path + "clean_unp_" + os.path.basename(p['I']).replace(".gz","") for p in interlaced]
            return all_reads
        else :
            all_reads = sum([r.values() for r in self.reads],[])
            return all_reads

                

class Assembly(object):

    def __repr__(self): return '<%s object %s with %i samples>' % \
        (self.__class__.__name__, self.name, len(self.samples) )
 
    
    def __init__(self, data_json, project):

        self.project = project
        
        self.json_path = data_json

        with open(self.json_path) as handle:
            self.run_data = json.load(handle)
    
        self.name = self.run_data['name']
        self.out_root = self.run_data['out_path']
        self.stats_out = self.out_root + "stats/"
        self.raw_ass = self.out_root + "assemblies/"
        self.sample_path = self.out_root +"samples/"
                        
        if not os.path.exists(self.raw_ass):
            os.makedirs(self.raw_ass)
        if not os.path.exists(self.stats_out):
            os.makedirs(self.stats_out)
        if not os.path.exists(self.out_root):
            os.makedirs(self.out_root)
        if not os.path.exists(self.sample_path):
            os.makedirs(self.sample_path)

        self.job_file = self.out_root + "job_ids.json"
        if not os.path.exists(self.job_file):
            self.job_ids={}
        else:
            with open(self.job_file) as handle:
                self.job_ids = json.load(handle)

                        
        self.samples = []
        for k in self.run_data['samples']:
            vals = self.run_data['samples'][k].copy()
            reads = vals['reads']
            del vals['reads']
            self.samples.append(Sample(k,reads, self.sample_path, vals))

    def get_all_reads(self, clean = True, kind = "all"):
        if kind == "all":
            return sum([s.get_all_reads() for s in self.samples],[])
        if kind == "paired":
            temp =  zip(*[s.get_pairs( clean=clean) for s in self.samples])
        else :
            temp = zip(*[s.get_singles( clean=clean) for s in self.samples])
        return [sum(t,[]) for t in temp]

    def get_clean_sample_dict(self):
        sample_dict = {}
        for s in self.samples:
            paired = s.get_pairs()
            unpaired = s.get_singles()
 
            sample_dict[s.name] = {}
            sample_dict[s.name]['1'] = paired['1']
            sample_dict[s.name]['2'] = paired['2']
            sample_dict[s.name]['U'] = sum(unpaired.values(),[])
        return  sample_dict
    
    def clean_all(self, submit = False):
        script = make_slurm_header(self.sample_path, "clean_all", self.project)

        libraries = sum([s.get_libraries() for s in self.samples], [])
        
        script += make_parallel_sickle_script(libraries)

        with open(self.sample_path +  "clean_all_script.sh","w") as handle:
            handle.writelines(script)
        if submit :
            clean_job = sh.sbatch(self.sample_path +  "clean_all_script.sh")
            self.job_ids['clean_all'] = clean_job.split()[-1]
            print "Launched clean_all"
            with open(self.job_file, 'w') as handle:
                json.dump(self.job_ids, handle)

                

    def megahit_assembly(self, name="megahit", max_read_len = 291, submit = False, reads = None, deps = None):
        if not deps:
            if self.job_ids.has_key('clean_all'):
                deps = [self.job_ids['clean_all']]
            else:
                deps = ['1']
        path = self.raw_ass + name +"/"
        if not reads:
            reads = self.get_all_reads()

        script = make_slurm_header(path, "assembly_" + name, self.project, constraint = "fat", time = "7-00:00:00", deps = deps)
        script += make_megahit_script(reads, path)
        with open(path +  "megahit_assembly_script.sh","w") as handle:
            handle.writelines(script)
        if submit :
            job = sh.sbatch(path + "megahit_assembly_script.sh")
            if not self.job_ids.has_key('assembly'):
                self.job_ids['assembly'] = {}
            self.job_ids['assembly'][name] = job.split()[-1]
            print "Launched megahit assembly named", name
            with open(self.job_file, 'w') as handle:
                json.dump(self.job_ids, handle)

    def minimus_scaffolder(self, name="megahit", submit = False):
        if self.job_ids.has_key('assembly'):
            deps = [self.job_ids['assembly'][name]]
        else:
            deps = ['1']
        
        path = self.raw_ass + name +"/"
        script = make_slurm_header(path, "scaffminimus_" + name, self.project, constraint = "medium", time = "4-00:00:00", deps = deps)
        script += make_minimus_script(path + "contigs.fasta", path + "scaffolds.fasta")
        with open(path +  "scaffminimus_script.sh","w") as handle:
            handle.writelines(script)
        if submit :
            job = sh.sbatch(path + "scaffminimus_script.sh")
            if not self.job_ids.has_key('scaff'):
                self.job_ids['scaff'] = {}
            self.job_ids['scaff'][name] = job.split()[-1]
            print "Launched megahit assembly named", name
            with open(self.job_file, 'w') as handle:
                json.dump(self.job_ids, handle)

    def map_samples(self, reference, bowtie_only = False, submit = False, deps = None):
        path = os.path.dirname(reference) + "/mapping_" + os.path.basename(reference) +"/"
        script = make_slurm_header(path, "mapping_samples_to_" + os.path.basename(reference), self.project,  time = "4-00:00:00", deps = deps)
        script += make_bbmapping_script(reference, path, self.get_clean_sample_dict(), only_bowtie = bowtie_only)
#        script += make_star_mapping_script(reference, path, self.get_clean_sample_dict())
#        script += make_mapping_script(reference, path, self.get_clean_sample_dict(), only_bowtie = bowtie_only)
        with open(path +  "mapping_script.sh","w") as handle:
            handle.writelines(script)
        if submit :
            job = sh.sbatch(path + "mapping_script.sh")
            if not self.job_ids.has_key('maps'):
                self.job_ids['maps'] = {}
            self.job_ids['maps'][reference] = job.split()[-1]
            print "Launched mapper on", reference
            with open(self.job_file, 'w') as handle:
                json.dump(self.job_ids, handle)

    def besst_scaffold(self, reference, bowtie_only = False, submit = False, deps = None):
        maps_path = os.path.dirname(reference) + "/mapping_" + os.path.basename(reference) +"/"
        path = os.path.dirname(reference) + "/BESST_scaffolding/"
        script = make_slurm_header(path, "BESST_scaff_of_" + os.path.basename(reference), self.project,  time = "4-00:00:00", deps = deps, cores=1)
        script += make_besst_script(reference, path, self.get_clean_sample_dict(), maps_path)
        with open(path +  "besst_script.sh","w") as handle:
            handle.writelines(script)
        if submit :
            job = sh.sbatch(path + "besst_script.sh")
            if not self.job_ids.has_key('scaffolds'):
                self.job_ids['scaffolds'] = {}
            self.job_ids['scaffolds'][reference] = job.split()[-1]
            print "Launched mapper on", reference
            with open(self.job_file, 'w') as handle:
                json.dump(self.job_ids, handle)

    def concocting(self, reference, coverage_csv, min_len = 1000, kmer=4, submit = False, deps = None):
        covfile = coverage_csv
        path = os.path.dirname(reference) + "/binning/"
        script = make_slurm_header(path, "concocting_of_" + os.path.basename(reference), self.project,  time = "4-00:00:00", deps = deps)
        script += make_concoct_script(covfile, reference, path, min_len, kmer)
        with open(path +  "concoct_script.sh","w") as handle:
            handle.writelines(script)
        if submit :
            job = sh.sbatch(path + "concoct_script.sh")
            if not self.job_ids.has_key('bins'):
                self.job_ids['bins'] = {}
            self.job_ids['bins'][reference] = job.split()[-1]
            print "Launched binning on", reference
            with open(self.job_file, 'w') as handle:
                json.dump(self.job_ids, handle)

                                
    def map_binannot(self, infasta, name, submit = False, deps = None):
        path = os.path.dirname(infasta) + "/annotation/" 
        script = make_slurm_header(path, "annotate_bin_" + name, self.project,  time = "1-00:00:00", deps = deps)
        script += make_bin_bmfa(infasta, path, name = name)
        with open(path +  "annotate_bin.sh","w") as handle:
            handle.writelines(script)
        if submit :
            job = sh.sbatch(path + "annotate_bin.sh")
            if not self.job_ids.has_key('annotate'):
                self.job_ids['annotate'] = {}
            self.job_ids['maps'][name] = job.split()[-1]
            print "Launched mapper on", infasta
            with open(self.job_file, 'w') as handle:
                json.dump(self.job_ids, handle)

    def annotate_all_bins(self, path, submit = False, deps = None, min_size = None):
        if min_size:
            bins = [".".join(b.split(".")[:-1]) for b in os.listdir(path) if "bin" in b and b[-5:]=="fasta"  and  os.path.getsize(path + b) > min_size]
        else:
            bins = [".".join(b.split(".")[:-1]) for b in os.listdir(path) if "bin" in b and b[-5:]=="fasta"]

        for b in bins:
            nproc = 1 
            if os.path.getsize(path + b + ".fasta") > 5000000:
                nproc = 16
            script = make_slurm_header(path + b + "/", "annotate_" + b, self.project,  time = "12:00:00", deps = deps, cores = nproc)
            script += make_bin_bmfa(path + b + ".fasta", path + b, name = b, threads = nproc)
            with open(path + b + "/" +  "annotate_bin.sh","w") as handle:
                handle.writelines(script)
            if submit :
                job = sh.sbatch(path + b + "/" + "annotate_bin.sh")
                if not self.job_ids.has_key('annotate'):
                    self.job_ids['annotate'] = {}
                print "Launched annotater on", b
                with open(self.job_file, 'w') as handle:
                    json.dump(self.job_ids, handle)

    def reconcoct_all_good_bins(self, covfile, path, submit = False, deps = None, min_contamination = 10):
        bins = [".".join(b.split(".")[:-1]) for b in os.listdir(path) if "bin" in b and "fasta" in b]
        bins = [b for b in bins if b != "bin_non-bin"]
        contams ={}
        for b in bins: 
            with open(path + b + "/checkm.txt") as handle: 
                lines =  handle.readlines()
            if len(lines) > 3:
                contams[b] = float([l.split()[13] for l in lines if b in l][0])

        bins = [b for b,c in contams.iteritems() if c > min_contamination]
        bins = [b for b in bins if not os.path.exists(path + b + "/" + b + "_non-bin.fasta")]
        for b in bins:
            
            nproc = 16
            script = make_slurm_header(path + b + "/", "reconcoct_" + b , self.project,  time = "12:00:00", deps = deps, cores = nproc)
            script += make_reconcoct_script(covfile, path + b + ".fasta",  b + "_" , cutoff = 5000)

            with open(path + b + "/" +  "reconcoct_bin.sh","w") as handle:
                handle.writelines(script)
            if submit :
                job = sh.sbatch(path + b + "/" + "reconcoct_bin.sh")
                print "Launched re-concocter on", b

    def diginorm(self, deps = None, submit = False, threads=4):
        path = self.sample_path
        script = make_slurm_header(path, "khmer_the_hell_" + self.name, self.project,  time = "4-00:00:00", deps = deps, cores = threads)
        script += make_khmer_script(self.get_all_reads(), path + "normalised_reads.fastq", threads)
        with open(path +  "diginorming.sh","w") as handle:
            handle.writelines(script)
        if submit :
            job = sh.sbatch(path + "diginorming.sh")
            if not self.job_ids.has_key('diginorming'):
                self.job_ids['diginorming'] = job.split()[-1]
            print "Launched diginorming on", self.name
            with open(self.job_file, 'w') as handle:
                json.dump(self.job_ids, handle)

    

#mendota = Assembly("/home/moritz/repos/assemblies/mendota.json","b2011138")
#troutbog = Assembly("/home/moritz/repos/assemblies/troutbog.json","b2011138")
#test = Assembly("/home/moritz/repos/assemblies/troutbog_onesample_test.json","b2011138")
#cami_hc = Assembly("/home/moritz/repos/assemblies/cami_h.json","b2013086")
#oilsands = Assembly("/home/moritz/repos/assemblies/oilsands.json","b2011032")
#deepbio_2 = Assembly("/home/moritz/repos/assemblies/deepbio_2.json","b2013086")
#deepbio_1 = Assembly("/home/moritz/repos/assemblies/deepbio_1.json","b2013086")

#tbe6 = Sample("TB-E6", [{ 'I' : "/home/moritz/b2011138/acI_enrichments/acI_TB-E6.fastq"}], troutbog.sample_path)

#TARA  = Assembly("/home/moritz/people/tara/TARA.json","b2014083")
