from __future__ import division

import sh
import os
import json
from  generic.scripts import *
from tqdm import tqdm
from pandas import DataFrame


class Sample(object):
    def __repr__(self):
        return '<%s object %s with %i libraries %i of which are paired>' % (self.__class__.__name__, self.name, len(self.reads), sum([len(r) == 2 for r in self.reads]) )        
        
    def __init__(self, name, reads, path, metadata = None):
        self.name = name
        self.reads = reads
        
        assert len([os.path.basename(r) for r in sum(self.reads, [])]) == len(set([os.path.basename(r) for r in sum(self.reads, [])])), "Don' t make my life a pain, and name your fastq with distinct names..."
        
        self.path = path + self.name + "/"
        self.metadata = metadata
        self.clean_path = self.path + "clean/"

        if not os.path.exists(self.clean_path):
            os.makedirs(self.clean_path)


    def get_pairs(self, clean = True):
        pairs = [p for p in self.reads if len(p) == 2]
        if clean:
            fwd = [self.clean_path + os.path.basename(p[0]).replace(".gz","") for p in pairs]
            rev = [self.clean_path + os.path.basename(p[1]).replace(".gz","") for p in pairs]
            single = [self.clean_path + "single_" + os.path.basename(p[0]).replace(".gz","") for p in pairs]
            return (rev, fwd, single)
        else :
            fwd = [p[0] for p in pairs]
            rev = [p[1] for p in pairs]
            return (rev, fwd)

    def get_singles(self, clean = True):
        pairs = [p for p in self.reads if len(p) == 2]
        singles = [p for p in self.reads if len(p) == 1]
        if clean:
            single = [self.clean_path + os.path.basename(p[0]).replace(".gz","") for p in singles]
            single_pairs = [self.clean_path + "single_" + os.path.basename(p[0]).replace(".gz","") for p in pairs]
            return (single, single_pairs)
        else :
            single = [p[0] for p in singles]
            return (single, [])

    def get_all_reads(self, clean = True):
        all_reads = sum(self.reads,[])
        if clean:
            all_reads = [self.clean_path + os.path.basename(p).replace(".gz","") for p in all_reads]
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
            sample_dict[s.name]['1'] = paired[0]
            sample_dict[s.name]['2'] = paired[1]
            sample_dict[s.name]['U'] = unpaired[0] + unpaired[1]
        return  sample_dict
    
    def clean_all(self, submit = False):
        script = make_slurm_header(self.sample_path, "clean_all", self.project)
        
        singles_in = self.get_all_reads( clean = False, kind = "single" )
        singles_out = self.get_all_reads( clean = True, kind = "single" )
        pairs_in = self.get_all_reads( clean = False, kind = "paired" )
        pairs_out = self.get_all_reads( clean = True, kind = "paired" )

        if len(singles_in) >0 :
            script += make_parallel_sickle_script(singles_in, singles_out, paired=False)

        if len(pairs_in[0]) >0 :
            script += make_parallel_sickle_script(pairs_in, pairs_out, paired=True)

        with open(self.sample_path +  "clean_all_script.sh","w") as handle:
            handle.writelines(script)
        if submit :
            clean_job = sh.sbatch(self.sample_path +  "clean_all_script.sh")
            self.job_ids['clean_all'] = clean_job.split()[-1]
            print "Launched clean_all"
            with open(self.job_file, 'w') as handle:
                json.dump(self.job_ids, handle)

                

    def megahit_assembly(self, name="megahit", max_read_len = 290, submit = False):
        if self.job_ids.has_key('clean_all'):
            deps = [self.job_ids['clean_all']]
        else:
            deps = ['1']
        path = self.raw_ass + name +"/"
        script = make_slurm_header(path, "assembly_" + name, self.project, constraint = "fat", time = "2-00:00:00", deps = deps)
        script += make_megahit_script(self.get_all_reads(), path)
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
        script += make_mapping_script(reference, path, self.get_clean_sample_dict(), only_bowtie = bowtie_only)
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


mendota = Assembly("mendota.json","b2011138")
troutbog = Assembly("troutbog.json","b2011138")
