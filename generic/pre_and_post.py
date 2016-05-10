from Bio import SeqIO, SeqUtils
from tqdm import tqdm
from pandas import DataFrame
from pandas import Index
import os

def clean_metabat_out(name, folder):
    list = [f for f in os.listdir(folder) if f[-3:] == ".fa"]
    list = { int(l.split(".")[3]) : l for l in list}
    fill = len(str(len(list)))
    for id,f in list.iteritems():
        clean_fasta_names(folder + f, folder + name +  "_bin_" + str(id).zfill(fill) + ".fasta",name + "_bin_" + str(id).zfill(fill) + "_")


def clean_fasta_names(infasta, outfasta, prefix = "seq_"):
    with open(infasta) as handle:
        seqs = [s for s in SeqIO.parse(handle,"fasta")]

    pads = len(str(len(seqs)))

    for i,s in enumerate(seqs):
        s.id = prefix + str(i).zfill(pads)
        s.description = ""
        s.name = s.id
        
    with open(outfasta,"w") as handle:
        SeqIO.write(seqs,handle,"fasta")
    
def filter_fasta_by_length(infasta, outfasta, min_length = 1000, max_length=float("inf")):
    with open(infasta) as handle:
        seqs = [s for s in SeqIO.parse(handle,"fasta") if len(s) >= min_length and len(s) < max_length]
        
    with open(outfasta,"w") as handle:
        SeqIO.write(seqs,handle,"fasta")
    

def coverage_merging(infasta, cover_dict, gff = None, ofile = None):
        with open(infasta) as handle:
            print "Opening contigs"
            covs = {s.id: dict({k: 0  for k in cover_dict.keys()}.items() + [('length',len(s)), ('GC', SeqUtils.GC(s.seq))]) for s in tqdm(SeqIO.parse(handle, "fasta")) }
        print "Opening", "coverage files"
        for k,v in cover_dict.iteritems():
            print "parsing", k
            with open(v) as handle:
                for l in tqdm(handle):
                    temp = l.split()
                    if covs.has_key(temp[0]):
                        covs[temp[0]][k] += int(temp[2])

        #contig_wise_dict = { c : {s: mean(v) for s,v in x.iteritems()} for c,x in covs.iteritems() }
        
#        t = DataFrame.from_dict(covs).transpose()
#        if ofile: t.to_csv(ofile)
        return covs

def make_bins(path, concoct_file, assembly, coverage_file, from_canopy = False, prefix = "bin", min_size=None):
        if from_canopy:
                contig_dict = DataFrame.from_csv(concoct_file, sep="\t", index_col=1, header=None)[0].to_dict()
        else :
                contig_dict = DataFrame.from_csv(concoct_file, header=None)[1].to_dict()
                
        mapi = DataFrame.from_csv(coverage_file, sep="\t")
        bin_list = list(set(contig_dict.values()))
        bins = { b: [] for b in bin_list}
        bins['non-bin'] = []
        with open(assembly) as handle:
                for s in tqdm(SeqIO.parse(handle,"fasta")):
                        if contig_dict.has_key(s.id):
                            bins[contig_dict[s.id]] += [s]
                        else:
                            bins['non-bin'] += [s]

        if not os.path.exists(path):
                os.makedirs(path)

        pads = len(str(max([len(c) for c in bins.values()])))
        bpads = len(str(len(bin_list)))
        name_dict ={}
        for i,b in tqdm(bins.iteritems()):
            if from_canopy:
                b_name = i 
            else:
                b_name = prefix + str(i).zfill(bpads) 
            for j,s in enumerate(b):
                name_dict[s.id] = b_name + "_"  + str(j).zfill(pads)
                s.id = name_dict[s.id]
                s.description = ""
                s.name = s.id

            if min_size:
                if sum([len(s.seq) for s in b]) > min_size:
                    with open(path + b_name +".fasta","w") as handle :
                        SeqIO.write(b,handle,"fasta")
                else:
                    bins['non-bin'] += b
            else:
                with open(path + b_name +".fasta","w") as handle :
                    SeqIO.write(b,handle,"fasta")

        if min_size:
                with open(path + prefix + "non-bin.fasta","w") as handle :
                        SeqIO.write(bins["non-bin"],handle,"fasta")

            
        mapi.index = Index([name_dict[v] for v in mapi.index])
        mapi.to_csv(path + "all_bins_mapping.csv")
        mapi.loc[[v for v in  mapi.index if not "non-bin" in v]].to_csv(path + "yes_bin_mapping.csv")
        [",".join([v,k]) for k,v in tqdm(name_dict.iteritems())]

        kk = set(name_dict.keys())
        nbined = [k for k in tqdm(kk) if name_dict[k] in kk]
        for k in  tqdm(nbined):
            t=name_dict[k]
            name_dict[k]=name_dict[name_dict[k]]
            del name_dict[t]
        return (mapi, name_dict)


def parse_all_gffs_for_ecs(annot_folder):
    all_gffs = [g[:-1] for g in sh.find(annot_folder) if ".gff" in g and not "genes.gff" in g]
    temp_dict = {}
    for gff in tqdm(all_gffs):
        bin = gff.split("/")[-3]
        with open(gff) as handle:
            lines = { l.split("\t")[0] : [e.split("=")[1] for e in l.split("\t")[8].split(";") if "eC_number" in e ][0]  for l in handle.readlines() if "eC_number" in l}
        temp_dict[bin] = lines
    all_ecs = set(sum([vs.values() for vs in temp_dict.values()],[]))
    out_dict ={}
    for gff in tqdm(all_gffs):
        bin = gff.split("/")[-3]
        tt = temp_dict[bin]
        out_dict[bin] = {}
        for e in all_ecs:
            out_dict[bin][e] = ";".join([g for g,ec in temp_dict[bin].iteritems() if ec == e] )
    return DataFrame.from_dict(out_dict)
