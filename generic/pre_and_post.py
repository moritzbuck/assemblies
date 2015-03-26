from Bio import SeqIO, SeqUtils
from tqdm import tqdm
from pandas import DataFrame
import os

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
    
def filter_fasta_by_length(infasta, outfasta, length = 1000):
    with open(infasta) as handle:
        seqs = [s for s in SeqIO.parse(handle,"fasta") if len(s) >= length]
        
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

def make_bins(path, concoct_file, assembly):
        contig_dict = DataFrame.from_csv(concoct_file, header=None)[1].to_dict()
        bin_list = list(set(contig_dict.values()))
        bins = { b: [] for b in bin_list}
        with open(assembly) as handle:
                for s in tqdm(SeqIO.parse(handle,"fasta")):
                        if contig_dict.has_key(s.id):
                                bins[contig_dict[s.id]] += [s]

        if not os.path.exists(path):
                os.makedirs(path)

        pads = len(str(max([len(c) for c in bins.values()])))
        bpads = len(str(len(bin_list)))
        name_dict ={}
        for i,b in tqdm(bins.iteritems()):
            for j,s in enumerate(b):
                name_dict[s.id] = "bin_" + str(i).zfill(bpads) + "_"  + str(j).zfill(pads)
                s.id = "bin_" + str(i).zfill(bpads) + "_"  + str(j).zfill(pads)
                s.description = ""
                s.name = s.id
#            with open(path + "bin_" + str(i).zfill(bpads) +".fasta","w") as handle :
#                    SeqIO.write(b,handle,"fasta")
        return name_dict
