from Bio import SeqIO, SeqUtils
from tqdm import tqdm 

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
