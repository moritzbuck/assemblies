import os

def make_slurm_header(path,name = "Monkey" , proj = "b2011138", time = "2-00:00:00", cores = 16, constraint = None, deps = None):
    if cores < 16:
        node_type = "core"
    else:
        node_type = "node"
    mem_size = {"fat": "mem512GB", "medium":  "mem256GB", "slim": "mem128GB" }
    
    if not os.path.exists(path):
        os.makedirs(path)

    if constraint == None :
        constraint_sbatch = ""
    else :
        constraint_sbatch = "\n#SBATCH -C " + mem_size[constraint]        

    if deps == None :
        deps_sbatch = ""
    else :
        deps_sbatch = "\n#SBATCH -d afterok:" + ":".join(deps)
    
    raw_slurm = """#!/bin/bash
#SBATCH -D %s
#SBATCH -J %s
#SBATCH -o %s
#SBATCH -e %s
#SBATCH -A %s
#SBATCH -t %s
#SBATCH -n %d
#SBATCH -p %s%s%s
#SBATCH --mail-user murumbii@gmail.com
#SBATCH --mail-type=ALL

""" % (path, name, path + name + ".err",  path + name + ".out", proj, time, cores, node_type, constraint_sbatch, deps_sbatch)

    return raw_slurm



raw_concoct = """
concoct --coverage_file  %s --composition_file %s -l %s -b %s --no_total_coverage
"""

raw_checkm = """
checkm lineage_wf -t 16 -x fasta . %s/checkm > %scheckm.txt
"""            


raw_hmm_pipe = """
module load hmmer
hmmsearch --cpu 16 %s %s > %s
""" 

raw_ray_single = """
module load openmpi
mpiexec -n 16 Ray231 -k %d -o %s  %s
"""



raw_newbler = """
/home/moritz/repos/metassemble/scripts/assembly/merge-asm-newbler.sh %s %s
"""

nr_threads = 16
executable = 'Ray231'


def make_parallel_sickle_script(library_list, threads = 16):
        interleaved =  [l for l in library_list if l.has_key('I')]
        paired = [l for l in library_list if l.has_key('1')]
        singles = [l for l in library_list if l.has_key('U')]

        raw = ""
        
        if len(interleaved) > 0 :
            raw += "\nparallel --xapply -j %i  split-paired-reads.py {1} ::: %s " %  (threads, " ".join([p['I'] for p in interleaved]))
            raw += "\nparallel --xapply -j %i sickle  pe -t sanger -f {1} -r {2} -o {3} -p {4} -s {5} ::: %s :::  %s ::: %s ::: %s ::: %s" % (threads, " ".join([os.path.basename(p['I']) + ".1" for p in interleaved]), " ".join([os.path.basename(p['I']) +".2" for p in interleaved]), " ".join([p['c1'] for p in interleaved]), " ".join([p['c2'] for p in interleaved]), " ".join([p['cU'] for p in interleaved]))
            raw += "\nparallel --xapply -j %i  rm {1} ::: %s " %  (threads, " ".join([os.path.basename(p['I']) +".1" for p in interleaved] + [os.path.basename(p['I']) + ".2" for p in interleaved]))
            
        if len(paired) > 0 :
            raw += "\nparallel --xapply -j %i sickle  pe -t sanger -f {1} -r {2} -o {3} -p {4} -s {5} ::: %s ::: %s ::: %s ::: %s ::: %s" % (threads, " ".join([p['1'] for p in paired]), " ".join([p['2'] for p in paired]), " ".join([p['c1'] for p in paired]), " ".join([p['c2'] for p in paired]), " ".join([p['cU'] for p in paired]))
            
        if len(singles) > 0 :
            raw += "\nparallel --xapply --j %i sickle  se -t sanger -g -f {1} -o {2} ::: %s ::: %s" % (threads, " ".join([p['U'] for p in singles]), " ".join([p['cU'] for p in singles]))
            
        return raw

def make_megahit_script(inreads, outfolder, mem = 500e9, max_read_len = 290):
    raw = """
/home/moritz/repos/megahit/megahit --cpu-only -m %d -l %d -o %s --input-cmd %s
ln -s %sfinal.contigs.fa %scontigs.fasta
"""
    raw = raw % (mem, max_read_len, outfolder, "\"cat " + " ".join(inreads) + "\"", outfolder,outfolder)
    return raw

def make_minimus_script(infasta, outfasta, maxtrim = 20):
    raw = """
toAmos -s %s -o minimus
minimus2 minimus -D %d
cat minimus.fasta minimus.singletons.seq > %s
rm minimus*
"""
    raw = raw % (infasta, maxtrim, outfasta)
    return raw

def make_mapping_script(infasta, path, sample_dictionary, only_bowtie = True, remove_duplicates=False, threads = 16 ):

    #sample dictionary is a dictionary of dictionaries, keys are the sample names, the dictionary contains list of fastq paths ['1' for first of paris. '2' second of pairs, and 'U' for unpaired fastqs]
    
    script = ""
    if not os.path.exists(infasta + ".1.bt2"):
        script += "bowtie2-build %s %s\n" % (infasta, infasta)

    for k,d in sample_dictionary.iteritems():
        unmapped = path + k + "/" + k + "_unmapped.fastq"
        sam = path + k + "/" + k + "_map.sam"
        err = path + k + "/" + k + "_map.stats"
        out = path + k + "/" + k + "_map.err"
        if len(d['1']) > 0 :
            pair1 = "-1 "+ ",".join(d['1'])
            pair2 = "-2 "+ ",".join(d['2'])
        else :
            pair1 = ""
            pair2 = ""
        if len(d['U']) >0:
            unpair = "-U " + ",".join(d['U'])
        else:
            unpair = ""
        if not os.path.exists(path + k ):
            os.makedirs(path + k )
        script += "\nbowtie2 -p %d -x %s %s %s %s  -S %s 2> %s > %s" % (threads, infasta, pair1, pair2, unpair, sam, err, out)

    script += "\n"

    if not only_bowtie:
        if not os.path.exists(infasta + ".fai"):
            script += "samtools faidx %s\n" % (infasta)
        sample_names = " ".join(sample_dictionary.keys())
        script += "parallel --xapply -j %d samtools view -bT %s %s{1}/{1}_map.sam '>' %s{1}/{1}_map.bam ::: %s\n" %( threads, infasta, path, path, sample_names)
        script += "parallel --xapply -j %d samtools sort %s{1}/{1}_map.bam %s{1}/{1}_map_sorted ::: %s\n" %( threads, path, path, sample_names)
        script += "parallel --xapply -j %d samtools index %s{1}/{1}_map_sorted.bam ::: %s\n" %( threads, path, sample_names)

        final_bam = "_map_sorted.bam"
        
        if remove_duplicates:
            script += "\nmodule load picard/1.92\n "
            java_opts = " -Xms2g -Xmx12g -XX:MaxPermSize=2g -XX:+CMSClassUnloadingEnabled "
            java_jar = " -jar $PICARD_TOOLS_DIR/MarkDuplicates.jar "
            java_MD_opts = " AS=TRUE VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 REMOVE_DUPLICATES=TRUE "
            io = " INPUT=%s{1}/{1}_map_sorted.bam INPUT=%s{1}/{1}_map_no_dups.bam :::  %s\n" %(path, path, sample_names)
            script += "parallel --xapply -j %d samtools java " %(threads) + java_opts + java_jar + java_MD_opts + io
            script += "parallel --xapply -j %d samtools sort %s{1}/{1}_map_no_dups.bam %s{1}/{1}_map_no_dups_sorted ::: %s\n" %( threads, path, path, sample_names)
            script += "parallel --xapply -j %d samtools index %s{1}/{1}_map_no_dups_sorted.bam ::: %s\n" %( threads, path, sample_names)
            script += "parallel --xapply -j %d samtools flagstat %s{1}/{1}_map_no_dups_sorted.bam '>' %s{1}/{1}_map.flagstat ::: %s\n" %( threads, path, path, sample_names)

            final_bam = "_map_no_dups_sorted.bam"

        script += "parallel --xapply -j %d genomeCoverageBed -ibam  %s{1}/{1}%s '>' %s{1}/{1}.coverage ::: %s\n" %( threads, path, final_bam, path, sample_names)
        script += "parallel --xapply -j %d python $METASSEMBLE_DIR/scripts/validate/map/gen_contig_cov_per_bam_table.py --isbedfiles %s %s{1}/{1}.coverage '>' %s{1}/{1}.coverage.percontig ::: %s\n" %( threads, infasta, path, path, sample_names)
        script += "\ncd %s\ncut -f1,2,3  %s/*.percontig  > tmp.txt" % (path,sample_names.split()[0])
        script += "\nfor l in `ls  %s`; do  cut -f4  $l | paste tmp.txt - > tmp2.txt; mv tmp2.txt tmp.txt ; done" %("/*.percontig ".join(sample_names.split()[::-1]) + "/*.percontig")
        script += """
echo "%s" >  %s/mapping.tsv
tail -n +2 tmp.txt >>  %s/mapping.tsv
""" % ("\t".join(["contig", "length", "GC"] + sample_names.split()), path, path)

    
    return script

def make_concoct_script(covfile, infasta, path, cutoff = 1000, kmer = 4, max_bins = 1000):
    script = """
cut -f1,4- %s > tmp
concoct --coverage_file  tmp  --composition_file %s -k %i -l %i -b %s -c %i --no_total_coverage
rm tmp
""" %(covfile, infasta, kmer, cutoff, path + "concoct_%i_kmer_%i/" %(cutoff, kmer), max_bins)
    return script


def make_bin_bmfa(infasta,  path, name, threads = 16, coverages = None, hmmdb="~/glob/data/pfam/Pfam-A.hmm"):
    raw = """module unload gcc
module load intel
module load hmmer
"""
    raw += """\nprokka --outdir %s/prokka/ --cpus %i --prefix %s --locustag %s --force %s""" % (path, threads, name, name,  infasta)
    raw += """\n#hmmsearch --cpu %d %s %s > %s""" % (threads, hmmdb, infasta, path + os.path.basename(hmmdb) +".raw")
    raw += """\ncheckm lineage_wf -t %i -x fna %s/prokka  %s/checkm > %s/checkm.txt""" % (threads, path, path, path)
    return raw



def make_besst_script(reference, path, sample_dictionary, maps_path):
    print "check for sorting of BAM-files"
    maps = [maps_path + k + "/" + k + "_map_sorted.bam" for k,d in sample_dictionary.iteritems() if len(d['1']) > 0 ]
    raw = """
runBESST  --orientation %s -d -z %s -e %s -w 1 -T %s -k %s -c %s -f %s -o %s
"""
    Ts = ["10000"] * len(maps)
    ks = ["250"] * len(maps)
    oris = " ".join(["fr"] *len(maps))
    zs = " ".join(["1000000"] *len(maps))
    es = " ".join(["2"] *len(maps))    
    #runBESST --orientation fr -d -z 100000 -e 1 -T 10000 -k 250 -w 1 -c megahit_out/final.contigs.fa -f mapped_sorted.bam -o besst_scaff_e_1
    raw = raw % (oris, zs, es, " ".join(Ts), " ".join(ks), reference, " ".join(maps) , path)
    return raw


def make_reconcoct_script(covfile, infasta, prefix, cutoff = 1000, kmer = 4):
    script = """
python /home/moritz/repos/Pyscratches/20150727_reconcoct/reconcoct.py -m %s -p %s %s %s
"""  %(str(cutoff), prefix, covfile, infasta)
    return script

def make_khmer_script(reads_in, reads_out, threads):
    script = """
normalize-by-median.py -k 20 -C 20 -N 4 -x 8e9 -s normC20k20.kh -o temp.fastq <(cat %s) 
filter-abund.py --threads %d -V normC20k20.kh temp.fastq -o temp.abundfilt.fastq
normalize-by-median.py -C 5 -k 20 -N 4 -x 8e9 temp.abundfilt.fastq -o %s
rm normC20k20.kh temp.abundfilt.fastq temp.fastq 
""" %(" ".join(reads_in),threads, reads_out)

    return script
    



def make_star_mapping_script(infasta, path, sample_dictionary, threads = 16 ):

    #sample dictionary is a dictionary of dictionaries, keys are the sample names, the dictionary contains list of fastq paths ['1' for first of paris. '2' second of pairs, and 'U' for unpaired fastqs]

    index = infasta + "_star_index/"
    script = ""
    if not os.path.exists(index):
        script += "mkdir %s" % (index)
        script += "/home/moritz/repos/STAR/source/STAR --runMode genomeGenerate --runThreadN %s --genomeDir %s --genomeFastaFiles %s --genomeChrBinNbits 13 --limitGenomeGenerateRAM 500397404373 \n" % (str(threads), index , infasta)
        
    for k,d in sample_dictionary.iteritems():
        print "Need to change code to all reads it bugs somewhere"
        unmapped = path + k + "/" + k + "_unmapped.fastq"
        sam = path + k + "/" + k + "_map"
        err = path + k + "/" + k + "_map.stats"
        out = path + k + "/" + k + "_map.err"
        if len(d['1']) > 0 :
            pair1 = "-1 "+ ",".join(d['1'])
            pair2 = "-2 "+ ",".join(d['2'])
        else :
            pair1 = ""
            pair2 = ""
        if len(d['U']) >0:
            unpair = ",".join(d['U'])
        else:
            unpair = ""
        if not os.path.exists(path + k ):
            os.makedirs(path + k )
        script +="\necho Mapping sample %s\n/home/moritz/repos/STAR/source/STAR --genomeLoad LoadAndKeep --limitBAMsortRAM 16000000000 --runThreadN %d --genomeDir %s --outSAMtype BAM SortedByCoordinate --readFilesIn %s --outFileNamePrefix %s 2> %s > %s" % (k, threads, index, unpair, sam, err, out)

    sample_names = " ".join(sample_dictionary.keys())

    script += "\n"

    script += "parallel --xapply -j %d samtools index %s{1}/{1}_mapAligned.sortedByCoord.out.bam ::: %s\n" %( threads, path, sample_names)
    script += "parallel --xapply -j %d genomeCoverageBed -ibam  %s{1}/{1}_mapAligned.sortedByCoord.out.bam '>' %s{1}/{1}.coverage ::: %s\n" %( threads, path, path, sample_names)
    script += "parallel --xapply -j %d python $METASSEMBLE_DIR/scripts/validate/map/gen_contig_cov_per_bam_table.py --isbedfiles %s %s{1}/{1}.coverage '>' %s{1}/{1}.coverage.percontig ::: %s\n" %( threads, infasta, path, path, sample_names)
    script += "\ncd %s\ncut -f1,2,3  %s/*.percontig  > tmp.txt" % (path,sample_names.split()[0])
    script += "\nfor l in `ls  %s`; do  cut -f4  $l | paste tmp.txt - > tmp2.txt; mv tmp2.txt tmp.txt ; done" %("/*.percontig ".join(sample_names.split()[::-1]) + "/*.percontig")
    script += """
echo "%s" >  %s/mapping.tsv
tail -n +2 tmp.txt >>  %s/mapping.tsv
""" % ("\t".join(["contig", "length", "GC"] + sample_names.split()), path, path)

    
    return script





def make_bbmapping_script(infasta, path, sample_dictionary, only_bowtie = True, remove_duplicates=False, threads = 16 ):

    #sample dictionary is a dictionary of dictionaries, keys are the sample names, the dictionary contains list of fastq paths ['1' for first of paris. '2' second of pairs, and 'U' for unpaired fastqs]
    
    script = "module load bbmap\n"
#    previous_refs = os.popen(" ".join(["find",os.path.dirname(infasta)+"/ref/", "-name" , "summary.txt"])).read().split()

#    if len(previous_refs) == 0 or (int(os.popen(" ".join([ "grep", "-c", infasta, " ".join(previous_refs)])).read()[:-1]) < 1):
#        script += "bbmap.sh ref=%s\n" % (infasta)

    for k,d in sample_dictionary.iteritems():
        unmapped = path + k + "/" + k + "_unmapped.fastq"
        sam = path + k + "/" + k + "_map.sam"
        err = path + k + "/" + k + "_map.stats"
        out = path + k + "/" + k + "_map.err"
        if len(d['1']) > 0 :
            pair1 = ",".join(d['1'])
            pair2 = ",".join(d['2'])
            script += "\nbbmap.sh ref=%s in=%s in2=%s out=%s 2> %s > %s \n" %(infasta, pair1, pair2, sam, err, out)
        elif len(d['U']) >0:
            unpair = " ".join(d['U'])
            script += "\ncat %s | bbmap.sh ref=%s in=stdin out=%s 2> %s > %s \n" %( unpair, infasta, sam, err, out)
        if not os.path.exists(path + k ):
            os.makedirs(path + k )

    script += "\n"

    if not only_bowtie:
        if not os.path.exists(infasta + ".fai"):
            script += "samtools faidx %s\n" % (infasta)
        sample_names = " ".join(sample_dictionary.keys())
        script += "parallel --xapply -j %d samtools view -bT %s %s{1}/{1}_map.sam '>' %s{1}/{1}_map.bam ::: %s\n" %( threads, infasta, path, path, sample_names)
        script += "parallel --xapply -j %d samtools sort %s{1}/{1}_map.bam %s{1}/{1}_map_sorted ::: %s\n" %( threads, path, path, sample_names)
        script += "parallel --xapply -j %d samtools index %s{1}/{1}_map_sorted.bam ::: %s\n" %( threads, path, sample_names)

        final_bam = "_map_sorted.bam"
        
        if remove_duplicates:
            script += "\nmodule load picard/1.92\n "
            java_opts = " -Xms2g -Xmx12g -XX:MaxPermSize=2g -XX:+CMSClassUnloadingEnabled "
            java_jar = " -jar $PICARD_TOOLS_DIR/MarkDuplicates.jar "
            java_MD_opts = " AS=TRUE VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 REMOVE_DUPLICATES=TRUE "
            io = " INPUT=%s{1}/{1}_map_sorted.bam INPUT=%s{1}/{1}_map_no_dups.bam :::  %s\n" %(path, path, sample_names)
            script += "parallel --xapply -j %d samtools java " %(threads) + java_opts + java_jar + java_MD_opts + io
            script += "parallel --xapply -j %d samtools sort %s{1}/{1}_map_no_dups.bam %s{1}/{1}_map_no_dups_sorted ::: %s\n" %( threads, path, path, sample_names)
            script += "parallel --xapply -j %d samtools index %s{1}/{1}_map_no_dups_sorted.bam ::: %s\n" %( threads, path, sample_names)
            script += "parallel --xapply -j %d samtools flagstat %s{1}/{1}_map_no_dups_sorted.bam '>' %s{1}/{1}_map.flagstat ::: %s\n" %( threads, path, path, sample_names)

            final_bam = "_map_no_dups_sorted.bam"

        script += "parallel --xapply -j %d genomeCoverageBed -ibam  %s{1}/{1}%s '>' %s{1}/{1}.coverage ::: %s\n" %( threads, path, final_bam, path, sample_names)
        script += "parallel --xapply -j %d python $METASSEMBLE_DIR/scripts/validate/map/gen_contig_cov_per_bam_table.py --isbedfiles %s %s{1}/{1}.coverage '>' %s{1}/{1}.coverage.percontig ::: %s\n" %( threads, infasta, path, path, sample_names)
        script += "\ncd %s\ncut -f1,2,3  %s/*.percontig  > tmp.txt" % (path,sample_names.split()[0])
        script += "\nfor l in `ls  %s`; do  cut -f4  $l | paste tmp.txt - > tmp2.txt; mv tmp2.txt tmp.txt ; done" %("/*.percontig ".join(sample_names.split()[::-1]) + "/*.percontig")
        script += """
echo "%s" >  %s/mapping.tsv
tail -n +2 tmp.txt >>  %s/mapping.tsv
""" % ("\t".join(["contig", "length", "GC"] + sample_names.split()), path, path)

    
    return script


def make_single_bbmapping_script(infasta, path, reads , only_bowtie = False, remove_duplicates=False, threads = 16 ):

    
    script = "module load bbmap\n"

    tpath = "/tmp/"
    sam = tpath +  "/map.sam"
    err = path +  "/map.stats"
    out = path +  "/map.err"

    #script += "\nbbmap.sh ref=%s in=%s in2=%s out=%s 2> %s > %s \n" %(infasta, pair1, pair2, sam, err, out)
    script += "\nunpigz -c  %s | bbmap.sh threads=%s nodisk=t ref=%s in=stdin out=%s 2> %s > %s \n" %(reads, str(threads), infasta, sam, err, out)
    if not os.path.exists(path ):
        os.makedirs(path)
    script += "grep -E \"@SQ|@HD|@PG|scaff|Chr\" %s/map.sam > %s/map_hits.sam " %( tpath, path if only_bowtie else tpath )
#    script += "sambamba_v0.6.0 view  -t %s -S -f bam  -T %s %s/map.sam > %s/map.bam \n" %( str(threads), infasta, tpath, path if only_bowtie else tpath )

    if not only_bowtie:

        if not os.path.exists(infasta + ".fai"):
            script += "samtools faidx %s\n" % (infasta)

        script += "samtools sort %s/map.bam %s/map_sorted \n" %(  tpath, tpath )
        script += "samtools index %s/map_sorted.bam \n" %(  tpath )

        final_bam = "map_sorted.bam"

        if remove_duplicates:
            script += "\nmodule load picard/1.92\n "
            java_opts = " -Xms2g -Xmx12g -XX:MaxPermSize=2g -XX:+CMSClassUnloadingEnabled "
            java_jar = " -jar $PICARD_TOOLS_DIR/MarkDuplicates.jar "
            java_MD_opts = " AS=TRUE VALIDATION_STRINGENCY=LENIENT MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 REMOVE_DUPLICATES=TRUE "
            io = " INPUT=%s/map_sorted.bam INPUT=%s/map_no_dups.bam \n" %(tpath, tpath )
            script += "samtools java " %(threads) + java_opts + java_jar + java_MD_opts + io
            script += "samtools sort %s/map_no_dups.bam %s/map_no_dups_sorted \n" %(  tpath, tpath )
            script += "samtools index %s/map_no_dups_sorted.bam \n" %(  tpath )
            script += "samtools flagstat %s/map_no_dups_sorted.bam > %s/map.flagstat \n" %(  tpath, tpath )

            final_bam = "map_no_dups_sorted.bam"

        script += "genomeCoverageBed -ibam  %s/%s > %s/coverage.tsv\n" %(  tpath, final_bam, path )
        script += "python $METASSEMBLE_DIR/scripts/validate/map/gen_contig_cov_per_bam_table.py --isbedfiles %s %s/coverage.tsv > %s/coverage.percontig.tsv \n" %(  infasta, path, path )
        script += "rm %s/map.sam"  %(tpath)
        script += "rm %s/map.bam" %(tpath)
    return script
