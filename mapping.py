#!/usr/bin/python

'''
Mapping reads from Atf3 paper to upload to GEO (Allison et al. Diabetes 2018)
'''

import subprocess
import glob
import datetime
import re
import os
import pandas as pd
import numpy as np

# mkdir function
def mkdir(directory):
    out = f"mkdir {directory}"
    print(out)
    subprocess.call(out, shell=True)

# grab software version
def grab_version(program):
    arg = f"{program} --version"
    output = subprocess.Popen(arg, stdout=subprocess.PIPE, shell=True).communicate()[0]
    output = output.decode('utf-8')
    return output

def combine(sample):
    '''combine FASTQ files into a single file per sample'''
    files = glob.glob(f"Run*/Project_mmyers/{sample}/*.fastq.gz")
    arg = f"cat {' '.join(files)} > {sample}.fastq.gz"
    print(arg)
    subprocess.call(arg, shell=True)

def filter(fastq):
    '''pipe FASTQ files through zcat and fastq_quality_filter'''
    out = f"zcat {fastq} | "
    out += "fastq_quality_filter -q 20 -z "
    out += f"-o filtered/{fastq}"
    time = f"{datetime.datetime.now().strftime('%c')}"
    print(f"{time}\n{out}\n")
    subprocess.call(out, shell=True)

def generate_genome(fasta, gtf):
    '''generate STAR genome index'''
    f"STAR {grab_version('star')}"
    if not os.path.isdir("data/genome"):
        mkdir('data/genome')
    out = "star "
    out += "--runMode genomeGenerate "
    out += "--genomeDir data/genome "
    out += "--runThreadN 8 "
    out += f"--genomeFastaFiles {fasta} "
    out += f"--sjdbGTFfile {gtf} "
    time = f"{datetime.datetime.now().strftime('%c')}"
    print(f"{time}\n{out}\n")
    subprocess.call(out, shell=True)

def count(fastq, gtf):
    '''map FASTQ files with STAR and generate count files'''
    sample = re.findall("Sample_[0-9]+", fastq)[0]
    out = "star "
    out += "--runThreadN 8 "
    out += "--genomeDir data/genome "
    out += f"--sjdbGTFfile {gtf} "
    out += "--readFilesCommand zcat "
    out += f"--readFilesIn {fastq} "
    out += f"--outFileNamePrefix data/STAR/{sample} "
    out += "--outSAMtype None "
    out += "--quantMode GeneCounts"
    # grab date and time
    time = f"{datetime.datetime.now().strftime('%c')}"
    print(f"{time}\n{out}\n")
    subprocess.call(out, shell=True)

def gunzip(file):
    out = f"gunzip {file}"
    print(out)
    subprocess.call(out, shell=True)

def gzip(file):
    out = f"gzip {file}"
    print(out)
    subprocess.call(out, shell=True)

def tx2gene(gtf):
    '''generate a tx2gene file from a GTF file'''
    print(f"\n----- Generating tx2gene CSV file from {gtf} -----\n")
    fields = ['transcript_id', 'transcript_name', 'gene_id', 'gene_name', 'gene_biotype']
    transcript_id = np.where(np.array(fields) == 'transcript_id')[0][0]
    # function to grab info associated with user input fields
    def grab_info(field):
        info = re.findall(field + ' \"(.*?)\"', line)
        if len(info) == 0:
            info = ""
        else:
            info = info[0]
        return info
    # generate dictionary from GTF file
    d = dict()
    with open(re.sub('.gz$', '', gtf), 'r') as f:
        for line in f:
            if line.startswith('#'): continue # skip header lines
            line = line.split(sep='\t')[8] # line with relevant info
            field_stash = list()
            for field in fields:
                hit = grab_info(field)
                field_stash.append(hit)
            # if the transcript_id is new, add it to the dictionary with its info
            if field_stash[transcript_id] not in d.keys():
                d[field_stash[transcript_id]] = np.array(field_stash)[np.arange(len(field_stash)) != transcript_id]
            else: continue
    # turn into data frame for easy output
    df = pd.DataFrame.from_dict(d, orient='index')
    df.reset_index(level=0, inplace=True)
    rename_col = {'index': 'transcript_id'}
    count = 0
    for field in fields:
        if field == 'transcript_id': continue
        else:
            rename_col[count] = field
            count += 1
    df.rename(columns=rename_col, inplace=True)
    # write to CSV file
    print("Saving file as data/tx2gene.csv")
    df.to_csv("data/tx2gene.csv", index=False)

def mv(file, destination):
    out = f"mv {file} {destination}"
    print(out)
    subprocess.call(out, shell=True)

def rm(dir, recursive=False):
    if recursive:
        out = f"rm -r {dir}"
    else:
        out = f"rm {dir}"
    print(out)
    subprocess.call(out, shell=True)

# - run -----------------------------------------------------------------------
if __name__ == '__main__':
    # read in run information
    metadata = pd.read_csv("samples.csv")
    # combine all the separated FASTQ files into a single file for each sample
    for i, df in metadata.iterrows():
        print(f"\nSample {i+1} of {len(metadata)}")
        combine(df.Sample_ID)
    # Filter FASTQ files before mapping
    fastq = glob.glob("*.fastq.gz")
    mkdir("filtered")
    for f in fastq:
        filter(f)
    # Map with STAR
    print("\n----- Mapping with STAR -----\n")
    gunzip("data/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz")
    gunzip("data/Mus_musculus.GRCm38.102.gtf.gz")
    print(f"STAR version {grab_version('star')}")
    generate_genome("data/Mus_musculus.GRCm38.dna.primary_assembly.fa",
                    "data/Mus_musculus.GRCm38.102.gtf")
    # map reads to genome index and count by gene
    mkdir('data/STAR')
    fastq = glob.glob("filtered/*.fastq.gz")
    for f in fastq:
        count(f, "data/Mus_musculus.GRCm38.102.gtf")
    # - Remove unecessary files --
    # first move count files to `counts`
    counts = glob.glob(f"data/STAR/*ReadsPerGene*")
    mkdir("counts")
    for count in counts:
        mv(count, "counts")
    # make tx2gene file for mapping gene_id and gene_name
    tx2gene("data/Mus_musculus.GRCm38.102.gtf")
    mv("data/tx2gene.csv", ".")
    # delete other unnecessary files
    print("\n----- Deleting unnecessary files -----\n")
    rm("data", recursive=True)
    rm("filtered", recursive=True)
