# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 00:13:01 2019

@author: tipputa
"""

from Bio import SeqIO
import os, sys
import pandas as pd

def get_gb_info_as_buffer(gb_file, buffer):
    for record in SeqIO.parse(gb_file, "genbank"):
        for feature in record.features:
            if feature.type in ["CDS", "rRNA", "tRNA"]:
                locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                product = feature.qualifiers.get("product", ["Unknown"])[0]
                gene = feature.qualifiers.get("gene", [""])[0]
                pos_start = str(feature.location.start).replace(">", "").replace("<", "")
                pos_end = str(feature.location.end).replace(">", "").replace("<", "")
                strand = str(feature.location.strand)
                nucleotide = str(feature.location.extract(record).seq)
                translation = feature.qualifiers.get("translation", [""])[0]
                buffer.append((record.id, locus_tag, feature.type, gene, product, pos_start, pos_end, strand, nucleotide, translation))            

def get_gb_info_as_buffer_for_takenaka191023(gb_file, buffer):
    for record in SeqIO.parse(gb_file, "genbank"):
        for feature in record.features:
            if feature.type in ["CDS", "rRNA", "tRNA"]:
                # use protein_id instead of locus_tag
                locus_tag = feature.qualifiers.get("protein_id", ["locus_tag", ""])[0]
                product = feature.qualifiers.get("product", ["Unknown"])[0]
                gene = feature.qualifiers.get("gene", [""])[0]
                pos_start = str(feature.location.start).replace(">", "").replace("<", "")
                pos_end = str(feature.location.end).replace(">", "").replace("<", "")
                strand = str(feature.location.strand)
                nucleotide = str(feature.location.extract(record).seq)
                translation = feature.qualifiers.get("translation", [""])[0]
                buffer.append((record.id, locus_tag, feature.type, gene, product, pos_start, pos_end, strand, nucleotide, translation))            


def write_gb_info(gb_file, output_path):
    buffer = []
    get_gb_info_as_buffer(gb_file, buffer)
    df = pd.DataFrame(buffer)
    df.to_csv(output_path, sep="\t", index=None, header=None)    

def write_gb_info_for_takenaka191023(gb_file, output_path):
    buffer = []
    get_gb_info_as_buffer_for_takenaka191023(gb_file, buffer)
    df = pd.DataFrame(buffer)
    df.to_csv(output_path, sep="\t", index=None, header=None)    

                    
def create_merged_csv(gb_dir, output_path, tag):
    files = os.listdir(gb_dir)
    num_gb = len([file for file in files if file.endswith(".gb")  or file.endswith('.gbk') or file.endswith(".gbff")])
    LoopCounter = 1
    
    buffer = [("acc","locusTag", "feature", "gene", "product", "start", "end", "strand", "nuc", "prot")]
    for file in files:
        if file.endswith(".gb") or file.endswith('.gbk'): 
            print("Reading: " + file + "  "+ str(LoopCounter) + "/" + str(num_gb))
            LoopCounter = LoopCounter + 1
            get_gb_info_as_buffer(gb_dir + "/" + file, buffer)
        elif file.endswith(".gbff"):
            print("Reading (gbff): " + file + "  " + str(LoopCounter) + "/" + str(num_gb))
            LoopCounter = LoopCounter + 1
            if(tag == "takenaka191023"):
                get_gb_info_as_buffer_for_takenaka191023(gb_dir + "/" + file, buffer)
    df = pd.DataFrame(buffer)
    df.to_csv(output_path, sep="\t", index=None, header=None)    

def write_each_csv(gb_dir, output_dir, tag):
    files = os.listdir(gb_dir)
    num_gb = len([file for file in files if file.endswith(".gb")  or file.endswith('.gbk') or file.endswith(".gbff")])
    LoopCounter = 1
    
    for file in files:
        if file.endswith(".gb") or file.endswith('.gbk'): 
            print("Reading: " + file + "  "+ str(LoopCounter) + "/" + str(num_gb))
            LoopCounter = LoopCounter + 1
            write_gb_info(gb_dir + "/" + file, output_dir + "/" + file + ".tsv")
        elif file.endswith(".gbff"):
            print("Reading (gbff): " + file + "  " + str(LoopCounter) + "/" + str(num_gb))
            LoopCounter = LoopCounter + 1
            if(tag == "takenaka191023"):
                write_gb_info_for_takenaka191023(gb_dir + "/" + file, output_dir + "/" + file + ".tsv")
                

def write_header(file_path):
    buffer = [("acc","locusTag", "feature", "gene", "product", "start", "end", "strand", "nuc", "prot")]
    df = pd.DataFrame(buffer)
    df.to_csv(file_path, sep="\t", index=None, header=None)    


if __name__ == '__main__':
    gb_dir =  "C:/Users/tipputa/Desktop/takenaka191023/gb/"    
    output_dir =  "C:/Users/tipputa/Desktop/takenaka191023/tsv/"    
    output_summary_file =  "C:/Users/tipputa/Desktop/takenaka191023/summary_gb.tsv"    

    tag = "takenaka191023"
    # create_merged_csv(gb_dir, output_summary_file, tag)

    if len(sys.argv)==3:
        gb_dir = sys.argv[1] + "/"
        output_dir = sys.argv[2] + "/"
        
    write_header(output_dir + "header.txt")
    write_each_csv(gb_dir, output_dir, tag)
    # run cat command in the all files
    