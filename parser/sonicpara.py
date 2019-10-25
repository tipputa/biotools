# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 23:46:01 2019

Parser of SonicParanoid

@author: tipputa
"""

import pandas as pd
import os, sys

def read_ortholog_group(work_dir):
    ortholog_group_file = work_dir + "/ortholog_groups.tsv"
    return pd.read_csv(ortholog_group_file, sep="\t")

def get_single_copy_genes(array, numThreshold):
    """ cluster size == num genomes, num genomes >= threshold"""
    return array.iloc[2] == array.iloc[1] and array.iloc[2] >= numThreshold

def getMap(df, key, value):
    return df.set_index(key)[value].to_dict()

def string2list(s, tagList):
    if(s != "*"):
        tagList += s.split(",")
        
def get_protein_from_locusTag(dic, locusTagList):   
    proteinList = [">" + locusTag + "\n" + dic[locusTag] for locusTag in locusTagList]         
    return "\n".join(proteinList)

def convert_asterisk_to_hyphen(arr):
    return [str(i).replace("*", "-") for i in arr]


def apply_consensusID_proteinFasta(series, locusTag_protein_dic, protein_dir):
    """ modified 191024 to adopt original file format """
    cId = str(series.iloc[0])
    locusTagList = []
    series.iloc[1:].apply(string2list, tagList = locusTagList)
    proteinFasta = get_protein_from_locusTag(locusTag_protein_dic, locusTagList)
    with open(protein_dir + "/" + cId + ".fas", "w") as f:
        f.write(proteinFasta)    
    print("fin: " + cId)


def convert_ortholog_groups_to_parsed_tsv(work_dir, tag):
    # output name of parsed sonicParanoid
    all_locusTag_fName = work_dir + "/all_locusTag_id.tsv" 
    
    ortholog_group_df = read_ortholog_group(work_dir);
    numCol = len(ortholog_group_df.columns)
    numGenomes = int((numCol - 4) / 2)    
    print("Num genomes: " + str(numGenomes))

    # change col names as genbank file names;
    if(tag == "takenaka191023"):
        ortholog_group_df.columns = [i.replace(".protein.faa", "").replace(".gbff_protein.faa", "") for i in ortholog_group_df.columns]

    # select useful columns
    is_locusTag = [True]+ [False for i in range(3)] + [True if i % 2 == 0 else False for i in range(numGenomes*2) ]
    ortholog_group_tag = ortholog_group_df.iloc[:, is_locusTag]

    c = ortholog_group_tag.columns.to_series()
    c[0] = "ConsensusID"
    ortholog_group_tag.columns = c
    ortholog_group_tag.to_csv(all_locusTag_fName, sep="\t", index=None)

def create_multi_fasta(work_dir, protein_dir, gb_summary_file):
    all_locusTag_fName = work_dir + "/all_locusTag_id.tsv" 
    all_locusTags_df = pd.read_csv(all_locusTag_fName, sep="\t")
    gb_df = pd.read_csv(gb_summary_file, sep="\t", engine="python")
    locusTag_protein_dic = getMap(gb_df, "locusTag", "prot")
    all_locusTags_df.apply(apply_consensusID_proteinFasta, axis=1, locusTag_protein_dic = locusTag_protein_dic, protein_dir = protein_dir)    

if __name__ == '__main__':
    work_dir = "C:/Users/tipputa/Desktop/takenaka191023"
    protein_dir = ""
    gb_summary_file = ""
    tag = "takenaka191023"

    if len(sys.argv)==4:
        work_dir = sys.argv[1] + "/"
        protein_dir = sys.argv[2] + "/"
        gb_summary_file = sys.argv[3]
    
    print("start process: ortholog_groups.tsv")
    convert_ortholog_groups_to_parsed_tsv(work_dir, tag)
    print("fin: convert from ortholog_groups.tsv to all_locusTag_id.tsv")
    
    print("start making protein fas")
    create_multi_fasta(work_dir, protein_dir, gb_summary_file)

    print("fin.")
    