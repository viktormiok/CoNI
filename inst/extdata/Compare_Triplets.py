#Regular expressions
import re

import random
import sys
import subprocess
from itertools import chain
#To search for files to read
from glob import glob
#Get the path
import os

###########################################################################################################
#################################################Functions#################################################
###########################################################################################################
"""
The input of this function are the tables used for network contrusction in CoNI.
The input are two tables, one for each treatment.
The output are the trios that are shared between treatments (lipid pair and protein affecting relationship)
"""

def compare_CoNI_results(file1,file2,outfileName="Shared_Metabolites_and_Genes.csv"):
    #Three dictionaries for shared and unique pairs
    dic_numberProteins_tr1={} #unique for treatment 1
    dic_numberProteins_tr2={} #unique for treatment 2
    dic_Proteins_Shared={} #shared between treatments

    #Get Metabolite Pairs and their gene content
    treat1=get_LipDict(file1) #treatment 1 chow_210_Dict
    treat2=get_LipDict(file2) #treatment 2 chow_17_Dict

    #Standarise key direction
    treat1=compare_keysDict(treat2,treat1)

    #L to store shared metabolite pairs
    lip_k_shared=[]

    #Open outfile to save shared edges with same gene content
    with open(outfileName,"w") as outfile:
        outfile.write("Vertex_1\tVertex_2\tEdge\n")

        #Get shared metabolites, and edges with same gene content
        for k, v in treat2.items():
            #Get metabolites
            lipids=k.split('/')
            lip1=lipids[0]
            lip2=lipids[1]
            #Check if metabolite pair is in treatment 1
            if k in treat1.keys():
                lip_k_shared.append(k) #Store shared lipid pair
                #Compare the gene content to find shared genes for the same lipid pair
                for k2 in treat2[k].keys():
                    protein=k2
                    if k2 in treat1[k].keys():
                        if k not in dic_Proteins_Shared.keys():
                            dic_Proteins_Shared[k]={}
                            dic_Proteins_Shared[k][k2]=1
                            string_toWrite="%s\t%s\t%s\n" % (lip1,lip2,protein)
                            outfile.write(string_toWrite)
                        else:
                            dic_Proteins_Shared[k][k2]=1
                            string_toWrite="%s\t%s\t%s\n" % (lip1,lip2,protein)
                            outfile.write(string_toWrite)
            else:
                #Store unique pairs for treatment 2
                dic_numberProteins_tr2[k]=v

        for k,v in treat1.items():
            if k not in treat2.keys():
                dic_numberProteins_tr1[k]=v
        #print(dic_Genes_Shared)
        return dic_Proteins_Shared

"""
This function gets the metabolite dictionary of a specific treatment. A metabolite pair is the key and the genes/proteins are the values
"""
def get_LipDict(file):
    dic_numberProteins={}
    with open(file,'r') as tbl:
        next(tbl)
        for line in tbl:
            line=line.strip("\n")
            line=line.replace('"','')
            ElementList=line.split(",")
            dictkey_fw="%s/%s" % (ElementList[0],ElementList[1]) #metabolite pair forward
            dictkey_rv="%s/%s" % (ElementList[1],ElementList[0]) #metabolite pair backwards
            proteins=ElementList[3]
            proteins=proteins.split(";")
            if dictkey_fw not in dic_numberProteins.keys():
                if dictkey_rv not in dic_numberProteins.keys():
                    dic_numberProteins[dictkey_fw]={}
                    #Create a function that takes as input the dictionary and the genes
                    add_count_proteins(dic_numberProteins,proteins,dictkey_fw)
                else:
                    add_count_proteins(dic_numberProteins,proteins,dictkey_rv)
            else:
                add_count_proteins(dic_numberProteins,proteins,dictkey_fw)
        return(dic_numberProteins)

"""
This function compares the keys of two dictionaries, testing the keys of the first one against the second one
It tests two versions, a foward key 'a_b' and a backward key 'b_a'. If the backward key is found in the second
dictionary it replaces it with the foward version. It returns the second dictionary.
"""

def compare_keysDict(dict1,dict2):
    for k, v in dict1.items():
        fwd_key=k
        rev_key="%s/%s" % (k.split("/")[1],k.split("/")[0])

        if fwd_key in dict2.keys():
            continue
        elif rev_key in dict2.keys(): #test if pair of metabolites is shared but name backwards
            print(fwd_key,rev_key)
            dict2[fwd_key] = dict2[rev_key] #keep consistent between treatments, keep fwd version
            del dict2[rev_key]#remove reverse version
    return(dict2)

"""
Function that takes as input a dictionary (where every key is a pair of metabolites), a list of genes, and
a metabolite pair key. It generates a dictionary per gene and assigns as value +1 (to count)
"""

def add_count_proteins(dic_numberProteins,proteins,dickey):
    for p in proteins:
        if p not in dic_numberProteins[dickey].keys():
            dic_numberProteins[dickey][p]=1
        else:
            dic_numberProteins[dickey][p]+=1
    return(dic_numberProteins)


###########################################################################################################
###########################################################################################################
###########################################################################################################

Arguments = sys.argv[1:]
file_Treat1 = Arguments[0]
file_Treat2 = Arguments[1]
outfileName = Arguments[2]

compare_CoNI_results(file_Treat1,file_Treat2,outfileName)
