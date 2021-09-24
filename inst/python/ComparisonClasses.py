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
This function gets a nested dictionary where the keys are the proteins and the values the lipid-class pairs.
The lipid-class pairs are keys that contain a single value, that is the number of times the combination was found
"""
def get_ClassNestedDict(file):
    dict_Proteins_ClassNested={}
    with open(file,'r') as tbl:
        next(tbl)
        for line in tbl:
            line=line.strip("\n")
            line=line.replace('"','')
            ElementList=line.split(",")

            #Get genes
            proteins=ElementList[3]
            proteins=proteins.split(";")
            #Metabolite Classes
            Lip1Class=ElementList[len(ElementList)-2]
            Lip2Class=ElementList[len(ElementList)-1]

            LipClassFwd="%s_%s" % (Lip1Class,Lip2Class)
            LipClassRev="%s_%s" % (Lip2Class,Lip1Class)

            for protein in proteins:
                if protein not in dict_Proteins_ClassNested.keys():
                    dict_Proteins_ClassNested[protein]={}
                    dict_Proteins_ClassNested[protein][LipClassFwd]=1
                elif LipClassFwd not in dict_Proteins_ClassNested[protein].keys() and LipClassRev not in dict_Proteins_ClassNested[protein].keys():
                    dict_Proteins_ClassNested[protein][LipClassFwd]=1
                elif LipClassFwd in dict_Proteins_ClassNested[protein].keys():
                    dict_Proteins_ClassNested[protein][LipClassFwd]+=1
                elif LipClassRev in dict_Proteins_ClassNested[protein].keys():
                    dict_Proteins_ClassNested[protein][LipClassFwd]+=1
    return(dict_Proteins_ClassNested)

"""
This function compares the lipid class-pairs per protein of two treatments
"""
def compare_CoNI_Classes(file1,file2,outfileName="ComparisonClasses.csv",groupName1="Treatment1",groupName2="Treatment2"):
    #Get classes per treatment
    TFNLipClassTr1=get_ClassNestedDict(file1)
    TFNLipClassTr2=get_ClassNestedDict(file2)

    with open(outfileName,"w") as outfile:
        header="EdgeFeature\tVertexClassPair\t%s\t%s\n" % (groupName1,groupName2)
        outfile.write(header)
        #Loop lipid pairs treatment 2
        for k in TFNLipClassTr1.keys():
            protein=k
            if protein in TFNLipClassTr2.keys():#Check if protein is shared
                for k2 in TFNLipClassTr1[protein].keys():
                    lipclass=k2
                    if lipclass in TFNLipClassTr2[protein].keys():
                        Tr1Count=TFNLipClassTr1[protein][lipclass]
                        Tr2Count=TFNLipClassTr2[protein][lipclass]
                        stringToprint="%s\t%s\t%d\t%d\n" % (protein,lipclass,Tr1Count,Tr2Count)
                        print(stringToprint)
                        outfile.write(stringToprint)
                    elif lipclass not in TFNLipClassTr2[protein].keys():
                        Tr1Count=TFNLipClassTr1[protein][lipclass]
                        Tr2Count="Pair Class Missing"
                        stringToprint="%s\t%s\t%d\t%s\n" % (protein,lipclass,Tr1Count,Tr2Count)
                        print(stringToprint)
                        outfile.write(stringToprint)
                for k3 in TFNLipClassTr2[protein].keys():
                    lipclass=k3
                    if lipclass not in TFNLipClassTr1[protein].keys():
                        Tr1Count="Pair Class Missing"
                        Tr2Count=TFNLipClassTr2[protein][lipclass]
                        stringToprint="%s\t%s\t%s\t%d\n" % (protein,lipclass,Tr1Count,Tr2Count)
                        print(stringToprint)
                        outfile.write(stringToprint)
###########################################################################################################
###########################################################################################################
###########################################################################################################

Arguments = sys.argv[1:]
file_Treat1 = Arguments[0]
file_Treat2 = Arguments[1]
outfileName = Arguments[2]
groupName1 = Arguments[3]
groupName2 = Arguments[4]

compare_CoNI_Classes(file_Treat1,file_Treat2,outfileName,groupName1,groupName2)
