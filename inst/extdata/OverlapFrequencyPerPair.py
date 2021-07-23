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

##################################################Functions##################################################
"""
This function makes shure a file exists
"""
def does_file_exist(name_file):
    try:
        open(name_file,'r')
    except IOError:
        print("The file %s can't be found" % name_file)
        sys.exit(1)

"""
This function reads the folders in a given directory and searches for specific files according to the user input string.
"""
def get_Files(directory,StringFile):
    #Get folders in current directory
    folders=glob("%s/*/" % directory)
    files_toReadWithPath=[]
    #loop folders
    for f in folders:
        #get files in folder
        files=os.listdir(f)
        #use regular expression to get the correct file
        fTablNetwork=list(filter(lambda x: re.search(r'%s' % StringFile, x), files))[0] #filter returns iterator, call list to get list
        files_toReadWithPath.append("%s%s" % (f,fTablNetwork))# add the path to the name file
    return(files_toReadWithPath)

"""
This function checks if the files stored in a list can be read
"""
def check_files(listFiles):
    for f in listFiles:
        try:
            open(f,'r')
        except IOError:
            print("File %s does not exist\n" % (f))
            print("Provide correct path or make sure files are in directory\n")

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

##################################################end##################################################

Arguments = sys.argv[1:]

if len(Arguments)==1: #User should specify the folder where the CoNI output is located. For every CoNI run there is a folder.
    cwd = Arguments[0]
    filesPath=get_Files(cwd,"TableFor")
else:
    cwd = os.getcwd()
    #Get files that were used to generate CONI networks
    filesPath=get_Files(cwd,"TableFor")

#Read the files
#Create dictionary where results will be stored
dic_numberProteins={}
#Save metabolite pair keys in a list to loop later in the same order
#Might not be necessary as I think the new feature of Python is to keep the keys in the same order as they were added
keys_toloop=[]
#Loop files
for f in filesPath:
    with open(f,'r') as tbl:
        next(tbl) #skip header
        for line in tbl:
            line=line.strip("\n")
            line=line.replace('"','') #R output files have "" for all strings, get rid of them
            ElementList=line.split(",")
            dictkey_fw="%s/%s" % (ElementList[0],ElementList[1]) #metabolite pair forward #metabolite pair forward
            dictkey_rv="%s/%s" % (ElementList[1],ElementList[0]) #metabolite pair backwards
            proteins=ElementList[3]
            proteins=proteins.split(";")
            if dictkey_fw not in dic_numberProteins.keys():
                if dictkey_rv not in dic_numberProteins.keys():
                    dic_numberProteins[dictkey_fw]={}
                    keys_toloop.append(dictkey_fw)
                    #Sum gene ocurrences for metabolite pair
                    add_count_proteins(dic_numberProteins,proteins,dictkey_fw)
                else:
                    #Sum gene ocurrences for metabolite pair
                    add_count_proteins(dic_numberProteins,proteins,dictkey_rv)
            else:
                #Sum gene ocurrences for metabolite pair
                add_count_proteins(dic_numberProteins,proteins,dictkey_fw)

with open("OverlapFrequency_perPair.csv","w") as outfile:
    outfile.write("Lipid1\tLipid2\tprotein\tfrequency\n")
    for k in keys_toloop:
        proteins=k.split('/')
        p1=proteins[0]
        p2=proteins[1]
        for key in dic_numberProteins[k]:
            protein=key
            freq=dic_numberProteins[k][protein]
            string_toPrint="%s\t%s\t%s\t%d\n" % (p1,p2,protein,freq)
            outfile.write(string_toPrint)



#Get current directory
