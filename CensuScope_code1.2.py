#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
                            ##Unique Accession Identifier##
""""


    Author(s); Lindsay Hopson and Charles Hadley King
    Date: 02-05-2020
    Title of program/source code: Unique Accession Identifier
    Code version:2.0
    Type: Python 3.8.0




    Title: Accession Search
    Author:<Charles Hadley King>
    Date:<August 27, 2019>
    Code version:<1.0>
    Availibility:<https://github.com/GW-HIVE/microbiome>


    This code reads through CensuScope output files and returns a file listing unique organisms (organsims that are not present in GFKB or the blacklist)
"""
################################################################################

#Declarations
import os, re, csv
my_dir = os.getcwd()
term = "dnaAccessionBasedResult" #identifier of all CensuScope output files
gutDB = '/Users/lindsayhopson/Documents/mouse_CensuScope_outputs/GutFeelingKnowledgeBase-v3.csv' #use your appropriate file path
thresh = 10
fileList = []
list_Census = []
blacklist = '/Users/lindsayhopson/Documents/mouse_CensuScope_outputs/blackList-v2.0.csv' #use your appropriate file path

#______________________________________________________________________________#
def createFileList(my_dir, term):
##     GOAL of method: Creates a list of file paths to each of the CensuScope outputs
##        Each file is identified by the ending term of the file name

   
    for subdir, dirs, files in os.walk(my_dir):
        for file in files:
            if re.search(re.escape(term), file):
                result = os.path.join(subdir, file)
                fileList.append(result)
    return fileList
#______________________________________________________________________________#
def readCensusresultFile (file, list_Census, blacklist, threshold):
##     GOAL of method: Takes an file path, empty list, the blackList as a csv, and threshold as an input and returns a list with the accession numbers that are unique to the file
##        A list called Black_List is created that only contains accession numbers from the blacklist
##        Accession numbers from the CensuScope output files that are below the threshold are ignored
##        Accession numbers from the CensuScope output files that meet the threshold are compared to the accession numbers in Black_List
##        The accessions that are not present in the Black_List are added to a new list called list_Census



    Black_List = []
    with open (blacklist, 'r') as read:
        reader = csv.reader(read)
        for cell in reader:
            if cell[0] == 'Accession': continue
            blk_accession = cell[0].split('.')[0]
            Black_List.append(blk_accession)

    
    sample = file.split('_')[-1].split('.')[0] #not really sure what this does but its neccesary for the code to work
                                                #the has something to so with identifying the file name
  
    with open (file, 'r') as read:
        reader = csv.reader(read)
        for cell in reader:
            if cell[0] == 'accession': continue
            if cell[0] == 'n/a': continue
            accession = cell[0].split('.')[0].split('*')[-1]
            if accession in Black_List: continue
            if int(cell[5]) < threshold: continue
            if accession in list_Census: continue
            else: list_Census.append(accession)

    return list_Census


###______________________________________________________________________________#
def getGutDB (file):
##    GOAL of method: Reads through GFKB and puts only the accession numbers into a list called gutDBList

    gutDBList = []

    with open (file, encoding ='ISO-8859-1') as read: 
        reader = csv.reader(read)
        for cell in reader:
            if cell[6] == 'GenBank Accession': continue
            
            accession = cell[6].split('.')[0].split('*')[0]
            gutDBList.append(accession)

    return gutDBList

#______________________________________________________________________________#

def gfkb_compare (newList, gfkb):
##     Goal of method: checks the content of a newList (list of CensuScope output accessions that are not in the Blacklist) against GFKB accession numbers.
##     Unique items are added to a new list called new_accession_List


    new_accession_List = []
    for accession in newList:
        if accession in gfkb: continue
        else:
            new_accession_List.append(accession)

    return new_accession_List

#______________________________________________________________________________#
def main ():
    
    fileList = createFileList(my_dir, term)
    gutDBList = getGutDB(gutDB)
    for i in fileList:
        sample_list = readCensusresultFile(i, list_Census, blacklist, thresh)
    final_List = gfkb_compare(sample_list, gutDBList)


#the final accession list in printed to a txt file called "newAccList1.2.txt"
    
    with open('newAccList1.2.txt', 'w') as file:
        for i in final_List:
            print (i)
            file.write(i+'\n')
#______________________________________________________________________________#


if __name__ == '__main__': main()


















