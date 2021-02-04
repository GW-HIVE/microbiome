#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
                            ##Unique Accession Identifier##
""""


    Author(s); Lindsay Hopson and Charles Hadley King
    Date: 02-04-2021
    Title of program/source code: Unique Accession Identifier
    Code version:2.1
    Type: Python 3.8.0




    Title: Accession Search
    Author:<Charles Hadley King>
    Date:<August 27, 2019>
    Code version:<1.0>
    Availibility:<https://github.com/GW-HIVE/microbiome>

    !!!IMORTANT: BEFORE RUNNING, COPE THIS FILE AND PASTE IT TO THE SAME FOLDER YOU SPECIFY ON LINE 39!!!
    
    This code will go through all your CensuScope output files, grab all unique accessions from the files, check to see if any of the unique accessions are in the Blacklist
    (you can read up on what the blacklist is in Hadley's paper he sent you and you can download the Blacklist from the GFBK website), the code removes all the accession that are in the blacklist, then outputs a list of accessions that make your new reference base that you will use for your Hexagon computations.
    Code 2.2 takes the list of accession created from Code 1.3 and grabs  their Organism name, taxonomy, ect from NCBI using Biopython. 
"""
################################################################################

#Declarations
import os, re, csv
my_dir = os.getcwd()
term = "dnaAccessionBasedResult" #identifier of all CensuScope output files
thresh = 10
fileList = []
list_Census = []
#blacklist = '/Users/username/Documents/folder/blackList-v2.0.csv' #use your appropriate file path and un-comment out line

#______________________________________________________________________________#
def createFileList(my_dir, term):
##GOAL of method: Creates a list of file paths to each of the CensuScope outputs. Each file is identified by the ending term of the file name
   
    for subdir, dirs, files in os.walk(my_dir):
        for file in files:
            if re.search(re.escape(term), file):
                result = os.path.join(subdir, file)
                fileList.append(result)
    return fileList
#______________________________________________________________________________#
def readCensusresultFile (file, list_Census, blacklist, threshold):
##     GOAL of method: Takes an file path, empty list, the blackList as a csv, and threshold as an input and returns a list with the accession numbers that are unique to the file
##        Accession numbers from the CensuScope output files that are below the threshold are ignored
##        Accession numbers from the CensuScope output files that meet the threshold are compared to the accession numbers in Black_List
##        The accessions that are not present in the Black_List are added to a new list called list_Census


    
    sample = file.split('_')[-1].split('.')[0] #not really sure what this does but its neccesary for the code to work
                                                #the has something to so with identifying the file name
  
    with open (file, 'r') as read:
        reader = csv.reader(read)
        for cell in reader:
            if cell[0] == 'accession': continue
            if cell[0] == 'n/a': continue
            accession = cell[0].split('.')[0].split('*')[-1]
            if accession in blacklist: continue
            if int(cell[5]) < threshold: continue
            if accession in list_Census: continue
            else: list_Census.append(accession)

    return list_Census

#______________________________________________________________________________#
def main ():
    
    fileList = createFileList(my_dir, term)
    for i in fileList:
        sample_list = readCensusresultFile(i, list_Census, blacklist, thresh)


#Outputs a file containing a list of accessions that were identified from all the mouse samples and that were above the 10 min match count threshold
#Output file is a .txt file named "'all_Acc_in_mice_samples1.3'." 
#Output file location will be within the same folder on your hard drive that this code is run from.
    with open('all_Acc_in_mice_samples1.3.txt', 'w') as file:
        for j in sample_list:
            print (j)
            file.write(j+'\n')
#______________________________________________________________________________#


if __name__ == '__main__': main()

















