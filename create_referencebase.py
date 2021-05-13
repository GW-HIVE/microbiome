#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
                            ##Create Referencebase##
""""


    Author(s); Lindsay Hopson and Charles Hadley King
    Date: 04-16-2021
    Title of program/source code: Create Referencebase
    Code version:3
    Type: Python 3.8.0
    Availibility:<https://github.com/GW-HIVE/microbiome>
    note: This version combined CensuScope_code1.3 and CensuScope_code2.2


    This code will go through all your CensuScope output files, grab all unique accessions from the files, check to
    see if any of the unique accessions are in the Blacklist (downloaded from https://hive.biochemistry.gwu.edu/gfkb),
    removes all the accession that are in the blacklist, then outputs a list of unique accessions (which will be used
    to create the reference base required for the Hexagon computations. Biopython package is then used to extract
    reference information from NCBI on the list of unique organisms. Non-bacteria organisms are first filtered out.
    The organims with incomplete genomes are filtered out and saved to a seperate CSV file for further manual curation.

    Output: This code will produce 3 CSV files. The first CSV file contains all unique organisms extracted from
    CensuScope output CSV files (include all organims (bacteria, viral, host, ect)). The second CSV file contains only
    bacteria with complete genomes (which is the filtered list). The third CSV file contains bacteria that do not have
    complete genomes. These bacteria will have to be manually curated to select its most representative complete genome.
"""
################################################################################
#Imports
import os, re, csv
import time
from Bio import Entrez

#Declarations
my_dir = os.getcwd()
term = "dnaAccessionBasedResult" #identifier of all CensuScope output files
gutDB = '/Users/lindsayhopson/Downloads/BMM/GutFeelingKnowledgeBase-v3.csv' #use your appropriate file path
thresh = 10
fileList = []
list_Census = []
blacklist = '/Users/lindsayhopson/Downloads/BMM/blackList-v2.0.csv' #use your appropriate file path
Entrez.email = 'lhopson15@gwu.edu'#input your email
bacteria_completeGenome = []
not_bacteria = []
bacteriaNot_completeGenome = []
#______________________________________________________________________________#

#GOAL of method: Creates a list of file paths to each of the CensuScope outputs. Each file is identified by the ending term of the file name

def createFileList(my_dir, term):
   
    for subdir, dirs, files in os.walk(my_dir):
        for file in files:
            if re.search(re.escape(term), file):
                result = os.path.join(subdir, file)
                fileList.append(result)
    return fileList
#______________________________________________________________________________#

#GOAL of method: Takes an file path, empty list, the blackList as a csv, and threshold as an input and returns a list with the accession numbers that are unique to the file
#Accession numbers from the CensuScope output files that are below the threshold are ignored
#Accession numbers from the CensuScope output files that meet the threshold are compared to the accession numbers in Black_List
#The accessions that are not present in the Black_List are added to a new list called list_Census

def readCensusresultFile (file, list_Census, blacklist, threshold):
    
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

#GOAL of method: Takes the unique list of organisms pulled from all the CensuScope output file and extracts their reference info from NCBI.
#The reference info is then put into a csv file.

def ncbi_ref_info(orgList):

    with open('unfiltered_unique_orgs.csv', 'w') as output_File:
        for i in orgList:
            handle = Entrez.esummary(db ="nucleotide", id=i)
            record = Entrez.read(handle)
            a = str(record[0]['Title'])
            b = str(record[0]['TaxId'])
            a = a.replace(',',' ')
            handle = Entrez.efetch( db ="taxonomy", id=b)
            record = Entrez.read(handle)
            c = record[0]['Lineage']


            output_File.write(str(i + ","+ a + ','+ c +'\n'))

    return 'unfiltered_unique_orgs.csv'

#--------------------------------------------------------#

#Goal of method: Filter through the unique accessions list and remove non-bacteria and incomplete genomes. Save the filtered list as csv
#The filtering process puts all removed organims in a list.
#This was done to that, if needed, one could retrieve the excat organims removes with few minor changes to this method.

def completeBacteria(fileIn):
    with open ('complete_bacteria.csv', 'w') as output_File:
        with open (fileIn, 'r') as read:
            reader = csv.reader(read)
            for cell in reader:
                Genbank_Accession = cell[0]
                Lineage = cell[2].split('; ')
                Organism_Name = cell[1]
    
                if cell[2] == 'Lineage': continue
                if cell[1] == "Organism Name": continue
           
            
                if "Bacteria" in Lineage:
                    if "genome" in Organism_Name:
                        if "complete genome" in Organism_Name:
                            bacteria_completeGenome.append(Genbank_Accession)
                        if "draft genome" in Organism_Name:
                            bacteriaNot_completeGenome.append(Genbank_Accession)
                        else: bacteria_completeGenome.append(Genbank_Accession)
                    else: bacteriaNot_completeGenome.append(Genbank_Accession)
                else: not_bacteria.append(Genbank_Accession)

            for i in bacteria_completeGenome:
                handle = Entrez.esummary(db ="nucleotide", id=i)
                record = Entrez.read(handle)
                a = str(record[0]['Title'])
                b = str(record[0]['TaxId'])
                a = a.replace(',',' ')
                handle = Entrez.efetch( db ="taxonomy", id=b)
                record = Entrez.read(handle)
                c = record[0]['Lineage']


                output_File.write(str(i + ","+ a + ','+ c +'\n'))

    return bacteriaNot_completeGenome

#-------------------------------------------------------------------------#

# GOAL of method: Produces a csv file containing bacteria with incomplete genomes, which will have to be manually curated to find most represenatitve genome.

def manualQC(not_complete_genomes):
    with open('bacteria_incomplete_genomes.csv', 'w') as output_File: 
        for i in not_complete_genomes:
            handle = Entrez.esummary(db ="nucleotide", id=i)
            record = Entrez.read(handle)
            a = str(record[0]['Title'])
            b = str(record[0]['TaxId'])
            a = a.replace(',',' ')
            handle = Entrez.efetch( db ="taxonomy", id=b)
            record = Entrez.read(handle)
            c = record[0]['Lineage']

            output_File.write(str(i + ","+ a + ','+ c +'\n'))

    return bacteriaNot_completeGenome
        


#--------------------------------------------------------------------------#
def main ():
    
    fileList = createFileList(my_dir, term)
    for i in fileList:
        sample_list = readCensusresultFile(i, list_Census, blacklist, thresh)

    in_file = ncbi_ref_info(sample_list)
    
    needQC_list = completeBacteria(in_file)

    manualQC(needQC_list)

    


        
          
#______________________________________________________________________________#


if __name__ == '__main__': main()


















