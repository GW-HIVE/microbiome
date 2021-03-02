#!/usr/bin/env python
# -*- coding: utf-8 -*-

#______________________________________________________________________________#
                            # Curated Accession List Generator #
""""
    Title: Curated Accession List
    Authors: <Charles Hadley King, Lindsay Hopson, and James Ziegler>
    Date: <02-11-2021>
    Code version: <1.0>
    Type: <Python 3.9.1>
    Availability:<https://github.com/GW-HIVE/microbiome>

    README:

    This code will go through all your CensuScope output files that include "dnaAccessionBasedResult" in the filename.
    Then, grab all unique accession number from the files, cross checking the numbers with the Blacklist to clean and remove any blacklisted accession numbers.
    Then, take the list of clean accession numbers and pull their organism name, taxonomy, ect. from NCBI using Biopython.
    Puts this list of accession numbers and relevant organism information into a .txt file that is ready for manual curation.

    This code is able to pull information from NCBI using 
    Biopython package (http://biopython.org/DIST/docs/tutorial/Tutorial.pdf). In order
    to utilize this package in Python, follow the intallation instructions provided in
    the README link: https://github.com/biopython/biopython/blob/master/README.rst
        
    
"""
#______________________________________________________________________________#

# Declarations
import os, re, csv
from Bio import Entrez
import time
my_dir = os.getcwd()
term = "dnaAccessionBasedResult" #identifier of all CensuScope output files
thresh = 10
fileList = []
list_Census = []

#!!!
blacklist = '/Users/username/folder/blackList-v2.0.csv'             # CHANGE ME: use your appropriate file path
Entrez.email = 'example@email.com'                                                  # CHANGE ME: input an email address
in_file = '/Users/username/folder/DELETE_ME_clean_acclist.txt' # CHANGE ME: Use same file path as for "blacklist"                                                     
#!!!

#______________________________________________________________________________#

def createFileList(my_dir, term):
# GOAL of method: Creates a list of file paths to each of the CensuScope outputs. Each file is identified by the ending term of the file name
   
    for subdir, dirs, files in os.walk(my_dir):
        for file in files:
            if re.search(re.escape(term), file):
                result = os.path.join(subdir, file)
                fileList.append(result)
    return fileList
#______________________________________________________________________________#

def readCensusresultFile (file, list_Census, blacklist, threshold):
# GOAL of method: Takes an file path, empty list, the blackList as a csv, and threshold as an input and returns a list with the accession numbers that are unique to the file
# Accession numbers from the CensuScope output files that are below the threshold are ignored
# Accession numbers from the CensuScope output files that meet the threshold are compared to the accession numbers in Black_List
# The accessions that are not present in the Black_List are added to a new list called list_Census
    
    sample = file.split('_')[-1].split('.')[0] #used to remove accession version numbers in order to pull the most up-to-date information
                                                                             
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
# GOAL of method: Outputs a file containing a list of accessions that were identified from all the mouse samples and that were above the 10 min match count threshold    
    
    fileList = createFileList(my_dir, term)
    for i in fileList:
        sample_list = readCensusresultFile(i, list_Census, blacklist, thresh)

    with open('DELETE_ME_clean_acclist.txt', 'w') as file: # in your file browser, delete this file once you have finished manual curation steps.
        for j in sample_list:
            print (j)
            file.write(j+'\n')
#______________________________________________________________________________#
# GOAL of method: create a .txt file with clean accession numbers

if __name__ == '__main__': main()

#______________________________________________________________________________#
# GOAL of method: Sleep for n seconds allowing the list with clean accession numbers to populate.

time.sleep(5)

#______________________________________________________________________________#
# GOAL of method: pull organism information from NCBI for each accession number

with open('curated_accession_list.txt', 'w') as output_File:

        with open(in_file, 'r') as input_File:
                id_list = input_File.readlines()
                count = 0
                for i in id_list:
                        handle = Entrez.esummary(db ="nucleotide", id=i)
                        record = Entrez.read(handle)
                        a = str(record[0]['Title'])
                        b = str(record[0]['TaxId'])
                        a = a.replace(',',' ')
                        handle = Entrez.efetch( db ="taxonomy", id=b)
                        record = Entrez.read(handle)
                        c = record[0]['Lineage']

                        output_File.write(str('\n'+i + '\t' + a + '\t'+ b+'\n'))
                        count = count + 1
                        print (i+'***'+a+'***'+ b+'***'+c)
                        time.sleep(1)
#______________________________________________________________________________#
                        
"""
Troubleshooting:

-If given the error: "urllib.error.URLError: <urlopen error 
[SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed: 
unable to get local issuer certificate (_ssl.c:1123)>" 
          
          And the output file is empty (0KB), 
          then complete the following steps and then rerun code:
          If you’re using macOS go to Macintosh HD > Applications > Python3.6 folder 
          (or whatever version of python you’re using) 
          > double click on “Install Certificates.command” file.
          
          See (https://stackoverflow.com/questions/50236117/scraping-ssl-certificate-verify-failed-error-for-http-en-wikipedia-org)
          for more details.
          
-If given the error: "ModuleNotFoundError: No module named 'Bio'"
          Be sure that Biopython is installed and up to date.        
          Use <pip3 install biopython> in command prompt for best results.
          See (https://biopython.org/wiki/Download) for more details!!
"""
