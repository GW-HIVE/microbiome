#!/usr/bin/env python
# -*- coding: utf-8 -*-
#####################################################################################
                ##Unique Accession Organism Information Retrieval##


"""

        Author(s); Lindsay Hopson and Charles Hadley King
        Date: 02-09-2021
        Title of program/source code: Unique Accession Organism Information Retrieval
        Code version: 2.1
        Type: Python 3.8.0


        Title: Accession to TaxID
        Author:<Charles Hadley King>
        Date:<August 27, 2019>
        Code version:<1.0>
        Availibility:<https://github.com/GW-HIVE/microbiome>
        
        !!!IMORTANT: BEFORE RUNNING, COPY THIS SCRIPT FILE AND PASTE IT TO THE SAME FOLDER YOU SPECIFY ON LINE 39!!!

        Takes a list of accessions and returnes the Accession, Organsism name, and TaxID

        This code is able to pull information from NCBI using 
        Biopython package (http://biopython.org/DIST/docs/tutorial/Tutorial.pdf). In order
        to utilize this package in Python, follow the intallation instructions provided in
        the README link: https://github.com/biopython/biopython/blob/master/README.rst
        
        !!Be sure that Biopython is installed and up to date.        
        Use <pip3 install biopython> in command prompt for best results.
        See (https://biopython.org/wiki/Download) for more details!!      
"""
######################################################################################

import time
from Bio import Entrez

Entrez.email = 'example@email.com'#input an email address

in_file = '/Users/username/directory/folder/filename.txt' #(Use your appropriate file path and
                                                          #the name of the file that was 
                                                          #generated after running CensuScope_code1.3)

with open('orgNames_newAccList2.2.txt', 'w') as output_File:

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
######################################################################################
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
"""
