#!/usr/bin/env python
# -*- coding: utf-8 -*-
#####################################################################################
                ##Unique Accession Organism Information Retrieval##


"""

        Author(s); Lindsay Hopson and Charles Hadley King
        Date: 02-05-2020
        Title of program/source code: Unique Accession Organism Information Retrieval
        Code version:2.0
        Type: Python 3.8.0


        Title: Accession to TaxID
        Author:<Charles Hadley King>
        Date:<August 27, 2019>
        Code version:<1.0>
        Availibility:<https://github.com/GW-HIVE/microbiome>


        Takes a list of accessions and returnes the Accession, Organsism name, and TaxID

        This code is able to pull information from NCBI using 
        Biopython package (http://biopython.org/DIST/docs/tutorial/Tutorial.pdf). In order
        to utilize this package in Python, follow the intallation instructions provided in
        the README link: https://github.com/biopython/biopython/blob/master/README.rst

"""
######################################################################################

import time
from Bio import Entrez

Entrez.email = 'lhopson15@gwu.edu'#input your email

in_file = '/Users/lindsayhopson/Documents/mouse_CensuScope_outputs/newAccList1.2.txt' #use your appropriate file path

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


