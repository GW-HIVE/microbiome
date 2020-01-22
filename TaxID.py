

#################################################################### ############
##Accession to TaxId##
"""Takes a list of accessions and returnes the TaxId and name"""
#################################################################### ############

import time
from Bio import Entrez

Entrez.email = 'tud13145@temple.edu'


with open(in_file) as file: 
	id_list = file.readlines()

id_list = [word.strip() for word in id_list] 

#print 'id_list loaded'
#print id_list

count = 0
for i in id_list:
	handle = Entrez.esummary( db ="nucleotide", id=i) record = Entrez.read(handle)
	a = str(record[0]['Title'])
	b = str(record[0]['TaxId'])
	a = a.replace(',','')
	handle = Entrez.efetch( db ="taxonomy", id=b) 
	record = Entrez.read(handle)
	c = record[0]['Lineage']
	#file.write(str(i + '\t' + a + '\t'+ b+'\n'))
	count = count + 1 
	print i+','+a+','+ b+','+c 
	time.sleep(1)
# # print "done"
