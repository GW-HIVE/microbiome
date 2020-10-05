import time
from Bio import Entrez

Entrez.email = 'lhopson15@gwu.edu'#input your email

in_file = '/Users/lindsayhopson/Documents/mouse_CensuScope_outputs/all_Acc_in_mice_samples1.3.txt' #use your appropriate file path

with open('names_all_orgs_lineage.txt', 'w') as output_File:

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
                        


                        output_File.write(str( i  +c + '\n'))
                        count = count + 1
                        print (i)
                        time.sleep(1)



