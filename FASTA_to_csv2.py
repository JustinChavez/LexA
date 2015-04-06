# -*- coding: utf-8 -*-
"""
Created on Tue Feb 24 14:25:49 2015

@author: nikhil
"""
from __future__ import division
import os
#from Bio import SeqIO
from Bio import Entrez
#import time
Entrez.email = "nikyesu1@umbc.edu"
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

#file4 = open("sample_fasta.fasta", "w") 
fasta_file = open("LexA_fasta.fasta")
#csv_file = open("lexA.csv" , "w")

#csv_file.write("gene name" + "," + "reference number" + "," + "blast gi" + "," + "organism" + "," + "protein sequence" + "," + "query sequence" + "," + "e value" + "," + "sequence length" + "," + "coverage" + "," + "blast score" + "," + "taxid" +"," +  "phylum" + "," + "class" + "," + "order" + "," + "family" + "," + "genus" + "\n")
                        

i=0
fasta_list = []
gi_list = []
flag = 0
count = 0
temp = ""
read = False
fasta_set = []
taxo_list = []
listp = []

def girefP(gi_data):
    print "parsing BLAST sequence info"
    filestring = gi_data
 
    par_list = []

    filestring = filestring +">$"
    for i in range(len(filestring)):
       
        
        if filestring[i] == ">":
            if filestring[i + 1] == "$":
                 print "end of entry \n"
            else:   
                
                count = 4 
                while(filestring[i+count] != "|"):
                    count += 1
            
                print filestring[i+4:i+ count]
                gi = filestring[i+4:i+ count]
                count2 = i+count + 5
                count3 = 0
            
                while(filestring[count2 + count3] != "|"):
                    count3 += 1
                print filestring[count2:count2 + count3]
                ref = filestring[count2:count2 + count3]
            
                count4 = count2 + count3 + 2
                count6 = 0
                
                try:                          
                
                    while (filestring[count4 + count6] != "[" ):
                        count6 += 1
                        
                except IndexError:
                    par_list = "Aberration in gene name"
                    return par_list
        
                print filestring[count4:count4 + count6]
                giName = filestring[count4:count4 + count6]               
                count7 = count4 + count6
                count5 = 0
                while (filestring[count7 + count5] != ">"):
                    count5 += 1
                print filestring[count7: count7 + count5]
                tempStr = filestring[count7: count7 + count5]
                tempStr = tempStr[1:]
                tempStr = tempStr[:-1]
                name = tempStr
                par_list.append([gi,ref,name, giName])
            print "blast sequence info parsed \n"
    return par_list
    
########################

def get_tax_id(species):
    """to get data from ncbi taxomomy, we need to have the taxid.  we can
    get that by passing the species name to esearch, which will return
    the tax id"""
    print "getting tax id"

    print species
    species = species.replace(" ", "+").strip()
    search = Entrez.esearch(term = species, db = "taxonomy", retmode = "xml")
    record = Entrez.read(search)

    
    
    try:
        taxid = record['IdList'][0]
    except IndexError:
        taxid = "No Tax ID found"

    print "tax id found \n"
    return taxid
    
    

def get_tax_data(taxmy):
    """once we have the tax """
    print "getting taxonomy info"
    print taxmy
    tax_list = []
    if taxmy == "No Tax ID found":
        tax_list == taxmy
        return tax_list
    
    #time.sleep(3)
    search = Entrez.read(Entrez.efetch(id = taxmy, db = "taxonomy", retmode = "xml"))[0]

    lineages = search['LineageEx']
    
    for lineage in lineages:
        i = 0
        
        print "getting lineage"
    
            
           
        if lineage['Rank'] == 'no rank':
            i += 1
            #print "no rank"            
      
        if lineage['Rank'] == 'superkingdom':
          i += 1
          #print "wrong rank"
                
        if lineage['Rank'] == 'subclass':
            i += 1
            #print "wrong rank"
                
        if lineage['Rank'] == 'suborder':
            i += 1
            #print "wrong rank"

        if lineage['Rank'] == 'phylum':
            phylum = lineage['ScientificName']           
            tax_list.append(phylum)
            i += 1
            
        if lineage['Rank'] == 'class':
            sclass = lineage['ScientificName']           
            tax_list.append(sclass)
            i += 1
            
        if lineage['Rank'] == 'order':
            order = lineage['ScientificName']
            tax_list.append(order)
            i += 1
            
        if lineage['Rank'] == 'family':
            family = lineage['ScientificName']
            tax_list.append(family)
            i += 1
            
        if lineage['Rank'] == 'genus':
            genus = lineage['ScientificName']
            tax_list.append(genus)
            i += 1
        
        else:
            i += 1
             
            for x in (range(len(tax_list))):
                
                print tax_list[x]
            
                  
    return tax_list
                
        
        
################################################################


#iterate through each line of the file
for line in fasta_file:
    print "reading fasta and assigning GIs \n"
#indicates the start of a new paragraph in the file
    if line[0:3] == ">gi":
#ensures we don't append when there is nothing there
        if read == True:
            flag == 0
            fasta_list.append(temp)
            count += 1
        temp = line
        read = True
        flag = 1
#gets the length of the gene id
        k = 4
        while line[k] != "|":
            k += 1
        gi_list.append(line[4:k])
#this indicates we are still reading the paragraph into the variable
    if flag == 1:
#makes sure we do not add the frist line twice
        if line[0] != ">":
            temp += line
fasta_list.append(temp)
print "fasta read and assigned to GIs"
        
#print "there are: " + str(count) + "fasta in the list"
#print gi_list

#########################################
#blast each GI entry
for gi in gi_list:
    if not os.path.exists(gi + 'blast.out'):
        result_handle1 = NCBIWWW.qblast("blastp", "nr", gi)
        save_file = open(gi +"blast.out", "w")
        save_file.write(result_handle1.read())
        save_file.close()
        result_handle1.close()
##########################################


#set E value limit and count
E_VALUE_THRESH = (10 ** -10)
count = 0
seq_hits = []

#open parse function
def blastparse(gi_list):
    count1 = 0
    count2 = 0
      
#for each gi, find blast output file and parse
    for gi in gi_list:
        count2 += 1
        print gi 
        if gi == "226366183":
            
            csv_file = open("lexA"+ gi +".csv" , "w")
            csv_file.write("gene name" + "," + "reference number" + "," + "blast gi" + "," + "organism" + "," + "protein sequence" + "," + "query sequence" + "," + "e value" + "," + "sequence length" + "," + "coverage" + "," + "blast score" + "," + "taxid" +"," +  "phylum" + "," + "class" + "," + "order" + "," + "family" + "," + "genus" + "\n")
        #print gi
        #print "GICOUNT:" + str(count1)
            count1 += 1
        #handle = Entrez.efetch(db="protein", id=gi, rettype="fasta", retmode="text")
            result_handle2 = open(gi +"blast.out")
            blast_record = NCBIXML.read(result_handle2)
        
#file_name = open(gi + "blast_eval.txt", "w")  
            print "opening file"



#for each blast alignment, filter based on e value and coverage   
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    print "filtering blast hits"
                    if (hsp.expect < E_VALUE_THRESH) and (((hsp.query_end - hsp.query_start)/len(hsp.sbjct)) > .9):
                
# assigns valuable blast info to variables
                    #print alignment.title + "\n"	
                        gi_data = alignment.title                    
                        prot_seq = hsp.sbjct
                        query_seq = hsp.query
                        length = alignment.length
                        e_value = hsp.expect
                        blast_score = hsp.score
                        coverage = (100*(hsp.query_end - hsp.query_start)/len(hsp.sbjct))
                    
#parse for blast_gi, organism name, ref_id, gene name
                    
                        par_list = girefP(gi_data)
                        
                        if par_list == "Aberration in gene name":
                            if "," in gi_data:
                                
                                gi_data = gi_data.replace(',', '')
                                
                                csv_file.write(gi_data + "," + " " + "," + " " + "," + " " + "," + prot_seq + "," + query_seq + "," + str(e_value) + "," + str(length) + "," + str(coverage) + "," + str(blast_score) + "\n")
                            
                        else:                        
                            for listp in par_list:
                                for p in listp:
                                    if "," in str(p):                                
                                        str(p).replace(',', '')
                                
                            
# assigns 
                                blast_gi = listp[0]                      
                                refnum = listp[1] 
                                organism = listp[2]
                                gene_name = listp[3]
                        
                        
                            
                                organism = organism[:-1]
                                
                                taxmy = str(get_tax_id(organism))
                                taxo_list = get_tax_data(taxmy)
                                
                                for x in (range(len(taxo_list))):
                                    print taxo_list[x]
                                             
                      
#write data to csv file
                                if "," in str(gene_name):                                
                                    gene_name = str(gene_name).replace(',', '')
                        
                                                 
                            
                        
                                csv_file.write(gene_name + "," + str(refnum) + "," + blast_gi + "," + organism + "," + prot_seq + "," + query_seq + "," + str(e_value) + "," + str(length) + "," + str(coverage) + "," + str(blast_score) + "," + str(taxmy))
                                for x in (range(len(taxo_list))):
                                    csv_file.write( "," + taxo_list[x])
                                csv_file.write("\n")
                    
                                 

blastparse(gi_list)


#gi|15988320|pdb|1JHH|A Chain A, Lexa S119a Mutant >gi|15988321|pdb|1JHH|B Chain B, Lexa S119a Mutant>

                    





