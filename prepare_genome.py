#!/usr/bin/env python
"""
@ERASMUS MC
Author: Diego Montiel 
Department: Genetic Identification
Preparare genome DNA and BLAST
"""
import os
import subprocess

#makeblastdb -in Homo_sapiens.GRCh38.cds.all.fa -dbtype nucl -out BLASTDB_human_cds/human_cds.db
#blastn -task blastn-short -query sequence.fasta  -db genome.db -outfmt 7 -evalue 0.05 -perc_identity 90 -out blastn_human_cds.out                                                                              
def menu():
    ans=True
    while ans:
        print ("""
        1.Make a database with BLAST
        2.Use BLASTN
        3.Exit/Quit
        """)
        ans=raw_input("What would you like to do? ") 
        if ans=="1":             
            print("\n Creating a database") 
            genome = raw_input("\n Please provide Genome sequence: ")
            database_name = raw_input("\n Please provide name for the database: ")
            print make_database(genome,database_name)
        elif ans=="2":
            print("\n BLASTN") 
            query = raw_input("\n Please provide query in multi-fasta: ")
            database = raw_input("\n Please provide database: ")
            output = raw_input("\n Please provide output name: ")
            print make_blast(query,database,output)
        elif ans=="3":
            print("\n Goodbye") 
            ans = None
        elif ans !="":
            print("\n Not Valid Choice Try again")

def make_database(genome,database_name):
 
    #Check if the output file already exists, if not, execute needle    
    if not os.path.isfile(database_name):       
        if subprocess.check_call('makeblastdb -in ' +genome + ' -out ' +database_name+ ' -dbtype nucl',\
            shell=True) == 0:                        
            return 'Success! database generated'
        else: 
            return 'Something went wrong, check parameters..'

def make_blast(query,database,output):
    #Check if the output file already exists, if not, execute needle    
    if not os.path.isfile(database):       
        if subprocess.check_call('blastn -query '+query + ' -db ' +database+ ' -out '+ output+\
            ' -task blastn-short -outfmt 7', shell=True) == 0:                        
            return 'Success! blast generated'
        else: 
            return 'Something went wrong, check parameters..'

if __name__ == "__main__":

    menu()
    
    
    
    
    