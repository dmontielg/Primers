#!/usr/bin/env python
"""
@ERASMUS MC
Author: Diego Montiel 
Department: Genetic Identification
PCR Identification and validation of primers in a Genome
"""
import sys
import itertools
import csv

def parse_blastfmt7(out_blast):
    """
    Parse a Blast output format 7 and return in a dictionary
    """
    dictionary = {} #Accession as key, sequence as value
    for entry in out_blast:
        if "#" not in entry:
            line = entry.split("\t")
            if line[0] in dictionary.keys():
                dictionary[line[0]].append(line[1:])
            else:
                dictionary[line[0]] = [line[1:]]
    return dictionary

def compare_primers(dictionary,cutoff_length):
    #cutoff_length = 5000
    dictionary_primers = []
    permut = list(itertools.permutations(dictionary,2))
    #Gives the possible combinations of the hists against the others hits
    #taking in account position seq[0] subjec tid    
    for iteration in permut:
        for seq1 in dictionary[iteration[0]]:
            for seq2 in dictionary[iteration[1]]:
                #If the query and seq are in the same position/chromosome/gene
                if seq1[0] == seq2[0]:
                    #If the they are in + or - strand
                    if int(seq1[8]) < int(seq2[7]):
                        size_seq1 = min(int(seq1[7]),int(seq1[8]))
                        size_seq2 = min(int(seq2[7]),int(seq2[8]))                        
                        length = size_seq2 - size_seq1    
                        #Check the length with the cutoff
                        if length <= cutoff_length:                                                        
                            #l.append((int(iteration[0]),int(iteration[1]))) 
                            l_tmp = (iteration[0],iteration[1])
                            dictionary_primers.append({"key":l_tmp,"seq1":seq1,\
                                "seq2":seq2,"length":length})                             
    return dictionary_primers
        
def write_report(dictionary_primers,fasta_file,output):
    
    #Sort the list of dictionaries by key 
    newlist = sorted(dictionary_primers, key=lambda k: k['key']) 
    with open(output, 'w') as csvfile:
        fieldnames = [  'Chromosome', 
                        'Primer_forward', 
                        'Primer_reverse', 
                        'Strand',                        
                        'Mismatches_Forward',
                        'Mismatches_Reverse',
                        'GC_content_Forward',
                        'GC_content_Reverse',                                             
                        'Location_start',
                        'Location_end',
                        'Length_KB',
                        'Evalue']
        
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for item in newlist:
            MIS_SEQ1 = str(int(len(fasta_file[str(item["key"][0])])) - int(item["seq1"][2]))
            MIS_SEQ2 = str(int(len(fasta_file[str(item["key"][1])])) - int(item["seq2"][2]))
            GC_SEQ1 = str(calculateGCContent(fasta_file[str(item["key"][0])]))
            GC_SEQ2 = str(calculateGCContent(fasta_file[str(item["key"][1])]))
            size_seq1 = min(int(item["seq1"][7]),int( item["seq1"][8]))
            size_seq2 = min(int(item["seq2"][7]),int( item["seq2"][8]))    
            
            #No more than 2 mismatches per sequence
            if int(MIS_SEQ1) < 3 and int(MIS_SEQ2) < 3:
                if size_seq1 < size_seq2:                    
                    location_start = size_seq1
                    location_end = int(size_seq2)-1
                else:
                    location_start = size_seq2
                    location_end = int(size_seq1)-1  
                    
                strand_seq1 = str(get_strand(size_seq1,size_seq2))
                strand_seq2 = str(get_strand(size_seq2,size_seq1))                
                writer.writerow({ 'Chromosome': item["seq1"][0] ,
                                'Primer_forward': item["key"][0] , 
                                'Primer_reverse': item["key"][1] , 
                                'Strand':strand_seq1+"/"+strand_seq2,                              
                                'Mismatches_Forward': MIS_SEQ1,
                                'Mismatches_Reverse': MIS_SEQ2,
                                'GC_content_Forward':GC_SEQ1+"%",
                                'GC_content_Reverse':GC_SEQ2+"%",
                                'Location_start': group(int(location_start)),
                                'Location_end': group(int(location_end)),
                                'Length_KB':item["length"] ,
                                'Evalue':item["seq1"][9]+" - "+item["seq2"][9] })            
    
    print "Report succesfully generated!"
    return True

def group(number):
    s = '%d' % number
    groups = []
    while s and s[-1].isdigit():
        groups.append(s[-3:])
        s = s[:-3]
    return s + ','.join(reversed(groups))

def get_strand(pos1,pos2):
     
    if pos1 < pos2:
        return "+"
    else:
        return "-"

def read_fasta_file(fastafile):
    "Get sequences from fasta file and store in dictionary"
    ###My idea:
    infile = open(fastafile).read()
    entries = infile.split(">")[1:]
    #print entries
    fastadict = {} #Accession as key, sequence as value
    for entry in entries:
        accession = entry.partition("\n")[0].replace("\r","")
        #print accession
        sequence = entry.partition("\n")[2].partition("\n")[0]
        #print sequence
        fastadict[accession] = sequence.rstrip()
    return fastadict

def calculateGCContent(sequence):    
    g = sequence.count('G')
    c = sequence.count('C')
    gc = float(g + c)
    gc_content = float((gc/len(sequence))*100)
    return round(gc_content,2)
    
if __name__ == "__main__":
    "Run this file from here as a script"
    #Check if parameters are provided; if not, exit with explanation
    if len(sys.argv) < 5:
        print "Please provide three parameters Blast output file output \
        file name and cutoff"
        sys.exit()
    else:

        query = sys.argv[1] #Primers
        output_blast = sys.argv[2] #Blastn output
        output_file = sys.argv[3] #output filename
        cutoff_length = int(sys.argv[4]) #a cufoff in kb
        fasta_file = read_fasta_file(query)
        #Separate file in lines
        lines = [line.rstrip('\n') for line in open(output_blast)]
        dictionary = parse_blastfmt7(lines)
        dictionary_primers = compare_primers(dictionary,cutoff_length)
        write_report(dictionary_primers,fasta_file, output_file)

        