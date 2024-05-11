import gffutils
import os, sys
from Bio import SeqIO

# ______________GENE ANNOTATION FUNCTIONS______________
class GeneAnnotation:
    def __init__(self, gtf_file_path, ref_fasta_file_path):
        self.gtf_db_file_path = gtf_file_path
        self.ref_fasta_file_path = ref_fasta_file_path
        self.index = {}
    
    
    def create_index(self, gtf_file_path, ref_fasta_file_path, gtf_db_file_path = 'temp.db'):
        db = gffutils.create_db(gtf_file_path, dbfn=self.gtb_db_file_path, force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
        record = SeqIO.read(self.ref_fasta_file_path, "fasta")
        # get the index number location, go to the transcript, get the sequence, store it in a dict 
        # dict = {gene strand -> gene id}
        


    
    def find_gene(self, gene_name):
        db = gffutils.FeatureDB(self.gtf_db_file_path, keep_order=True)
        return db[gene_name]


# ______________PSEUDOALIGNMENT FUNCTIONS______________
class Pseudoalignment:
    
    def __init__(self, fasta_file_path):
        self.fasta_file_path = fasta_file_path
        self.debrujin = {'forward':{}, 'reverse':{}}
        self.equiv_classes = {}
    
    def construct_de_brujin(self, kmer_len, ga_instance):
        # Use Biopython's parse function to process individual
        # FASTA records (thus reducing memory footprint)
        for record in SeqIO.parse(self.fasta_file_path, 'fasta'):

            # Extract individual parts of the FASTA record
            identifier = record.id
            # description = record.description
            sequence = record.seq
            reverse_seq = sequence.reverse_complement()
            
            # get all the k_mers
            for i in range(0, len(sequence) - kmer_len - 1):
                kmer = sequence[i:i+kmer_len]
                kmer_next = sequence[i + 1:i + 1+kmer_len]
                r_kmer = reverse_seq[i:i+kmer_len]
                r_kmer_next = reverse_seq[i + 1:i + 1+kmer_len]
                if (kmer in self.debrujin['forward']):
                    self.debrujin['forward'][kmer]['ids'].add(identifier)
                    self.debrujin['forward'][kmer]['next'].add(kmer_next)
                elif kmer in self.debrujin['reverse']:
                    self.debrujin['reverse'][kmer]['ids'].add(identifier)
                    self.debrujin['reverse'][kmer]['next'].add(kmer_next)
                else:
                    self.debrujin['forward'][kmer]['ids'] = {identifier}
                    self.debrujin['reverse'][kmer]['ids'] = {identifier}
                    self.debrujin['forward'][kmer]['next'] = {kmer_next}
                    self.debrujin['reverse'][kmer]['next'] = {kmer_next}
                    
                if r_kmer in self.debrujin['forward']:
                    self.debrujin['forward'][r_kmer]['ids'].add(identifier)
                    self.debrujin['forward'][r_kmer]['next'].add(r_kmer_next)
                elif r_kmer in self.debrujin['reverse']:
                    self.debrujin['reverse'][r_kmer]['ids'].add(identifier)
                    self.debrujin['reverse'][r_kmer]['next'].add(r_kmer_next)
                else:
                    self.debrujin['forward'][r_kmer]['ids'] = {identifier}
                    self.debrujin['reverse'][r_kmer]['ids'] = {identifier}
                    self.debrujin['forward'][r_kmer]['next'] = {r_kmer_next}
                    self.debrujin['reverse'][r_kmer]['next'] = {r_kmer_next}

         # now take the union of all the paths in order to get the equivalence classes using the index
        for gene in ga_instance.index:
            pass