from Bio import SeqIO
from collections import defaultdict
from itertools import product
from tqdm import tqdm

class Kallisto:
    def __init__(self, index_file, reads_file, kmer_len, single_forgiving = False, debug = False):
        self.index_file = index_file
        self.reads_file = reads_file
        self.kmer_len = kmer_len
        self.single_forgiving = single_forgiving
        self.debug = debug

        self.debrujin = defaultdict(set)
        self.kmer_to_id = defaultdict(set)
        self.equivalence_counts = {}

    def generate_kmers(self, sequence):
        kmers = [sequence[i: i+self.kmer_len] for i in range(len(sequence) - self.kmer_len + 1)]
        return kmers

    def build_debrujin_graph(self):
        for record in tqdm(SeqIO.parse(self.index_file, 'fasta')):

            # Extract individual parts of the FASTA record
            identifier = record.id
            sequence = str(record.seq)
            reverse_seq = str(record.seq.reverse_complement())

            forward_kmers = self.generate_kmers(sequence)
            reverse_kmers = self.generate_kmers(reverse_seq)
            
            for kmer in forward_kmers + reverse_kmers:
                self.kmer_to_id[kmer].add(identifier)
                if kmer not in self.debrujin:
                    self.debrujin[kmer] = []
                prefix = kmer[:-1] #AC-
                suffix = kmer[1:] #-CT
                self.debrujin[prefix].add(suffix)
        print('Successfully created index.')
    
    def process_reads_file(self, early_stop=0, skipping=False):
        count = 0
        for record in tqdm(SeqIO.parse(self.reads_file, 'fasta')):
            identifier = record.id
            read = record.seq
            reverse_seq = read.reverse_complement()
            
            equiv_classes = self.align_read(read)
            r_equiv_classes = self.align_read(reverse_seq)
            
            if len(equiv_classes) > 0:
                equiv_classes = (str(equiv_classes)[1:-1]).replace('\'', '')
            else:
                equiv_classes = "NA"
            if len(r_equiv_classes) > 0:
                r_equiv_classes = (str(r_equiv_classes)[1:-1]).replace('\'', '')
            else:
                r_equiv_classes = "NA"

            if equiv_classes in self.equivalence_counts:
                self.equivalence_counts[equiv_classes] += 1
            else:
                self.equivalence_counts[equiv_classes] = 1
            if r_equiv_classes in self.equivalence_counts:
                self.equivalence_counts[r_equiv_classes] += 1
            else:
                self.equivalence_counts[r_equiv_classes] = 1
            if self.debug:
                print(f"Forward strand: {equiv_classes}\nReverse strand: {r_equiv_classes}")
                print(identifier)
            count += 1
            if count == early_stop:
                break
        print("")

    
    def align_read(self, read, skipping):
        """
        for getting a set of equivalence classes for a read 
        @param read (string)
        @return (set) set of equivalence classes
        """
        equivalence_classes = []
        kmers = self.generate_kmers(read)

        for kmer in kmers:
            possible_kmers = self.process_kmer(kmer)
            current_classes = set()
            #take the union to get all the possible equivalence classes
            for p_kmer in possible_kmers:
                prefix = p_kmer[:-1]
                suffix = p_kmer[1:]

                if prefix in self.debrujin and suffix in self.debrujin[prefix]:
                    if self.debug: 
                        print(f"Kmer: {p_kmer} Ids: {self.kmer_to_id[p_kmer]}")
                    current_classes.update(self.kmer_to_id[p_kmer])
            if not equivalence_classes:
                equivalence_classes = current_classes
            else:
                equivalence_classes &= current_classes
            
        #get intersection of all sets in equiv_classes array
        return equivalence_classes
    
    def process_kmer(self, kmer):
        """
        for all possibilities of kmers
        @param kmer (string): kmer string
        @returns 
            set of all kmer combinations
        """
        #kmer contains an N, potentially can have many possibilities, need to keep track of all paths
        kmers = set()
        if 'N' in kmer:
            positions = [i for i, base in enumerate(kmer) if base == 'N']
            nucleotides = ['A', 'T', 'G', 'C']
            for replacements in product('ATGC', repeat=len(positions)):
                new_kmer = list(kmer)
                for pos, replacement in zip(positions, replacements):
                    new_kmer[pos] = replacement
                kmers.add(''.join(new_kmer))
        else:
            kmers.add(kmer)
        if self.single_forgiving:
            nucleotides = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G'}
            for i in range(self.kmer_len):
                replaced_kmer = None
                if i == self.kmer_len - 1:
                    replaced_kmer = kmer[:i] + nucleotides[kmer[i]] + kmer[i+1:]
                else:
                    replaced_kmer = kmer[:i] + nucleotides[kmer[i]] + kmer[i+1:]
                kmers.add(replaced_kmer)
        return kmers
    