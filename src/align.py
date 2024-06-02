from Bio import SeqIO
from collections import defaultdict
import copy

class Naive:
    def __init__(self, index_file, reads_file, kmer_len):
        self.index_file = index_file
        self.reads_file = reads_file
        self.kmer_len = kmer_len

        self.debrujin = defaultdict(list)
        self.kmer_to_id = defaultdict(list)
        self.equivalence_counts = {}

    def generate_kmers(self, sequence):
        kmers = [sequence[i: i+self.kmer_len] for i in range(len(sequence) - self.kmer_len + 1)]
        return kmers

    def create_index(self):
        for record in SeqIO.parse(self.index_file, 'fasta'):

            # Extract individual parts of the FASTA record
            identifier = record.id
            sequence = str(record.seq)
            reverse_seq = str(record.seq.reverse_complement())

            forward_kmers = self.generate_kmers(sequence)
            reverse_kmers = self.generate_kmers(reverse_seq)
            
            # get all the k_mers
            for i in range(0, len(sequence) - self.kmer_len - 1):
                kmer = sequence[i:i+self.kmer_len]
                kmer_next = sequence[i + 1:i + 1+self.kmer_len]
                r_kmer = reverse_seq[i:i+self.kmer_len]
                r_kmer_next = reverse_seq[i + 1:i + 1+self.kmer_len]
                if kmer in self.debrujin['forward']:
                    self.debrujin['forward'][kmer]['ids'].add(identifier)
                    self.debrujin['forward'][kmer]['next'].add(kmer_next)
                elif kmer in self.debrujin['reverse']:
                    self.debrujin['reverse'][kmer]['ids'].add(identifier)
                    self.debrujin['reverse'][kmer]['next'].add(kmer_next)
                else:
                    self.debrujin['forward'][kmer] = {'ids':None, 'next': None}
                    self.debrujin['reverse'][kmer] = {'ids':None, 'next': None}
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
                    self.debrujin['forward'][r_kmer] = {'ids':None, 'next': None}
                    self.debrujin['reverse'][r_kmer] = {'ids':None, 'next': None}
                    self.debrujin['forward'][r_kmer]['ids'] = {identifier}
                    self.debrujin['reverse'][r_kmer]['ids'] = {identifier}
                    self.debrujin['forward'][r_kmer]['next'] = {r_kmer_next}
                    self.debrujin['reverse'][r_kmer]['next'] = {r_kmer_next}
        print('Successfully created index.')
    
    def process_reads_file(self, early_stop=0):
        count = 0
        for record in SeqIO.parse(self.reads_file, 'fasta'):
            read = record.seq
            reverse_seq = read.reverse_complement()
            equiv_classes = self.align_read(read, 0, 'forward')
            r_equiv_classes = self.align_read(reverse_seq, 0, 'reverse')
            equiv_classes = str(set.union(*equiv_classes))[1:-1]
            r_equiv_classes = str(set.union(*r_equiv_classes))[1:-1]
            if equiv_classes in self.equivalence_counts:
                self.equivalence_counts[equiv_classes] += 1
            else:
                self.equivalence_counts[equiv_classes] = 1
            if r_equiv_classes in self.equivalence_counts:
                self.equivalence_counts[r_equiv_classes] += 1
            else:
                self.equivalence_counts[r_equiv_classes] = 1
            if count % 50 == 0:
                print(f"Iteration: {count}")
            count += 1
            if count == early_stop:
                break

    
    def align_read(self, read, index, direction, path = [], equiv_classes=[]):
        """
        for getting a set of equivalence classes for a read (recursive)
        @param path (array): array of equivalence classes (a set)
        """
        #base case
        kmer = read[index:index + self.kmer_len]
        if index == len(read) - self.kmer_len:
            possible_kmers = self.process_kmer(kmer, direction)
            equiv_class = []
            for p_kmer in possible_kmers:
                curr_equiv_classes = copy.deepcopy(path)
                if p_kmer in self.debrujin[direction]:
                    curr_equiv_classes.append(self.debrujin[direction][p_kmer]['ids'])
                    equiv_class.append(set.intersection(*curr_equiv_classes))
            equiv_classes.append(set.union(*equiv_class))
        else:
            next_nucleotide = read[index+self.kmer_len]
            #need to get the first position of it. then subsequently check only the next values of the current kmer
            possible_kmers = self.process_kmer(kmer, direction)
            for p_kmer in possible_kmers:
                if p_kmer in self.debrujin[direction]:
                    next_p_kmer = p_kmer[1:] + next_nucleotide
                    if next_p_kmer in self.debrujin[direction]:
                        path.append(self.debrujin[direction][p_kmer]['ids'])
                        self.align_read(read, index + 1, direction, path, equiv_classes)
                        path.pop()
        return equiv_classes
            
        #get union of all sets in equiv_classes array
        # return set.union(*equiv_classes)
    
    def process_kmer(self, kmer, direction):
        """
        for filtering errors, either one is returned (true match), or multiple ('N' or single mismatch), 
        or 0 matches are found
        @param kmer (string): kmer string
        @param direction (string): 'forward' or 'reverse'
        @returns 
            array of kmers that are in the debrujin graph (recall there is only one of kmer in graph)
        """
        #kmer contains an N, potentially can have many possibilities, need to keep track of all paths
        if kmer in self.debrujin[direction]:
            return [kmer]
        kmers = []
        if 'N' in kmer and kmer != "N" * self.kmer_len:
            nucleotides = ['A', 'T', 'G', 'C']
            kmers = [kmer.replace('N', nucleotide) for nucleotide in nucleotides if kmer.replace('N', nucleotide) in self.debrujin[direction]]
        else:
            nucleotides = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G'}
            kmers = []
            for i in range(self.kmer_len):
                replaced_kmer = None
                if i == self.kmer_len - 1:
                    replaced_kmer = kmer[:i] + nucleotides[kmer[i]] + kmer[i+1:]
                else:
                    replaced_kmer = kmer[:i] + nucleotides[kmer[i]] + kmer[i+1:]
                if replaced_kmer in self.debrujin[direction]:
                    kmers.append(replaced_kmer)
        return kmers
    