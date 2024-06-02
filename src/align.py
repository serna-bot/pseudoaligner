from Bio import SeqIO

class Naive:
    def __init__(self, index_file, reads_file, kmer_len):
        self.index_file = index_file
        self.reads_file = reads_file
        self.kmer_len = kmer_len

        self.debrujin = None

    def create_index(self):
        for record in SeqIO.parse(self.index_file, 'fasta'):

            # Extract individual parts of the FASTA record
            identifier = record.id
            sequence = record.seq
            reverse_seq = sequence.reverse_complement()
            
            # get all the k_mers
            for i in range(0, len(sequence) - self.kmer_len - 1):
                kmer = sequence[i:i+self.kmer_len]
                kmer_next = sequence[i + 1:i + 1+self.kmer_len]
                r_kmer = reverse_seq[i:i+self.kmer_len]
                r_kmer_next = reverse_seq[i + 1:i + 1+self.kmer_len]
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
    
    def process_reads_file(self):
        for record in SeqIO.parse(self.reads_file, 'fasta'):
            read = record.seq
            equiv_classes = self.align_read(read)
    
    def align_read(self, read):
        reverse_seq = read.reverse_complement()
        equiv_classes = []

        #need to get the first position of it. then subsequently check only the next values of the current kmer
        kmer = read[:self.kmer_len]
        r_kmer = reverse_seq[:self.kmer_len]
        
        for i in range(1, len(read) - self.kmer_len - 1):
            next_kmer = read[i:i+self.kmer_len]
            next_r_kmer = reverse_seq[i:i+self.kmer_len]
            equiv_class = []

            f_positions = self.process_kmer(kmer)
            r_positions = self.process_kmer(r_kmer)

            for f_pos in f_positions: #needs a queue or something to process it in order
                if next_kmer in f_pos['next']:
                    equiv_class.append(f_pos['ids'])
            for r_pos in r_positions:
                if next_r_kmer in r_pos['next']:
                    equiv_class.append(r_pos['ids'])

            if not equiv_class: 
                return -1
            equiv_classes.append(equiv_class)
            kmer = read[i:i+self.kmer_len]
            r_kmer = reverse_seq[i:i+self.kmer_len]
        
        #get union of all sets in equiv_classes array
        return set.union(*equiv_classes)
    
    def process_kmer(self, kmer):
        #kmer contains an N, potentially can have many possibilities, need to keep track of all paths
        if self.find_kmer_in_graph(kmer):
            return [kmer]
        elif 'N' in kmer and kmer != "N" * self.kmer_len:
            nucleotides = ['A', 'T', 'G', 'C']
            kmers = [kmer.replace('N', nucleotide) for nucleotide in nucleotides if self.find_kmer_in_graph(kmer.replace('N', nucleotide))]
            if kmers:
                return kmers
            else: 
                return -1
        else:
            nucleotides = {'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G'}
            kmers = []
            for i in range(self.kmer_len):
                replaced_kmer = None
                if i == self.kmer_len - 1:
                    replaced_kmer = kmer[:i] + nucleotides[kmer[i]] + kmer[i+1:]
                else:
                    replaced_kmer = kmer[:i] + nucleotides[kmer[i]] + kmer[i+1:]
                found_kmers = self.find_kmer_in_graph(replaced_kmer)
                kmers.extend(found_kmers)
            if kmers:
                return kmers
            else: 
                return -1
    
    def find_kmer_in_graph(self, kmer):
        positions = []
        if (kmer in self.debrujin['forward']):
            positions.append(self.debrujin['forward'][kmer])
        elif (kmer in self.debrujin['reverse']):
            positions.append(self.debrujin['reverse'][kmer])
        return positions