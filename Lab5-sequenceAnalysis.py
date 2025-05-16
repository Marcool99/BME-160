class NucParams:
    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'   # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    def __init__ (self, inString=''): #takes an optional sequence of nucleotides from inString {ACGTUN}
        self.addSequence(inString)
        self.codon_comp = {
            codon: 0 for codon in [
                'UUU', 'UCU', 'UAU', 'UGU', 'UUC', 'UCC', 'UAC', 'UGC',
                'UUA', 'UCA', 'UAA', 'UGA', 'UUG', 'UCG', 'UAG', 'UGG',
                'CUU', 'CCU', 'CAU', 'CGU', 'CUC', 'CCC', 'CAC', 'CGC',
                'CUA', 'CCA', 'CAA', 'CGA', 'CUG', 'CCG', 'CAG', 'CGG',
                'AUU', 'ACU', 'AAU', 'AGU', 'AUC', 'ACC', 'AAC', 'AGC',
                'AUA', 'ACA', 'AAA', 'AGA', 'AUG', 'ACG', 'AAG', 'AGG',
                'GUU', 'GCU', 'GAU', 'GGU', 'GUC', 'GCC', 'GAC', 'GGC',
                'GUA', 'GCA', 'GAA', 'GGA', 'GUG', 'GCG', 'GAG', 'GGG'
            ]
        }
        self.nucleotide_comp = {'A':0, 'C':0, 'G':0, 'T':0, 'U':0, 'N':0}
        self.amino_comp = {aa: 0 for aa in "ACDEFGHIKLMNPQRSTVWY-"}
        
            
        
    def addSequence (self, inSeq):
        self.sequence = inSeq.strip().upper()

        #putting nucleotides into nucleotide_comp
        for nucleotide in self.sequence:
            if nucleotide in self.nucleotide_comp: #for each nucleutide in sequence, if it is also in self.nucleotide_comp it is real and 
                self.nucleotide_comp[nucleotide] += 1 #increment the value of the nucleotide
        
        #putting codons in codon codon_comp
        self.rna_seq = inSeq.replace('T', 'U') #replace T with U temporarily so that we can use RNA codon table exclusively for codons and amino acids
        for n in range(0, len(self.rna_seq)-2, 3): #starting at the beginning of rna_seq until the end by 3s to get the codons
            if self.rna_seq[n: n+3] in self.codon_comp: #if the codon is in the codon_comp, increment the count
                self.codon_comp[self.rna_seq[n: n+3]] +=1

        #putting amino acids in amino_comp
        for n in range(0, len(self.rna_seq)-2, 3): #starting at the beginning of rna_seq until the end by 3s to get the amino acids corresponding with the codons
            if self.rna_seq[n: n+3] in self.rnaCodonTable.keys(): #if the codon is in the codon_comp, get the corresponding amino acid by using the codon and 
                codon = self.rna_seq[n: n+3]
                self.amino_comp[self.rnaCodonTable[codon]] +=1 #then increment the amino acid 
        return self.sequence
    
    def aaComposition(self):
        return self.amino_comp
    def nucComposition(self):
        return self.nucleotide_comp
    def codonComposition(self):
        return self.codon_comp
    def nucCount(self):
        counter = 0
        for nucleotide in self.nucleotide_comp.keys():
            counter += self.nucleotide_comp[nucleotide] #adds the count of each nucleotide
        return counter
    
import sys
class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=None):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is None:
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence
#!/usr/bin/env python3
# Name: Marco Lanza (mnlanza)
# Group Members: Adrian Flores (aflor115), Noah Kirsch (nkirsch), Chris


class ProteinParam :
    ''' proteinParam class has different functions, allowing you to get the count of 
    amino acids, the theoretical isoelectric point, the composition of amino acids in 
    the protein, the charge of the protein, the molar extinction coefficient,
     the mass extinction coefficient, and the molecular weight of the protein.'''
    # These tables are for calculating:
    #     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
    #     absorbance at 280 nm (aa2abs280)
    #     pKa of positively charged Amino Acids (aa2chargePos)
    #     pKa of negatively charged Amino acids (aa2chargeNeg)
    #     and the constants aaNterm and aaCterm for pKa of the respective termini
    #  Feel free to move these to appropriate methods as you like

    # As written, these are accessed as class attributes, for example:
    # ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__ (self, protein): #done?
        '''initializes the protein object and makes a dictionary amino_composition
        with how much of each amoino acid the protein has'''
        self.protein = protein.upper()
        self.amino_composition = {
            'A':0, 'C':0, 'D':0, 'E':0, 'F':0, 'G':0, 'H':0, 'I':0, 
            'L':0, 'K':0, 'M':0, 'N':0, 'P':0, 'Q':0, 'R':0, 'S':0, 
            'T':0, 'V':0, 'Y':0, 'W':0
        }
        for char in protein:
            if char in self.amino_composition.keys():
                self.amino_composition[char] += 1
            

    def aaCount (self): #done?
        '''returns a single integer count of valid amino acid characters found'''
        count = 0
        for char in self.protein:
            if char in self.amino_composition:
                count += 1
        return count

    def pI (self, precision=2):
        '''takes the protein as an input and outputs the theoretical isoelectric point of the protein'''
        low = 0.0
        high = 14.0
        
        # Define the threshold for zero charge, based on the precision
        threshold = 10 ** (-precision)
        
        # Perform the binary search
        while (high - low) > threshold:
            mid = (low + high) / 2.0
            charge = self._charge_(mid)
        
            if charge > 0:
                # If charge is positive, look in the higher pH range
                low = mid
            elif charge < 0:
                # If charge is negative, look in the lower pH range
                high = mid
            else:
                # If charge is exactly 0, we've found the pI
                return round(mid, precision)

        # Return the midpoint rounded to the specified precision
        return round((low + high) / 2.0, precision)

    def aaComposition (self) : #done?
        ''' this function returns the counts of each amino acid in the protein''' 
        return self.amino_composition

    def _charge_ (self, pH): #done?
        '''takes charge and protein as inputs and outputs the net charge of the protein'''
        net_charge = 0.0
        net_charge += 10**(self.aaNterm - pH) / (10**(self.aaNterm - pH) + 1)
        
        # C-terminus: -1 group
        net_charge -= 1.0 / (10**(self.aaCterm - pH) + 1)
    
        # Now handle side chains
        for aa, count in self.amino_composition.items():
            if count == 0:
                continue
            # Basic side chains (K, R, H) => add +1 * fraction
            if aa in self.aa2chargePos:
                pKa = self.aa2chargePos[aa]
                fraction = 10**(pKa - pH) / (10**(pKa - pH) + 1)
                net_charge += count * fraction
            
            # Acidic side chains (D, E, C, Y) => add -1 * fraction
            if aa in self.aa2chargeNeg:
                pKa = self.aa2chargeNeg[aa]
                fraction = 1.0 / (10**(pKa - pH) + 1)
                net_charge -= count * fraction
        
        return net_charge

    
    def molarExtinction (self, cystine=True):
        '''outputs the molar extinction coefficient'''
        coefficient = 0
        nY = self.amino_composition['Y']  # count of Tyr
        nW = self.amino_composition['W']  # count of Trp
        nC = self.amino_composition['C']  # count of Cys
        if cystine:
            coefficient = nY*self.aa2abs280['Y'] + nW*self.aa2abs280['W'] + nC*self.aa2abs280['C']
        else:
            coefficient = nY*self.aa2abs280['Y'] + nW*self.aa2abs280['W']
        return coefficient


    def massExtinction (self, cystine=True):
        '''Outputs the mass extinction coefficient given the protein using the molarExtinction function'''
        myMW = self.molecularWeight()
        if myMW == 0:
            return 0.0
        return self.molarExtinction(cystine=cystine) / myMW

    def molecularWeight (self): #done?
        '''outputs the mollecular weight of the protein'''
        sum = 0
        for aa, aa_count in self.amino_composition.items():
            sum += (self.aa2mw[aa]*aa_count - self.mwH2O*aa_count)
        sum += self.mwH2O
        return sum
    

    """
    Finds Open Reading Frames (ORFs) in DNA sequences.
    
    Parameters:
    - lG (bool): If True, only keeps the longest gene for each stop codon
    - mG (int): Minimum gene length to report
    - s (list): List of start codons
    - t (list): List of stop codons
    """
    def __init__(self, lG=False, mG=100, s=None, t=None):
        self.longest_gene = lG
        self.minimum_gene = mG
        self.start_codons = s if s is not None else ['ATG']
        self.stop_codons = t if t is not None else ['TAG', 'TGA', 'TAA']
        self.orfs = []  # List containing all found ORFs
    
    def reverse_complement(self, sequence):
        """Returns the reverse complement of a DNA sequence"""
        comp_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(comp_map[base] for base in reversed(sequence))
    
    def format_orf_data(self, strand, frame, start_pos, stop_pos, seq_length):
        """
        Format ORF data with standardized coordinates
        
        Returns:
        [frame, start, stop, length] where:
        - frame is 1,2,3 for forward strand, -1,-2,-3 for reverse strand
        - start/stop are 1-based coordinates
        - length is the ORF length in nucleotides
        """
        orf_length = stop_pos - start_pos
        
        if strand == '+':
            out_frame = frame + 1
            out_start = start_pos + 1
            out_stop = stop_pos
        else:  # strand == '-'
            out_frame = -1 * (frame + 1)
            out_start = (seq_length - stop_pos) + 1
            out_stop = (seq_length - start_pos)
            
        return [out_frame, out_start, out_stop, orf_length]
    
    def find_orfs(self, sequence):
        """
        Find all possible ORFs in a DNA sequence in both strands and all reading frames.
        Results are stored in self.orfs
        """
        # Clean and prepare the sequence
        sequence = sequence.upper().replace(' ', '')
        
        # Process both strands
        strands = {'+': sequence, '-': self.reverse_complement(sequence)}
        
        self.orfs = []  # Reset the ORFs list
        
        # Process each strand and reading frame
        for strand, seq in strands.items():
            for frame in range(3):
                self._find_orfs_in_frame(seq, strand, frame)
                
        # Filter and sort the ORFs according to settings
        self._filter_orfs()
    
    def _find_orfs_in_frame(self, sequence, strand, frame):
        """Find ORFs in a specific strand and reading frame"""
        seq_length = len(sequence)
        empty_frame = True
        stop_upstream = False
        su_hanging_start = False
        
        # Scan the sequence in the current frame
        for i in range(frame, seq_length, 3):
            current_codon = sequence[i:i+3]
            
            # Process stop codons
            if len(current_codon) == 3 and current_codon in self.stop_codons:
                stop_upstream = True
                empty_frame = False
                
                # Check for hanging stop
                hanging_stop = True
                for k in range(i, -1, -3):
                    check_codon = sequence[k:k+3]
                    # Not a hanging stop if preceded by another stop or if first codon is a start
                    if (k != i) and ((check_codon in self.stop_codons) or (sequence[0:3] in self.start_codons)):
                        hanging_stop = False
                        break
                
                if hanging_stop:
                    orf_data = self.format_orf_data(strand, frame, 0, i+3, seq_length)
                    self.orfs.append(orf_data)
            
            # Process start codons
            if len(current_codon) == 3 and current_codon in self.start_codons:
                empty_frame = False
                hanging_start = True
                
                # Look for matching stop codon
                for j in range(i, seq_length, 3):
                    check_stop = sequence[j:j+3]
                    if len(check_stop) == 3 and check_stop in self.stop_codons:
                        hanging_start = False
                        # Regular ORF
                        orf_data = self.format_orf_data(strand, frame, i, j+3, seq_length)
                        self.orfs.append(orf_data)
                        
                        # Check for possible upstream hanging start
                        if not stop_upstream and sequence[0:3] not in self.start_codons:
                            orf_data = self.format_orf_data(strand, frame, 0, j+3, seq_length)
                            self.orfs.append(orf_data)
                        break
                
                # Handle hanging start (no stop found)
                if hanging_start:
                    # Regular hanging start
                    orf_data = self.format_orf_data(strand, frame, i, seq_length, seq_length)
                    self.orfs.append(orf_data)
                    
                    # Check for upstream hanging start
                    if (not stop_upstream and 
                        sequence[0:3] not in self.start_codons and 
                        not su_hanging_start):
                        su_hanging_start = True
                        orf_data = self.format_orf_data(strand, frame, 0, seq_length, seq_length)
                        self.orfs.append(orf_data)
        
        # Handle completely empty frames (no start/stop codons)
        if empty_frame:
            orf_data = self.format_orf_data(strand, frame, 0, seq_length, seq_length)
            self.orfs.append(orf_data)
    
    def _filter_orfs(self):
        """Filter ORFs based on minimum length and longest gene setting"""
        # Sort by length (descending)
        sorted_orfs = sorted(self.orfs, key=lambda x: x[3], reverse=True)
        
        filtered_orfs = []
        reference_positions = set()
        
        for orf in sorted_orfs:
            frame, start, stop, length = orf
            
            # Skip if shorter than minimum length
            if length < self.minimum_gene:
                continue
                
            if self.longest_gene:
                # For positive frames, track by stop position
                if frame > 0:
                    if stop not in reference_positions:
                        filtered_orfs.append(orf)
                        reference_positions.add(stop)
                # For negative frames, track by start position
                else:
                    if start not in reference_positions:
                        filtered_orfs.append(orf)
                        reference_positions.add(start)
            else:
                # If not longest_gene, keep all that meet minimum length
                filtered_orfs.append(orf)
                
        self.orfs = filtered_orfs
    
    def clear_data(self, filename=''):
        """Clear the contents of an output file"""
        with open(filename, 'w') as f:
            pass
    
    def write_data(self, header=''):
        """Print ORF data with formatting"""
        print(header)
        for orf in self.orfs:
            frame, start, stop, length = orf
            print('{:+d} {:>5d}..{:>5d} {:>5d}'.format(frame, start, stop, length))



    """
    Finds Open Reading Frames (ORFs) in DNA sequences.
    
    An ORF is a region that starts with a start codon and ends with a stop codon,
    with both codons in the same reading frame. This class handles finding ORFs
    in all six possible reading frames (3 on forward strand, 3 on reverse strand).
    
    If no start or stop codons are found in a given reading frame, the entire
    sequence in that frame is reported as an ORF.
    
    Parameters:
    - lG (bool): If True, only reports the longest gene in each ORF. Default is False.
    - mG (int): Minimum gene length to report (including start and stop codons). Default is 100.
    - s (list): List of valid start codons. Default is ['ATG'].
    - t (list): List of valid stop codons. Default is ['TAG', 'TGA', 'TAA'].
    """
    def __init__(self, lG=False, mG=100, s=['ATG'], t=['TAG', 'TGA', 'TAA']):
        self.longestGene = lG
        self.minimumGene = mG
        self.startCodons = s
        self.stopCodons = t
        self.orfs = []  # List to store found ORFs
    
    def reverse_complement(self, sequence):
        """
        Generate the reverse complement of a DNA sequence.
        
        Args:
            sequence (str): The DNA sequence to convert
            
        Returns:
            str: The reverse complement of the input sequence
        """
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(complement.get(base, base) for base in reversed(sequence))
    
    def find_orfs(self, sequence):
        """
        Find all ORFs in a DNA sequence for all six reading frames.
        
        Args:
            sequence (str): The DNA sequence to analyze
        """
        # Clear previous results
        self.orfs = []
        
        # Normalize sequence
        sequence = sequence.upper().replace(' ', '')
        
        # Process both forward and reverse strands
        for direction, seq in {'+': sequence, '-': self.reverse_complement(sequence)}.items():
            # Process each reading frame (0, 1, 2)
            for frame in range(3):
                self._find_orfs_in_frame(seq, direction, frame)
        
        # Apply filtering and sorting
        self._filter_and_sort_orfs()
    
    def _find_orfs_in_frame(self, sequence, direction, frame):
        """
        Find ORFs in a specific reading frame of a sequence.
        
        Args:
            sequence (str): The DNA sequence to analyze
            direction (str): '+' for forward strand, '-' for reverse strand
            frame (int): Reading frame (0, 1, or 2)
        """
        seq_len = len(sequence)
        
        # Check if there are any start or stop codons in this frame
        has_start = False
        has_stop = False
        
        for i in range(frame, seq_len - 2, 3):
            codon = sequence[i:i+3]
            if len(codon) == 3:  # Make sure it's a complete codon
                if codon in self.startCodons:
                    has_start = True
                if codon in self.stopCodons:
                    has_stop = True
                    
                if has_start and has_stop:
                    break
        
        # If no start or stop codons, report the entire frame as an ORF
        if not has_start and not has_stop:
            orf_data = self._format_orf_data(direction, frame, 0, seq_len, seq_len)
            self.orfs.append(orf_data)
            return
        
        # List to hold potential genes in this frame
        potential_genes = []
        
        # Check for gene fragment at 5' end (no start codon)
        leading_fragment = self._check_leading_fragment(sequence, frame)
        if leading_fragment:
            potential_genes.append(leading_fragment)
        
        # Find all internal genes (with both start and stop)
        internal_genes = self._find_internal_genes(sequence, frame)
        potential_genes.extend(internal_genes)
        
        # Check for gene fragment at 3' end (no stop codon)
        trailing_fragment = self._check_trailing_fragment(sequence, frame)
        if trailing_fragment:
            potential_genes.append(trailing_fragment)
        
        # If no genes found (could happen if all codons are too short), report the entire frame
        if not potential_genes:
            orf_data = self._format_orf_data(direction, frame, 0, seq_len, seq_len)
            self.orfs.append(orf_data)
            return
        
        # Format and add all found genes to the results
        for start_pos, end_pos in potential_genes:
            orf_data = self._format_orf_data(direction, frame, start_pos, end_pos, seq_len)
            self.orfs.append(orf_data)
    
    def _check_leading_fragment(self, sequence, frame):
        """
        Check for gene fragment at the 5' end of the sequence (fragment with stop but no start).
        
        Args:
            sequence (str): The DNA sequence to analyze
            frame (int): Reading frame (0, 1, or 2)
            
        Returns:
            tuple or None: (0, stop_pos) if fragment found, None otherwise
        """
        # Start from the first codon in this frame
        for i in range(frame, len(sequence), 3):
            codon = sequence[i:i+3]
            if len(codon) < 3:  # Incomplete codon at end
                break
                
            # If we hit a start codon first, no leading fragment
            if codon in self.startCodons:
                return None
                
            # If we hit a stop codon, we have a leading fragment
            if codon in self.stopCodons:
                return (0, i+3)
        
        return None
    
    def _find_internal_genes(self, sequence, frame):
        """
        Find all internal genes (with both start and stop codons) in a reading frame.
        
        Args:
            sequence (str): The DNA sequence to analyze
            frame (int): Reading frame (0, 1, or 2)
            
        Returns:
            list: List of (start_pos, end_pos) tuples for all genes found
        """
        genes = []
        
        # Track possible start positions
        start_positions = []
        
        # Scan the sequence in this frame
        for i in range(frame, len(sequence), 3):
            if i + 3 > len(sequence):  # Incomplete codon at end
                break
                
            codon = sequence[i:i+3]
            
            # Record position if start codon
            if codon in self.startCodons:
                start_positions.append(i)
            
            # If stop codon found, create genes from all recorded start positions
            elif codon in self.stopCodons and start_positions:
                end_pos = i + 3
                
                # Add a gene for each start position
                for start_pos in start_positions:
                    genes.append((start_pos, end_pos))
                
                # Clear start positions if longestGene is True
                # Otherwise keep them for potential overlapping genes
                if self.longestGene:
                    start_positions = []
        
        return genes
    
    def _check_trailing_fragment(self, sequence, frame):
        """
        Check for gene fragment at the 3' end of the sequence (fragment with start but no stop).
        
        Args:
            sequence (str): The DNA sequence to analyze
            frame (int): Reading frame (0, 1, or 2)
            
        Returns:
            tuple or None: (start_pos, len(sequence)) if fragment found, None otherwise
        """
        # Scan from end toward beginning to find the last start codon with no stop after it
        last_start = None
        
        for i in range(frame, len(sequence), 3):
            if i + 3 > len(sequence):  # Incomplete codon at end
                break
                
            codon = sequence[i:i+3]
            
            if codon in self.startCodons:
                last_start = i
            elif codon in self.stopCodons:
                last_start = None
        
        if last_start is not None:
            return (last_start, len(sequence))
        
        return None
    
    def _format_orf_data(self, direction, frame, start_pos, end_pos, seq_len):
        """
        Format ORF data with correct coordinates.
        
        Args:
            direction (str): '+' for forward strand, '-' for reverse strand
            frame (int): Reading frame (0, 1, or 2)
            start_pos (int): Start position in the sequence (0-based)
            end_pos (int): End position in the sequence (0-based)
            seq_len (int): Length of the sequence
            
        Returns:
            list: [frame, start_pos, end_pos, length] with 1-based coordinates
        """
        length = end_pos - start_pos
        
        if direction == '+':
            out_frame = frame + 1
            out_start = start_pos + 1
            out_end = end_pos
        else:  # direction == '-'
            out_frame = -(frame + 1)
            out_start = seq_len - end_pos + 1
            out_end = seq_len - start_pos
        
        return [out_frame, out_start, out_end, length]
    
    def _filter_and_sort_orfs(self):
        """
        Filter ORFs based on criteria and sort by length (decreasing) and position.
        """
        
        if self.longestGene:
            # Track ORFs by their end coordinates to ensure we only keep the longest
            unique_orfs = {}
            
            for orf in self.orfs:
                frame, start, end, length = orf
                
                # Key is the end position for positive frames, start position for negative frames
                key = end if frame > 0 else start
                
                # Keep this ORF if it's longer than any previously seen ORF with the same key
                if key not in unique_orfs or length > unique_orfs[key][3]:
                    unique_orfs[key] = orf
            
            # Convert back to list
            self.orfs = list(unique_orfs.values())
        
        # Filter by minimum length
        self.orfs = [orf for orf in self.orfs if orf[3] >= self.minimumGene]
        
        # Sort by length (descending), then by start position (ascending)
        self.orfs.sort(key=lambda x: (-x[3], x[1]))
        
    
    def write_data(self, header='', orf_data=None):
        """
        Write ORF data to standard output.
        
        Args:
            header (str): Header for the sequence
            orf_data (list): List of ORF data to write. If None, uses self.orfs
        """
        if orf_data is None:
            orf_data = self.orfs
            
        print(header)
        for orf in orf_data:
            frame, start, stop, length = orf
            print('{:+d} {:>5d}..{:>5d} {:>5d}'.format(frame, start, stop, length))
    
    def clear_data(self, fileName=''):
        """
        Clear the contents of a file.
        
        Args:
            fileName (str): Name of the file to clear
        """
        with open(fileName, 'w') as f:
            pass

class OrfFinder:
    """
    Finds Open Reading Frames (ORFs) in DNA sequences.
    
    Parameters:
    - lG (boolean): If True, only keeps the longest gene for each stop codon 
    - mG (integer): Minimum gene length to report
    - s (list): List of start codons
    - t (list): List of stop codons
    """
    def __init__(self, lG=False, mG=100, s=None, t=None): # longest gene default false and mG defalut false
        self.longest_gene = lG
        self.minimum_gene = mG
        self.start_codons = s if s is not None else ['ATG'] # default start codons ATG
        self.stop_codons = t if t is not None else ['TAG', 'TGA', 'TAA'] #defalut stop codons ['TAG', 'TGA', 'TAA'] 
        self.orfs = []  # List containing all found ORFs
    
    def reverse_complement(self, sequence):
        """Returns the reverse complement of a DNA sequence (bottom strand of DNA)"""
        comp_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(comp_map[base] for base in reversed(sequence))
    
    def format_orf_data(self, strand, frame, start_pos, stop_pos, seq_length):
        """
        Format ORF data with standardized coordinates
        
        Returns:
        [frame, start, stop, length] where:
        - frame is 1,2,3 for forward strand, -1,-2,-3 for reverse strand
        - length is the ORF length in nucleotides (8->13 is 6 nucleotides long)
        """
        orf_length = stop_pos - start_pos
        
        if strand == '+':
            # Frame is 1-based (1, 2, 3)
            out_frame = frame + 1
            # Convert to 1-based coordinates
            out_start = start_pos + 1
            out_stop = stop_pos
        else:  # strand == '-'
            # Frame is negative 1-based (-1, -2, -3)
            out_frame = -1 * (frame + 1)
            # Convert to 1-based coordinates in original sequence
            out_start = (seq_length - stop_pos) + 1
            out_stop = (seq_length - start_pos)
            
        return [out_frame, out_start, out_stop, orf_length]
    
    def find_orfs(self, sequence):
        """
        Find all possible ORFs in a DNA sequence in both strands and all reading frames.
        Results are stored in self.orfs
        """
        # Clean and prepare the sequence
        sequence = sequence.upper().replace(' ', '')
        
        # Process both strands
        strands = {'+': sequence, '-': self.reverse_complement(sequence)}
        
        self.orfs = []  # Reset the ORFs list
        
        # Process each strand and reading frame
        for strand, seq in strands.items():
            for frame in range(3):
                self._find_orfs_in_frame(seq, strand, frame)
                
        # Filter and sort the ORFs according to settings
        self._filter_orfs()
    
    def _find_orfs_in_frame(self, sequence, strand, frame):
        """Find ORFs in a specific strand and reading frame"""
        seq_length = len(sequence)
        empty_frame = True
        stop_upstream = False
        su_hanging_start = False
        
        # Scan the sequence in the current frame
        for i in range(frame, seq_length, 3):
            current_codon = sequence[i:i+3]
            
            # Skip incomplete codons at the end
            if len(current_codon) < 3:
                continue
            
            # Process stop codons
            if current_codon in self.stop_codons:
                stop_upstream = True
                empty_frame = False
                
                # Check for hanging stop
                hanging_stop = True
                for k in range(i, -1, -3):
                    if k < 0 or k + 3 > len(sequence):
                        continue
                    check_codon = sequence[k:k+3]
                    # Not a hanging stop if preceded by another stop or if first codon is a start
                    if (k != i) and ((check_codon in self.stop_codons) or (k == frame and check_codon in self.start_codons)):
                        hanging_stop = False
                        break
                
                if hanging_stop:
                    orf_data = self.format_orf_data(strand, frame, 0, i+3, seq_length)
                    self.orfs.append(orf_data)
            
            # Process start codons
            if current_codon in self.start_codons:
                empty_frame = False
                hanging_start = True
                
                # Look for matching stop codon
                for j in range(i, seq_length, 3):
                    if j + 3 > seq_length:
                        break
                    check_stop = sequence[j:j+3]
                    if check_stop in self.stop_codons:
                        hanging_start = False
                        # Regular ORF
                        orf_data = self.format_orf_data(strand, frame, i, j+3, seq_length)
                        self.orfs.append(orf_data)
                        
                        # Check for possible upstream hanging start
                        if not stop_upstream and (frame > 0 or sequence[frame:frame+3] not in self.start_codons):
                            orf_data = self.format_orf_data(strand, frame, 0, j+3, seq_length)
                            self.orfs.append(orf_data)
                        break
                
                # Handle hanging start (no stop found)
                if hanging_start:
                    # Regular hanging start
                    orf_data = self.format_orf_data(strand, frame, i, seq_length, seq_length)
                    self.orfs.append(orf_data)
                    
                    # Check for upstream hanging start
                    if (not stop_upstream and 
                        (frame > 0 or sequence[frame:frame+3] not in self.start_codons) and 
                        not su_hanging_start):
                        su_hanging_start = True
                        orf_data = self.format_orf_data(strand, frame, 0, seq_length, seq_length)
                        self.orfs.append(orf_data)
        
        # Handle completely empty frames (no start/stop codons)
        if empty_frame:
            orf_data = self.format_orf_data(strand, frame, 0, seq_length, seq_length)
            self.orfs.append(orf_data)
    
    def _filter_orfs(self):
        """Filter ORFs based on minimum length and longest gene setting"""
        # Sort by length (descending)
        sorted_orfs = sorted(self.orfs, key=lambda x: x[3], reverse=True)
        
        filtered_orfs = []
        reference_positions = set()
        
        for orf in sorted_orfs:
            frame, start, stop, length = orf
            
            # Skip if shorter than minimum length
            if length < self.minimum_gene:
                continue
                
            if self.longest_gene:
                # For positive frames, track by stop position
                if frame > 0:
                    ref_key = (frame, stop)
                    if ref_key not in reference_positions:
                        filtered_orfs.append(orf)
                        reference_positions.add(ref_key)
                # For negative frames, track by start position
                else:
                    ref_key = (frame, start)
                    if ref_key not in reference_positions:
                        filtered_orfs.append(orf)
                        reference_positions.add(ref_key)
            else:
                # If not longest_gene, keep all that meet minimum length
                filtered_orfs.append(orf)
                
        self.orfs = filtered_orfs
    
    def clear_data(self, filename=''):
        """Clear the contents of an output file using write (in file)"""
        with open(filename, 'w') as f:
            pass
    
    def write_data(self, header='', orf_data=None):
        """Print ORF data with formatting"""
        if orf_data is None:
            orf_data = self.orfs
            
        print(header)
        for orf in orf_data:
            frame, start, stop, length = orf
            print('{:+d} {:>5d}..{:>5d} {:>5d}'.format(frame, start, stop, length))
