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
