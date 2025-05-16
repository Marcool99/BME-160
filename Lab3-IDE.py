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

# Please do not modify any of the following.  This will produce a standard output that can be parsed
    
import sys
def main():
    '''prints all of the data from the functions above in the correct format'''
    inString = input('protein sequence?')
    while inString :
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        print ("Amino acid composition:")
        
        if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present 
        
        for aa,n in sorted(myParamMaker.aaComposition().items(), 
                           key= lambda item:item[0]):
            print ("\t{} = {:.2%}".format(aa, n/myAAnumber))
    
        inString = input('protein sequence?')
#input = 'VLSPADKTNVKAAW'
if __name__ == "__main__":
    main()