from sequenceAnalysis import FastAreader, NucParams 

def main (fileName=None):
    """
    If fileName is None, FastAreader automatically reads from sys.stdin;
    otherwise, from the specified file.
    """
    myReader = FastAreader(fileName) 
    myNuc = NucParams()

    for head, seq in myReader.readFasta() :
        myNuc.addSequence(seq)
    
    #get and print length of sequence
    nucCount = 0
    for nucI in myNuc.nucleotide_comp:
        nucCount += myNuc.nucleotide_comp[nucI]
    sequenceLength = nucCount/1_000_000
    print(f"sequence length = {sequenceLength:.2f} Mb\n", end="\n")

    #get and print the percentage of G and C are in the sequence
    gcCount = myNuc.nucleotide_comp.get('G', 0) + myNuc.nucleotide_comp.get('C', 0)
    gcContent = (gcCount / nucCount) * 100 if nucCount else 0
    print(f"GC content = {gcContent:.1f}%\n")

    # sort codons in alpha order, by Amino Acid and then uses codon counts
    # Create a dictionary sorted by amino acid, then by codon
    sorted_codons = dict(sorted(myNuc.codon_comp.items(), key=lambda item: (myNuc.rnaCodonTable.get(item[0], ''), item[0])))

    
    # calculate relative codon usage for each codon and print
    for codon, codon_counts in sorted_codons.items():
        aa = myNuc.rnaCodonTable[codon]
        total_for_aa = myNuc.aaComposition().get(aa, 0)
        val = (codon_counts / total_for_aa) if total_for_aa else 0
        print("{:s} : {:s} {:5.1f} ({:6d})".format(codon, aa, val * 100, codon_counts))

if __name__ == "__main__":
    main() # make sure to change this in order to use stdin


'''
Inspection: 
I did all of my code basically in my init of sequence analysis and simply returned the dictionaries in the functions. 
I talked to my group and decided on this because it made more sense to me. 
'''
    