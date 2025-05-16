from sequenceAnalysis import FastAreader#, OrfFinder 

#!/usr/bin/env python3
# Name: Marco Lanza (mnlanza)
# Group Members: Noah Kirsch, Adrian Flores, Chris Quach


########################################################################
# CommandLine
########################################################################
class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
 
    '''
    
    def __init__(self, inOpts=None) : #don't touch command line
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does', 
                                             epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1[default] <input >output'
                                             )
        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=False, help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices= (0,100,200,300,500,1000), default=100, action = 'store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action = 'append', default = ['ATG'],nargs='?', 
                                 help='start Codon') #allows multiple list options
        self.parser.add_argument('-t', '--stop', action = 'append', default = ['TAG','TGA','TAA'],nargs='?', help='stop Codon') #allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

########################################################################
# Main
# Here is the main program
# 
#
########################################################################
   
from sequenceAnalysis import OrfFinder as OrfFinder
from sequenceAnalysis import FastAreader as FastAreader

def main(inFile=None, options=None):
    thisCommandLine = CommandLine(options)
    
    # Get arguments from command line
    lG = thisCommandLine.args.longestGene
    mG = thisCommandLine.args.minGene
    s = thisCommandLine.args.start
    t = thisCommandLine.args.stop
    
    # Create OrfFinder instance object
    myOrf = OrfFinder(lG, mG, s, t)
    
    # Read FASTA file and process each sequence
    reader = FastAreader(inFile)
    for head, seq in reader.readFasta():
        # Reset orfs for each sequence
        myOrf.orfs = []
        
        # Find ORFs
        myOrf.find_orfs(seq)
        
        # Write output to STDOUT
        myOrf.write_data(head, myOrf.orfs)

if __name__ == "__main__":
    main()
    #main(inFile = 'tass2.fa', options = ['-mG=300', '-lG']) # delete this stuff if running from commandline



'''
Inspection:
After talking with my groupmates, I ended up using a lot of helper functions to make my code more readable and understandable for me. 
These functions allowed me to deal with edge cases one by one and get through the project.
'''




