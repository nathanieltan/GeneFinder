# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: YOUR NAME HERE

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    # TODO: implement this
    if nucleotide == 'A':
        return 'T'
    if nucleotide == 'T':
        return 'A'
    if nucleotide == 'G':
        return 'C'
    if nucleotide == 'C':
        return 'G'
    pass


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    # TODO: implement this
    dnaLength = len(dna)
    reverseComplement = ""
    for x in range(0,dnaLength):
        reverseComplement = reverseComplement + get_complement(dna[dnaLength - x-1])
    return reverseComplement


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    # TODO: implement this
    dnaLength = len(dna)
    result=""
    stopLocation = -1 #Location of stop codon, -1 means no stop codon
    for a in range(0,int(dnaLength/3)):
        codon = ""     
        for x in range(a*3,a*3+3):
            codon = codon+dna[x]
        if codon == "TAG" or codon =="TAA" or codon == "TGA":
            stopLocation = a  
            ORF = ""
            for b in range(0,a*3):
                ORF = ORF + dna[b]
            return ORF
         
    if stopLocation==-1:
        return dna  
    pass


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    # TODO: implement this
    myList = []
    while rest_of_ORF(dna) != dna:
        dnaLength = len(dna)
        myList.append(rest_of_ORF(dna))
        ORF = rest_of_ORF(dna)
        dna = dna[len(ORF)+3:dnaLength]
        if rest_of_ORF(dna) == dna:
            myList.append(dna)
    return myList
   


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    # TODO: implement this
    myList = [] 
    for x in range(0,3):
        dnaCopy = dna
        dnaCopy = dnaCopy[x:len(dnaCopy)]
        a = 0
        startIndex = 0
        found = False
        while a < int(len(dnaCopy)/3) and found == False:
            codon = dnaCopy[a*3:a*3+3]
            if codon == "ATG":
                startIndex = a*3
                found=True
            else:
                a = a+1
        if found == True:
            dnaCopy = dnaCopy[startIndex:len(dnaCopy)]
            if rest_of_ORF(dnaCopy) == dnaCopy:
                myList.append(dnaCopy)
            else: 
                while rest_of_ORF(dnaCopy) != dnaCopy and rest_of_ORF(dnaCopy) != "" and found == True:
                    dnaLength = len(dnaCopy)
                    myList.append(rest_of_ORF(dnaCopy))
                    ORF = rest_of_ORF(dnaCopy)
                    dnaCopy = dnaCopy[len(ORF)+3:dnaLength]
                    found = False
                    while a < int(len(dnaCopy)/3) and found == False:
                         codon = dnaCopy[a*3:a*3+3]
                         if codon == "ATG":
                             startIndex = a*3
                             found=True
                         else:
                             a = a+1
                    
                    if rest_of_ORF(dnaCopy) == dnaCopy and rest_of_ORF(dnaCopy) !="" and found == True:
                        myList.append(dnaCopy)
    return myList   



def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # TODO: implement this
    myList = []
    myList = myList+find_all_ORFs(dna)+find_all_ORFs(get_reverse_complement(dna))
    return myList
    pass


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    myList = find_all_ORFs_both_strands(dna)
    return max(myList) 


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    myList = []
    for x in range(0, num_trials):
        newDNA = shuffle_string(dna)
        myList.append(longest_ORF(newDNA))
    return len(max(myList)) 


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    acidString = ""
    length = len(dna)
    for x in range(0,int(length/3)):
        acidString = acidString + aa_table[dna[x*3:x*3+3]]
    return acidString

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    threshold = longest_ORF_noncoding(dna, 1500)
    print(threshold)
    myList = []
    ORFList = find_all_ORFs_both_strands(dna)  
    for x in range(0, len(ORFList)):
        if len(ORFList[x])>threshold:
            myList.append(coding_strand_to_AA(ORFList[x]))
    return myList
if __name__ == "__main__":
    import doctest
    doctest.run_docstring_examples(coding_strand_to_AA, globals())
   # doctest.testmod()
