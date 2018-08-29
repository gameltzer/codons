
# coding: utf-8

# In[ ]:




# In[3]:

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Alphabet import IUPAC




class NotFastaOrFastqFileError(Exception): 
    """This class is an exception will be generated if the file type is incorrect.
       The instantiation is overriden and a value (the file extension) is passed to an instance variable.
       This is instance variable 'value' is returned when a string representation is requested. 
    """
    def __init__(self, value):
        self.value = value
    
    def __repr__(self):
        return repr(self.value)

"""This function gets the file extension from the filePath. If no file 
    extension is found, returns the message "Not a file". This is what will
    be passed to the exception handling if it is not a FASTA or FASTQ file.
    
"""
def getFileExtension(filePath):
    if '.' in filePath:
        """Reverses the filePath and then finds the first (i.e., the last)'.' in the reversed
            filePath. """
        extensionPeriodLocation = filePath[::-1].find('.')

        """The variable fileExtension is only true for the reverse string. This code makes 
            sure it applies to the standard orientation. Because in a slice the execution
            never reaches the item in the middle position, we must add 1 to it so we get 
            the period. Since we are moving backwards, we are converting it to a negative
            number. We then use a second slice operation to get the result facing foward.
        """
        fileExtension = filePath[-1:-(extensionPeriodLocation +1):-1][::-1]
        """In case there is a period in one of the sections between the 
            slashes, but the portion after the last slash contains no periods.
            In this cae, it is probably a directory if it is a valid path at
            all. """
        if fileExtension != '':
            """The lowercase function is used in case the Biopython parse 
            function requires lowercase arguments for the file type. Not sure
            how much this matters but good just to be safe. 
            """
            return fileExtension.lower()
        
        else: return "Not a file"
    # Has no file extension, so cannot  be a file
    else:
        return "Not a file"

def formatAminoAcids(count_dictionary):
    ### Amino acid names  gottens from http://www.chem.qmul.ac.uk/iupac/AminoAcid/AA1n2.html 
    # key of new dictionary is amino acid name, value is amino acid count
    formatted_dictionary = {'Alanine': count_dictionary['A']}
    formatted_dictionary['Arginine'] = count_dictionary['R']
    formatted_dictionary['Asparagnine'] = count_dictionary['N']
    formatted_dictionary['Aspartic acid'] = count_dictionary['D']
    formatted_dictionary['Cysteine'] = count_dictionary['C']
    ### 5 done
    formatted_dictionary['Glutamine'] = count_dictionary['Q']
    formatted_dictionary['Glutamic acid'] = count_dictionary['E']
    formatted_dictionary['Glycine'] = count_dictionary['G']
    formatted_dictionary['Histidine'] = count_dictionary['H']
    formatted_dictionary['Isoleucine'] = count_dictionary['I']
    # 10 done
    formatted_dictionary['Leucine'] = count_dictionary['L']
    formatted_dictionary['Lysine'] = count_dictionary['K']
    formatted_dictionary['Methionine'] = count_dictionary['M']
    formatted_dictionary['Phenylalanine'] = count_dictionary['F']
    formatted_dictionary['Proline'] = count_dictionary['P']
    ### 15 done'
    formatted_dictionary['Serine'] = count_dictionary['S']
    formatted_dictionary['Threonine'] = count_dictionary['T']
    formatted_dictionary['Tryptophan'] = count_dictionary['W']
    formatted_dictionary['Tyrosine'] = count_dictionary['Y']
    formatted_dictionary['Valine'] = count_dictionary['V']
    return formatted_dictionary

""" This function actually processes data from the FASTA or FASTQ file. This is pretty much just a loop
    that iterates through all the records of the SeqIO object. It converts the input into the necessary
    forms required to count the amino acids, and then writes the relevant information to the file specified
    by "output." "iterator" should be the iterator returned by parsing a SEQIO object.
    Using the IUPACAmbiguousDNA Alphabet... our example has N in it, 
    and it is good because it allows for flexibility in future inputs, as well. This should include
    all DNA files we input. It also includes GATC, which are  the only letters  in 
    IUPACUnambiguousDNA. IUPACAmbiguousDNA allows for a greater subset. """
def processFile(iterator,output):
    ### This loop prints the protein count of all the record###
    for record in iterator:
        ### gets the sequence from the record
        thisSeq = record.seq
        '''This prevents a Biopython warning from showing up if there are incomplete codons.
            Appends an N (wildcard) to the sequence until it is divisible by 3. This is essentially what is 
            recommended in the Biopython warning message.'''
        if (len(thisSeq) % 3 != 0):
            leftoverNucCount = len(thisSeq) % 3 
            for i in range(3 - leftoverNucCount):
                thisSeq = Seq((str(thisSeq)+ "N"), thisSeq.alphabet)
        ## translates the record as a new amino acid/peptide sequence
        translated_sequence = thisSeq.translate()
        ### This changes the sequence to a ProteinAnalysis object, which lets us call the needed methods on it. ###
        analyzed_sequence = ProteinAnalysis(str(translated_sequence))
        output.write("Name: {0}\nDescription: {1}\nAnnotations: {2}".format(record.name, record.description,
                                                                            record.annotations))
        ### gets the amino acid count
        aminoAcidCountDictionary = formatAminoAcids(analyzed_sequence.count_amino_acids())
        ### prints the amino acid count!
        output.write("\n\nThis is the amino acid count of record {0}:".format(record.id) +"\n\n")
        ### splits the output so that each Amino Acid gets it's own line
        for aminoAcid, count in aminoAcidCountDictionary.items():
            output.write(aminoAcid + ": " + str(count)+"\n")
        ###print(aminoAcidCountDictionary)
        ### turns the sequence into  RNA
        thisSeqRNA= thisSeq.transcribe()
        analyzed_RNAsequence = ProteinAnalysis(str(thisSeqRNA.translate()))

        """ Since the RNA is the same as the DNA with the exception of one nucleotide, getting an amino acid
        count from the RNA should be the same as the amino acid count from the DNA.
        """
        output.write("\n\nThis is the amino acid count of the protein sequence derived from the RNA resulting from "
        "the DNA sequence in the file. It should be the same as the previous amino acid count: \n\n")
        rnaAcidCountDictionary = formatAminoAcids(analyzed_RNAsequence.count_amino_acids())

        for aminoAcid, count in rnaAcidCountDictionary.items():

            output.write(aminoAcid+ ": "+ str(count)+"\n")
        output.write('\n**************************************\n\n')


while True:
    try:
    # this gets input from the Python interpreter... We can then pass this to SeqIo to read it into Biopython
        fileInputPath = raw_input("Enter FASTA or FASTQ file to read. If this request was a mistake, " 
        "press enter with no other input. Otherwise, enter a file containing DNA sequences: \n")
        if fileInputPath == '':
            """This should make it easier to run in Anaconda so that we don't have to keep 
               restarting the kernel. If enter is pressed but no input is given, we end
               the progam.
            """      
            print("No input provided. Leaving program....")
            break
    # Uses the written getFileExtension function to get the extension from the end of the file
        fileExtension = getFileExtension(fileInputPath)
    #raises NotFastaOrFastqException if file doesn't end in fasta or fastq
        if fileExtension != 'fasta' and fileExtension != 'fastq':
            raise NotFastaOrFastqFileError(fileExtension)
        try:
            """Parsing in SeqIO does not appear to return a FileNotFoundError if the file does 
                not exist. Instead the exception is thrown later. I want to avoid using a lot of try
                statements within loops; it seems difficult to keep track of, so below, I am 
                testing the file beforehand. So a file handle is created beforehand and then closed """
            fileTest = open(fileInputPath)
            fileTest.close()
            # if the above statement is successful, we have what we need. We can leave the loop!
            break
            # if the file is not found, lets us know and allows us the opportunity to try again/.
        except FileNotFoundError as e:
            print("Try again. File not found! ",e)

    ### If the file extension is not right, prints out a message, and allows us to try again. 
    except NotFastaOrFastqFileError as e:
        print("Try again. This is not a fasta or fastq file, but rather a file of type: ", e)
#If we don't specify an input path, we shouldn't bother trying to specify output path. 
if (fileInputPath != ''):
# This allows us to specify wher are outputting the file to. 
    while True:
        try:
            # This allows us to specify where we are outputting the file to. 
            fileOutputPath = raw_input('Write the path of the file you want to output the information to!'
                                   " If this was a mistake, presss enter. Otherwise, provide the file path: \n'")

            if fileOutputPath == '':
                """This should make it easier to run in Anaconda so that we don't have to keep 
                   restarting the kernel. If enter is pressed but no input is given, we end the program. 
                """      
                print("No output path provided. Leaving program....")
                break
            else:
                output_file = open(fileOutputPath,'a') 
                break
        except IOError:
            print(fileOutputPath + "is not a valid file path. Try again")
### We should stil define fileOutputPath even if we don't have an input path to prevent errors.
else:
    fileOutputPath = ''
"""we only need to check the fileOutputPath since it is set to empty if fileInputPath is empty.
If the fileInputPath is entered, we want to know if the outputPath has been entered. If it has, it is not empty and 
we can proceed. If it hasn't, it will be empty and we should reach the end of the program. If the fileInputPath is 
not entered fileOutputPath is set to empty."""
if (fileOutputPath != ''):
    print("Outputting information to {0}...".format(fileOutputPath))
            
    """ This creates an iterator that will return SeqRecord objects from the file we have selected, via the path we input from
    the interpreter..
    This will place less demands on the memory than turning the results of the parse into a list. 
    fileInputPath and fileExtension are passed from the while loop above after having been
    determined to be valid input.
    """
    record_iterator = SeqIO.parse(fileInputPath, fileExtension, IUPAC.ambiguous_dna)
    processFile(record_iterator, output_file)
    ### This closes the output file    
    output_file.close()
    print("Done!")


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:



