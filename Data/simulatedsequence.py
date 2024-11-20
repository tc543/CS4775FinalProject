import random

currstate = None
u = 0.01

#make hidden state sequence

#initialization
firststate = random.choices(["h","l"], [0.5,0.5], k=1)
currstate = firststate[0]

hiddensequence = ""
#loop
for i in range(1,1001):
    hiddensequence += currstate
    if currstate == "h":
        newstate = random.choices(["h","l"], [1-u, u], k=1)
        currstate = newstate[0]
    
    if currstate == "l":
        newstate = random.choices(["l","h"], [1-u, u], k=1)
        currstate = newstate[0]

#now make DNA sequence based on hidden sequence
dnasequence = ""
for i in range (0, len(hiddensequence)):
    if hiddensequence[i] == "h":
        nucleotide = random.choices(["A","C","G","T"], [0.13,0.37,0.37,0.13], k=1)
        dnasequence+= nucleotide[0]

    if hiddensequence[i] == "l":
        nucleotide = random.choices(["A","C","G","T"], [0.32,0.18,0.18,0.32], k=1)
        dnasequence+= nucleotide[0]


#put in file.

header = ">simulatedsequence5"

with open('simulatedsequence5.fasta', 'w') as file:
    # Write the header and sequence to the file
    file.write(header + '\n')  # Write the header with a newline
    file.write(dnasequence + '\n')  # Write the sequence with a newline
        


