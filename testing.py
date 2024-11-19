import algorithms
import numpy as np

'''Reads the fasta file and outputs the sequence to analyze.
Arguments:
	filename: name of the fasta file
Returns:
	s: string with relevant sequence
'''
def read_fasta(filename):
    with open(filename, "r") as f:
        s = ""
        s += f.readlines()[0].strip()
    return s

def test1():
    mu = 0.05
    transition_probabilities = {
        'h': {'h': np.log(1 - mu), 'l': np.log(mu)},
        'l': {'h': np.log(mu), 'l': np.log(1 - mu)}
    }
    emission_probabilities = {
        'h': {'A': np.log(0.13), 'C': np.log(0.37), 'G': np.log(0.37), 'T': np.log(0.13)},
        'l': {'A': np.log(0.32), 'C': np.log(0.18), 'G': np.log(0.18), 'T': np.log(0.32)}
    }
    initial_probabilities = {'h': np.log(0.5), 'l': np.log(0.5)}
    s = read_fasta(
        "C:\\Users\\wzhyt\\Downloads\\dont_delete\\CS4775FinalProject\\CS4775FinalProject\\Data\\SampleDNAData")
    a, b = algorithms.viterbi(s, transition_probabilities, emission_probabilities, initial_probabilities)
    print(a)
    print(b)
