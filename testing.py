import algorithms
import numpy as np
import time

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
    start_time = time.time_ns()
    a, b = algorithms.viterbi(s, transition_probabilities, emission_probabilities, initial_probabilities)
    end_time = time.time_ns()
    time_elapsed = (end_time - start_time) / 1000000
    # print(a)
    # print(b)
    print("Time it takes to run test1", time_elapsed, "seconds")
