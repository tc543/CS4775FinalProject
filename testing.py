import numpy as np
import time
import bellman
import dijkstra
import viterbi



# Read the file
def read_fasta(filename):
    with open(filename, "r") as f:
        s = ""
        s += f.readlines()[1].strip()
    return s


# List of all probabilities
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


# All the files
files = []
files.append(read_fasta("C:\\Users\\wzhyt\\Downloads\\dont_delete\\CS4775FinalProject\\CS4775FinalProject\\Data\\simulatedsequence1.fasta"))
files.append(read_fasta("C:\\Users\\wzhyt\\Downloads\\dont_delete\\CS4775FinalProject\\CS4775FinalProject\\Data\\simulatedsequence2.fasta"))
files.append(read_fasta("C:\\Users\\wzhyt\\Downloads\\dont_delete\\CS4775FinalProject\\CS4775FinalProject\\Data\\simulatedsequence3.fasta"))
files.append(read_fasta("C:\\Users\\wzhyt\\Downloads\\dont_delete\\CS4775FinalProject\\CS4775FinalProject\\Data\\simulatedsequence4.fasta"))
files.append(read_fasta("C:\\Users\\wzhyt\\Downloads\\dont_delete\\CS4775FinalProject\\CS4775FinalProject\\Data\\simulatedsequence5.fasta"))


# All the data of times
viterbiTime = []
bellmanTime = []
dijkstraTime = []


# Running all DNA seq on viterbi
def generateViterbiTime():
    for file in files:
        start_time = time.time_ns()
        a, b = viterbi.viterbi(file, transition_probabilities, emission_probabilities, initial_probabilities)
        end_time = time.time_ns()
        time_elapsed = (end_time - start_time) / 1000000
        viterbiTime.append(time_elapsed)


# Running all DNA seq on bellman ford
def generateBellmanTime():
    for file in files:
        start_time = time.time_ns()
        a, b = bellman.bellman_ford(file, transition_probabilities, emission_probabilities, initial_probabilities)
        end_time = time.time_ns()
        time_elapsed = (end_time - start_time) / 1000000
        bellmanTime.append(time_elapsed)


# Running all DNA seq on dijkstra
def generateDijkstraTime():
    for file in files:
        start_time = time.time_ns()
        a, b = dijkstra.dijkstra(file, transition_probabilities, emission_probabilities, initial_probabilities)
        end_time = time.time_ns()
        time_elapsed = (end_time - start_time) / 1000000
        dijkstraTime.append(time_elapsed)


# Running all the algorithms
def testRunTime():
    generateViterbiTime()
    generateBellmanTime()
    generateDijkstraTime()
