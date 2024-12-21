import numpy as np
import time
import math
import bellman
import bidijkstra
import dijkstra
import viterbi
import tracemalloc


# Read the file
def read_fasta(filename):
    with open(filename, "r") as f:
        s = ""
        s += f.readlines()[1].strip()
    return s


# Read the multilines file
def read_fasta_multilines(filename):
    with open(filename, "r") as f:
        s = ""
        skip = True
        for line in f.readlines():
            if skip:
                skip = False
                continue
            s += line.strip()
        realS = ""
        for c in s:
            if c == 'c':
                realS += "C"
            elif c == "a":
                realS += "A"
            elif c == "t":
                realS += "T"
            elif c == "g":
                realS += "G"
            else:
                realS += c
    return realS


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
bidijkstraTime = []
viterbiMem = []
bellmanMem = []
dijkstraMem = []
bidijkstraMem = []
realVTime = []
realBTime = []
realDTime = []
realBiDTime = []
realVMem = []
realBMem = []
realDMem = []
realBiDMem = []


# Running all DNA seq on viterbi
def generateViterbiTime():
    for file in files:
        start_time = time.time_ns()
        a, b = viterbi.viterbi(file, transition_probabilities, emission_probabilities, initial_probabilities)
        end_time = time.time_ns()
        time_elapsed = (end_time - start_time) / (1000 * 1000)
        viterbiTime.append(math.log10(time_elapsed))
        tracemalloc.start()
        a, b = viterbi.viterbi(file, transition_probabilities, emission_probabilities, initial_probabilities)
        viterbiMem.append(tracemalloc.get_traced_memory()[1] / (1024 * 1024))
        tracemalloc.stop()


# Running all DNA seq on bellman ford
def generateBellmanTime():
    for file in files:
        start_time = time.time_ns()
        a, b = bellman.bellman_ford(file, transition_probabilities, emission_probabilities, initial_probabilities)
        end_time = time.time_ns()
        time_elapsed = (end_time - start_time) / 1000000
        bellmanTime.append(math.log10(time_elapsed))
        tracemalloc.start()
        a, b = bellman.bellman_ford(file, transition_probabilities, emission_probabilities, initial_probabilities)
        bellmanMem.append(tracemalloc.get_traced_memory()[1] / (1024 * 1024))
        tracemalloc.stop()


# Running all DNA seq on dijkstra
def generateDijkstraTime():
    for file in files:
        start_time = time.time_ns()
        a, b = dijkstra.dijkstra(file, transition_probabilities, emission_probabilities, initial_probabilities)
        end_time = time.time_ns()
        time_elapsed = (end_time - start_time) / 1000000
        dijkstraTime.append(math.log10(time_elapsed))
        tracemalloc.start()
        a, b = dijkstra.dijkstra(file, transition_probabilities, emission_probabilities, initial_probabilities)
        dijkstraMem.append(tracemalloc.get_traced_memory()[1] / (1024 * 1024))
        tracemalloc.stop()


# Running all DNA seq on bidijkstra
def generateBidijkstraTime():
    for file in files:
        start_time = time.time_ns()
        bidijkstra.BiDirectionalDijkstra()
        a, b = bidijkstra.BiDirectionalDijkstra.dijkstra(file, transition_probabilities, emission_probabilities, initial_probabilities)
        end_time = time.time_ns()
        time_elapsed = (end_time - start_time) / 1000000
        dijkstraTime.append(math.log10(time_elapsed))
        tracemalloc.start()
        a, b = dijkstra.dijkstra(file, transition_probabilities, emission_probabilities, initial_probabilities)
        dijkstraMem.append(tracemalloc.get_traced_memory()[1] / (1024 * 1024))
        tracemalloc.stop()


# Running the two real DNA seq on all algorithms
def generateRealSeqTime():
    realFile = []
    realFile.append(read_fasta_multilines("C:\\Users\\wzhyt\\Downloads\\dont_delete\\CS4775FinalProject\\CS4775FinalProject\\chimpanzee.txt"))
    realFile.append(read_fasta_multilines("C:\\Users\\wzhyt\\Downloads\\dont_delete\\CS4775FinalProject\\CS4775FinalProject\\human_chrome.txt"))
    for file in realFile:
        start_time = time.time_ns()
        a, b = viterbi.viterbi(file, transition_probabilities, emission_probabilities, initial_probabilities)
        end_time = time.time_ns()
        time_elapsed = (end_time - start_time) / 1000000
        realVTime.append(math.log10(time_elapsed))
        tracemalloc.start()
        a, b = viterbi.viterbi(file, transition_probabilities, emission_probabilities, initial_probabilities)
        realVMem.append(tracemalloc.get_traced_memory()[1] / (1024 * 1024))
        tracemalloc.stop()
        print(time_elapsed)
        # start_time = time.time_ns()
        # a, b = bellman.bellman_ford(file, transition_probabilities, emission_probabilities, initial_probabilities)
        # end_time = time.time_ns()
        # time_elapsed = (end_time - start_time) / 1000000
        # realBTime.append(math.log10(time_elapsed))
        # tracemalloc.start()
        # a, b = bellman.bellman_ford(file, transition_probabilities, emission_probabilities, initial_probabilities)
        # realBMem.append(tracemalloc.get_traced_memory()[1] / (1024 * 1024))
        # tracemalloc.stop()
        # # print(time_elapsed)
        start_time = time.time_ns()
        a, b = dijkstra.dijkstra(file, transition_probabilities, emission_probabilities, initial_probabilities)
        end_time = time.time_ns()
        time_elapsed = (end_time - start_time) / 1000000
        realDTime.append(math.log10(time_elapsed))
        tracemalloc.start()
        a, b = dijkstra.dijkstra(file, transition_probabilities, emission_probabilities, initial_probabilities)
        realDMem.append(tracemalloc.get_traced_memory()[1] / (1024 * 1024))
        tracemalloc.stop()
        print(time_elapsed)


# Running all the algorithms
def testRunTime():
    # generateViterbiTime()
    # generateBellmanTime()
    # generateDijkstraTime()
    generateBidijkstraTime()
    # generateRealSeqTime()
