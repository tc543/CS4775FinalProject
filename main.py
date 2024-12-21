import math

import testing
import pandas as pd


import numpy as np
import matplotlib.pyplot as plt


def main():
    testing.testRunTime()
    # d = {
    #     'Viterbi': testing.viterbiTime,
    #     'Bellman': testing.bellmanTime,
    #     'Dijkstra': testing.dijkstraTime
    # }
    # df = pd.DataFrame(data=d)
    # d2 = {
    #     'Viterbi': [10 ** x for x in testing.viterbiTime],
    #     'Bellman': [10 ** x for x in testing.bellmanTime],
    #     'Dijkstra': [10 ** x for x in testing.dijkstraTime]
    # }
    # df2 = pd.DataFrame(data=d2)
    # d1 = {
    #     'Viterbi': testing.viterbiMem,
    #     'Bellman': testing.bellmanMem,
    #     'Dijkstra': testing.dijkstraMem
    # }
    # df1 = pd.DataFrame(data=d1)
    #
    #
    # xaxis = [i for i in range(len(testing.viterbiTime))]
    # barWidth = 0.25
    # br1 = np.arange(len(testing.viterbiTime))
    # br2 = [x + barWidth for x in br1]
    # br3 = [x + barWidth for x in br2]
    # plt.bar(br1, testing.viterbiTime, color ='r', width = barWidth,
    #         edgecolor ='grey', label ='Viterbi')
    # plt.bar(br2, testing.bellmanTime, color ='g', width = barWidth,
    #         edgecolor ='grey', label ='Bellman Ford')
    # plt.bar(br3, testing.dijkstraTime, color ='b', width = barWidth,
    #         edgecolor ='grey', label ='Dijkstra')
    # plt.xticks(xaxis)
    # plt.xlabel('DNA sequence', fontweight ='bold', fontsize = 15)
    # plt.ylabel('Time it takes to run in log10 scale (ms)', fontweight='bold', fontsize=15)
    # plt.legend()
    # # plt.savefig('Run_Time.png')
    #
    # plt.clf()
    # plt.bar(br1, testing.viterbiMem, color='r', width=barWidth,
    #         edgecolor='grey', label='Viterbi')
    # plt.bar(br2, testing.bellmanMem, color='g', width=barWidth,
    #         edgecolor='grey', label='Bellman Ford')
    # plt.bar(br3, testing.dijkstraMem, color='b', width=barWidth,
    #         edgecolor='grey', label='Dijkstra')
    # plt.xticks(xaxis)
    # plt.xlabel('DNA sequence', fontweight='bold', fontsize=15)
    # plt.ylabel('Peak Memory Usage (MB)', fontweight='bold', fontsize=15)
    # # plt.savefig('Memory_Usage.png')

    xaxis = [i for i in range(len(testing.realVTime))]
    barWidth = 0.25
    br1 = np.arange(len(testing.realVTime))
    br2 = [x + barWidth for x in br1]
    br3 = [x + barWidth for x in br2]
    plt.clf()
    plt.bar(br1, testing.realVTime, color='r', width=barWidth,
            edgecolor='grey', label = 'Viterbi')
    plt.bar(br2, testing.realBTime, color='g', width=barWidth,
            edgecolor='grey', label='Bellman Ford')
    plt.bar(br3, testing.realDTime, color='b', width=barWidth,
            edgecolor='grey', label='Dijktra')
    plt.xticks(xaxis)
    plt.xlabel('Real DNA sequence', fontweight='bold', fontsize=15)
    plt.ylabel('Time it takes to run (ms)', fontweight='bold', fontsize=15)
    plt.savefig('Real_Run_Time.png')

    plt.clf()
    plt.bar(br1, testing.realVMem, color='r', width=barWidth,
            edgecolor='grey', label='Viterbi')
    plt.bar(br2, testing.realBMem, color='g', width=barWidth,
            edgecolor='grey', label='Bellman Ford')
    plt.bar(br3, testing.realDMem, color='b', width=barWidth,
            edgecolor='grey', label='Dijktra')
    plt.xticks(xaxis)
    plt.xlabel('Real DNA sequence', fontweight='bold', fontsize=15)
    plt.ylabel('Peak Memory Usage (MB)', fontweight='bold', fontsize=15)
    plt.savefig('Real_Memory_Usage.png')

    # print(df)
    # print(df1)
    # print(df2)


if __name__ == "__main__":
    main()