import testing
import pandas as pd
def main():
    testing.testRunTime()
    d = {
        'Viterbi': testing.viterbiTime,
        'Bellman': testing.bellmanTime,
        'Dijkstra': testing.dijkstraTime
    }
    df = pd.DataFrame(data=d)
    print(df)

if __name__ == "__main__":
    main()