
def testingAround():
    print("test")

def viterbi(obs, trans_probs, emiss_probs, init_probs):
    ''' Complete this function. '''

    # m is a dictionary that has two keys, 'h' and 'l', the values of each key being a list with the length of the observations.
    # The value list for key 'h' should be storing the log-probabilities of the optimal Viterbi path to state 'h' up to observation x_t for each t.
    # The value list for key 'l' should be storing the log-probabilities of the optimal Viterbi path to state 'l' up to observation x_t for each t.
    # We provided the correct m in this format (in JSON format) for test case 1 & 2 to help you debug, but m will not be part of the grade.
    # which means that you don't have to use m as specified here if you find it more natural using other implementations.
    m = {
            'h' : [],
            'l' : [],
        }
    m['h'].append(init_probs['h'] + emiss_probs['h'][obs[0]])
    m['l'].append(init_probs['l'] + emiss_probs['l'][obs[0]])
    backtrack = { # 0 back to h, 1 back to l, 2 no pointer
        'h' : [],
        'l' : [],
    }
    for i in range (1, len(obs)):
        hh = m['h'][-1] + trans_probs['h']['h'] + emiss_probs['h'][obs[i]]
        lh = m['l'][-1] + trans_probs['l']['h'] + emiss_probs['h'][obs[i]]
        hl = m['h'][-1] + trans_probs['h']['l'] + emiss_probs['l'][obs[i]]
        ll = m['l'][-1] + trans_probs['l']['l'] + emiss_probs['l'][obs[i]]
        if hh > lh:
            m['h'].append(hh)
            backtrack['h'].append(0)
        else:
            m['h'].append(lh)
            backtrack['h'].append(1)
        if hl > ll:
            m['l'].append(hl)
            backtrack['l'].append(0)
        else:
            m['l'].append(ll)
            backtrack['l'].append(1)
    prob = max(m['h'][-1], m['l'][-1])
    seq = ""
    flip = 1
    if m['h'][-1] > m['l'][-1]:
        flip = 0
    idx = len(backtrack['h']) - 1
    while idx > -1:
        if flip == 0:
            seq += 'h'
            flip = backtrack['h'][idx]
            idx -= 1
        else:
            seq += 'l'
            flip = backtrack['l'][idx]
            idx -= 1
    if flip == 0:
        seq += 'h'
    else:
        seq += 'l'
    finalSeq = seq[::-1]
    return finalSeq, prob