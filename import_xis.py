import numpy as np

def load_xis(datapath):
    f = open(datapath + 'xis.txt', "r")
    liste = []

    for line in f:
        liste.append(line.split())
    f.close()

    xis = np.zeros(len(liste))
    N_xi = len(liste)
    for i in range(N_xi):
        xis[i] = liste[i][0]
    
    return xis