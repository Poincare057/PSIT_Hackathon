from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

f = '1.nc'

file = Dataset(f)

file.variables.keys()

lat = np.asarray(file.variables["latitude"])
long = np.asarray(file.variables["longitude"])
time = np.asarray(file.variables["time"])

def sup(M, T, X, Y):
    res = 0
    for t in range(T):
        for x in range(X):
            for y in range(Y):
                if M[t][x][y] > res:
                    res = M[t][x][y]
    return(res)

def supT(M, T):
    res = 0
    for t in range(T):
        if M[t] > res:
            res = M[t]
    return(res)

T = 2184
X = 5
Y = 5
t2m = np.asarray(file.variables["t2m"])
t2m = t2m/sup(t2m, T, X, Y)
t2m = t2m - np.mean(t2m)
tcc = np.asarray(file.variables["tcc"])
tcc = tcc/sup(tcc, T, X, Y)
ssrd = np.asarray(file.variables["ssrd"])
ssrd = ssrd/sup(ssrd, T, X, Y)

def pssrd(x, y):
    res = []
    for t in range(2184):
        res += [ssrd[t][x][y]]
    return(np.array(res))

def ptcc(x, y):
    res = []
    for t in range(2184):
        res += [tcc[t][x][y]]
    return(np.array(res))

def pt2m(x, y):
    res = []
    for t in range(2184):
        res += [t2m[t][x][y]]
    return(np.array(res))

def plot_at(x, y):
    plt.plot(pssrd(x, y))
    plt.plot(ptcc(x, y))
    plt.plot(pt2m(x, y))
    plt.legend(["ssrd", "tcc", "t2m"])
    plt.xlabel("Location is (" + str(x) + ", " + str(y) + ")")
    plt.show()
    
def phase_find0(x, y):
    res = 1000000
    soln = 0
    for t in range(T):
        k = np.mean(abs(pssrd(x, y) - pt2m(x, y) - np.array(list(ptcc(x, y)[t:]) + list(np.zeros(t)))))
        if k < res:
            res = k
            soln = t
    return(soln)

def KS_ssrd_ptcc_t2m():
    corrs = np.zeros((5, 5))
    for i in range(5):
        for j in range(5):
            corrs[i][j] = supT(pssrd(i, j) - pt2m(i, j) - ptcc(i, j), T)
    return corrs 

def plot_KS():
    plt.contourf(np.arange(0, 5, 1), np.arange(0, 5, 1), KS_ssrd_ptcc_t2m())
    plt.show()

g = "2.nc"

file = Dataset(g)

file.variables.keys()

lat1 = np.asarray(file.variables["latitude"])
long1 = np.asarray(file.variables["longitude"])
time1 = np.asarray(file.variables["time"])

T = 252
X = 141
Y = 141
t2m1 = np.asarray(file.variables["t2m"])
t2m1 = t2m1/sup(t2m1, T, X, Y)
t2m1 = t2m1 - np.mean(t2m1)
tcc1 = np.asarray(file.variables["tcc"])
tcc1 = tcc1/sup(tcc1, T, X, Y)
ssrd1 = np.asarray(file.variables["ssrd"])
ssrd1 = ssrd1/sup(ssrd1, T, X, Y)

def pssrd1(x, y):
    res = []
    for t in range(252):
        res += [ssrd1[t][x][y]]
    return(np.array(res))

def ptcc1(x, y):
    res = []
    for t in range(252):
        res += [tcc1[t][x][y]]
    return(np.array(res))

def pt2m1(x, y):
    res = []
    for t in range(252):
        res += [t2m1[t][x][y]]
    return(np.array(res))

def plot_at1(x, y):
    plt.plot(pssrd1(x, y))
    plt.plot(ptcc1(x, y))
    plt.plot(pt2m1(x, y))
    plt.legend(["ssrd", "tcc", "t2m"])
    plt.xlabel("Location is (" + str(x) + ", " + str(y) + ")")
    plt.show()

def phase_find(x, y):
    res = 1000000
    soln = 0
    for t in range(252):
        k = np.mean(abs(pssrd1(x, y) - pt2m1(x, y) - np.array(list(ptcc1(x, y)[t:]) + list(np.zeros(t)))))
        if k < res:
            res = k
            soln = t
    return(soln)
        
def plot_corr_1(x, y):
    plt.plot(pssrd1(x, y) - pt2m1(x, y))
    plt.plot(ptcc1(x, y)[phase_find(x, y):] + 0.2)
    plt.xlabel("Location is (" + str(x) + ", " + str(y) + ")")
    plt.legend(["ssrd - pt2m", "ptcc"])
    plt.show()

def spacorr_ssrd_ptcc1():
    corrs = np.zeros((X, Y))
    for i in range(X):
        for j in range(Y):
            corrs[i][j] = supT(pssrd1(i, j) - ptcc1(i, j), T)
    return corrs

def spacorr_ssrd_pt2m1():
    corrs = np.zeros((X, Y))
    for i in range(X):
        for j in range(Y):
            corrs[i][j] = supT(pssrd1(i, j) - pt2m1(i, j), T)
    return corrs

def KS_ssrd_ptcc_t2m1():
    corrs = np.zeros((X, Y))
    for i in range(X):
        for j in range(Y):
            corrs[i][j] = supT(pssrd1(i, j) - pt2m1(i, j) - ptcc1(i, j), T)
    return corrs 

def plot_KS_1():
    plt.contourf(np.arange(0, 141, 1), np.arange(0, 141, 1), KS_ssrd_ptcc_t2m1())
    plt.show()
