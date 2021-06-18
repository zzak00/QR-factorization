#!/usr/bin/env python
# coding: utf-8
# Author: zzak00

import numpy as np

def Hv(x):
    n=len(x)
    e1=np.zeros(n)
    e1[0]=1
    v=x+x[0]/abs(x[0])*np.linalg.norm(x)*e1
    H=np.eye(n,n)-2*np.outer(v,v)/v.dot(v)
    return H


def HouseHolder_decomp(A):
    m,n=A.shape
    R=A
    Q=np.eye(m,m)
    H=np.eye(m,m)
    for i in range(n):
        H[i:,i:]=Hv(R[i:,i])
        print(H)
        R=H.dot(R)
        Q=Q.dot(H)
        H=np.eye(m,m)
    return Q,R   



def givens_rotation(A):
    m,n = np.shape(A)

    Q = np.eye(m)
    R = np.copy(A)

    (rows, cols) = np.tril_indices(m, -1, n)
    for (row, col) in zip(rows, cols):
        if R[row, col] != 0:
            (c, s) = gmre(R[col, col], R[row, col])
            G = np.eye(m)
            G[col, col] = c
            G[row, row] = c
            G[row, col] = s
            G[col, row] = -s
            R = np.dot(G, R)
            Q = np.dot(Q, G.T)
    return (Q, R)


def gmre(a, b):
    r = np.sqrt(a**2+b**2)
    c = a/r
    s = -b/r
    return (c, s)

def Gram_Schmidt_process(A):
    m,n=A.shape
    R=np.zeros([m,n])
    Q=np.zeros([m,m])
    R[0,0]=np.linalg.norm(A[:,0])
    Q[:,0]=A[:,0]/R[0,0]
    for j in range(1,n):
        qt=A[:,j]
        for i in range(j):
            R[i,j]=qt.dot(Q[:,i])
            qt=qt-R[i,j]*Q[:,i]
        R[j,j]=np.linalg.norm(qt)
        if R[j,j]==0:
            break
        else:
            Q[:,j]=qt/R[j,j]
    R=R[~np.all(R == 0, axis=1)]
    Q=Q.T
    Q=Q[~np.all(Q == 0, axis=1)].T
    return Q,R
    



