# -*- coding:utf-8 -*-
import scipy.io as sio
from numpy import *
from scipy.linalg import eig,eigh

set_printoptions(precision=5, suppress=True)
#precision: number of digits after floating point
#suppress: whether or not suppress printing of small floating point values using scientific notation

def stiffness(m=0,trunc=1,vp_id=0): #RETURNS EIGENVECTOR OF MATRIX m CORRESPONDING TO THE vp_id'th EIGENVALUE
    parms = sio.loadmat('conformation data/stiffs.mat')
    if m==0 : #stiffness matrix
        #parms = sio.loadmat('conformation data/parms_basejhb_TATA.mat')
        stiff = matrix(parms['scTATA'])
    elif m==1 : #truncated stiffness
        stiff = matrix(parms['tr66TATA'])
    elif m==2 :
        stiff = matrix(parms['tr55TATA'])
    elif m==3 :
        stiff = matrix(parms['tr53TATA'])
    stiff = around(stiff, 9) #to ensure it's symmetric
    print all(stiff-stiff.T == 0) #check if symmetric
    nbbases = stiff.shape[0]
    vp, vecp = eig(stiff) #eigenvalues, eigenvectors ###EIGH ???????????????????????
    ind = argsort(vp)
    vp = vp[ind]
    vecp = vecp[:,ind]
    #eigens = dict(zip(vp,vecp)) #Creates a dictionary with keys from vp and values set to vecp
    vecp_chosen = vecp[:,vp_id] #eigenvectors are columns...
    #print "Eigenvector:\n",vecp_chosen
    #print "Rotations:\n", vecp_chosen.reshape(nbbases/6,6)[:,0:3]
    #print "Translations:\n", vecp_chosen.reshape(nbbases/6,6)[:,3:6]
    #print "Eigenvalues:\n", vals
    #print "Smallest eigenvalue: ", vals[vals.shape[0]-1]
    #print "Largest eigenvalue: ", vals[0]
    print "Eigenvalue: ", vp[vp_id]
    neg=0
    for i in range(vp.shape[0]):
        if vp[i]<0 : neg += 1
    print "Number of negative eigenvalues", neg
    return vecp_chosen

def ACP(w,vecp,mult=1): #RETURNS THE CONFIG w MODIFIED BY AN (EIGEN)VECTOR vecp
    w, sequence = w
    shape = w.shape
    w = w.ravel() #flattens
    wmod = w + mult*vecp
    wmod = wmod.reshape(shape)
    return wmod, sequence
