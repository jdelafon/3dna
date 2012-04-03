# -*- coding:utf-8 -*-
from numpy import *

def cayley(v):
    v = array(v)
    m = matrix(zeros((3,3)))
    I = matrix([[1,0,0],[0,1,0],[0,0,1]])
    m[0,1] = -v[2]; m[0,2] =  v[1]
    m[1,0] =  v[2]; m[1,2] = -v[0]
    m[2,0] = -v[1]; m[2,1] =  v[0]
    Lambda = I + (4./(4.+dot(v,v)))*(m+0.5*(m*m))
    return Lambda

def rectangle(r, dir1, dir2):
    long = 4. ; larg = 3.
    decalage = 2*long*dir2     
    c1 = r + long*dir2 + larg*dir1 +decalage
    c2 = r + long*dir2 - larg*dir1 +decalage
    c3 = r - long*dir2 - larg*dir1 +decalage
    c4 = r - long*dir2 + larg*dir1 +decalage
    rect = array([c1,c2,c3,c4,c1])
    return rect

def adenine(r, d): #purine
    ad = array(
    [[-1.291, 4.498, 0.000], #N9
     [ 0.024, 4.897, 0.000], #C8
     [ 0.877, 3.902, 0.000], #N7
     [ 0.071, 2.771, 0.000], #C5
     [ 0.369, 1.398, 0.000], #C6
     [-0.668, 0.532, 0.000], #N1
     [-1.912, 1.023, 0.000], #C2
     [-2.320, 2.290, 0.000], #N3
     [-1.267, 3.124, 0.000], #C4
     [-1.291, 4.498, 0.000]])#N9
    #creation des vecteurs de position dans la base d[a]
    vad = zeros((9,3))
    coord = zeros((10,3))
    coord[0] = ad[0]
    for n in range(9):
        vad[n] = ad[n+1] - ad[n]
    for n in range(9):
        coord[n+1] = coord[n]+vad[n]
    coord = dot(coord,d)+r #changement de base
    return coord #10 coordonnees, C8->C8
    
def cytosine(r, d): #pyrimidine
    cy = array(
    [[-1.285, 4.542, 0.000], #N1
     [-1.472, 3.158, 0.000], #C2
     [-0.391, 2.344, 0.000], #N3
     [ 0.837, 2.868, 0.000], #C4
     [ 1.056, 4.275, 0.000], #C5
     [-0.023, 5.068, 0.000], #C6
     [-1.285, 4.542, 0.000]])#N1
    #creation des vecteurs de position dans la base d[a]
    vcy = zeros((6,3))
    coord = zeros((7,3))
    coord[0] = cy[0]
    for n in range(6):
        vcy[n] = cy[n+1] - cy[n]
    for n in range(6):
        coord[n+1] = coord[n]+vcy[n]
    coord = dot(coord,d)+r  
    return coord #7 coordonnees, C6->C6
