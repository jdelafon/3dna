# -*- coding:utf-8 -*-
from numpy import *
from bases import adenine, cytosine

def data3d(r,d,sequence): #BUILDS THE VTK SOURCES AND COMPUTES THE BASES' SHAPES
    nbbases = r.shape[0]/2
    points = zeros((2*nbbases*10,3))
    triangles = zeros((2*nbbases*7,3))
    centres = zeros((2*nbbases,3,3))
    frames = zeros((2*nbbases,3,3))
    i=0; j=0 #indexes of the arrays points and triangles
    repere=[[0,0]] #indexe i et j pour les retrouver ensuite
    scalars=[]

    #REFERENCE STRAND
    for a in range(nbbases):
        #DESSINE LES BASES
        if (sequence[a] == "A") | (sequence[a] =="G"): #purine
            base = adenine(r[a], d[a])
            c = 0.5*(base[3]+base[6]) #arbitraire
            points[i:i+10] = base
            triangles[j:j+7] = \
                array(((i,i+1,i+2),(i,i+2,i+8),(i+2,i+3,i+8),\
                (i+3,i+4,i+8),(i+4,i+7,i+8),(i+4,i+6,i+7),\
                (i+4,i+5,i+6)))
            i+=10; j+=7
            repere.append([i,j])
        elif (sequence[a] == "T") | (sequence[a] =="C"): #pyrimidine
            base = cytosine(r[a], d[a])
            c = 0.5*(base[0]+base[3])
            points[i:i+7] = base
            triangles[j:j+4] = \
                array(((i,i+1,i+2),(i,i+2,i+5),(i+2,i+3,i+5),(i+3,i+4,i+5)))
            i+=7; j+=4
            repere.append([i,j])
        else : print "sequence incorrecte"; break
        #DEFINES SCALARS FOR THE COLORS
        if (sequence[a] == "A"):
            for x in range(10): scalars.append(1)
        elif (sequence[a] =="G"):
            for x in range(10): scalars.append(2)
        elif (sequence[a] =="T"):
            for x in range(7): scalars.append(3)
        elif (sequence[a] =="C"):
            for x in range(7): scalars.append(4)
        #DESSINE LES REPERES d
        centres[a] = vstack((c,c,c)).T
        frames[a] = d[a].T

    #COMPLEMENTARY STRAND
    for a in range(nbbases):
        #DRAWS FRAMES
        if (sequence[a] == "T") | (sequence[a] =="C"):
            #aussi inverser en x ???
            base = adenine(r[a+nbbases], [d[a+nbbases,0],-d[a+nbbases,1],d[a+nbbases,2]])
            c = 0.5*(base[3]+base[6])
            points[i:i+10] = base
            triangles[j:j+7] = \
                array(((i,i+1,i+2),(i,i+2,i+8),(i+2,i+3,i+8),\
                (i+3,i+4,i+8),(i+4,i+7,i+8),(i+4,i+6,i+7),\
                (i+4,i+5,i+6)))
            i+=10; j+=7
            repere.append([i,j])
        elif (sequence[a] == "A") | (sequence[a] =="G"):
            base = cytosine(r[a+nbbases], [d[a+nbbases,0],-d[a+nbbases,1],d[a+nbbases,2]])
            c = 0.5*(base[0]+base[3])
            points[i:i+7] = base
            triangles[j:j+4] = \
                array(((i,i+1,i+2),(i,i+2,i+5),\
                (i+2,i+3,i+5),(i+3,i+4,i+5)))
            i+=7; j+=4
            repere.append([i,j])
        else : print "sequence incorrecte"; break
        #DEFINES SCALARS FOR THE COLORS
        if (sequence[a] == "A"):
            for x in range(7): scalars.append(3)
        elif (sequence[a] =="G"):
            for x in range(7): scalars.append(4)
        elif (sequence[a] =="T"):
            for x in range(10): scalars.append(1)
        elif (sequence[a] =="C"):
            for x in range(10): scalars.append(2)
        #DRAWS FRAMES d
        centres[a+nbbases] = vstack((c,c,c)).T
        frames[a+nbbases] = d[a+nbbases].T

    return points, triangles, centres, frames, i,j, repere, scalars




