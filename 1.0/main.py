# -*- coding:Utf-8 -*-
from numpy import *
from menu import *
import scipy.io as sio

set_printoptions(edgeitems=8, precision=3, linewidth=250)

def readParameters(): #READS THE CONFIG 
    #from text files
    """
    etav = open('conformation data/etav_intra_rot.txt','r')
    wav = open('conformation data/wav_intra_trans.txt','r')
    uav = open('conformation data/uav_inter_rot.txt','r')
    vav = open('conformation data/vav_inter_trans.txt','r')
    seq = open('conformation data/sequence.txt','r')
    sequence = str(seq.readline())
    nbbases = len(sequence) #number of bases
    w = zeros(12*nbbases-6).reshape(4*nbbases-2,3)
    intra = arange(0,nbbases*4,4)
    inter = arange(2,nbbases*4-2,4)
    for i in intra:
        w[i] = array(etav.readline().rsplit(), dtype=float)     #intra rotations
        w[i+1] = array(wav.readline().rsplit(), dtype=float)    #intra translations
    for j in inter:
        w[j] = array(uav.readline().rsplit(), dtype=float)      #inter rotations
        w[j+1] = array(vav.readline().rsplit(), dtype=float)    #inter translations
    wav.close(); etav.close(); uav.close(); vav.close(); seq.close
    """
    #from matlab file
    conformation = sio.loadmat('conformation data/parms_basejhb_TATA.mat')
    etav = matrix(conformation['etav']) #intra rotations
    wav = matrix(conformation['wav'])   #intra translations
    uav = matrix(conformation['uav'])   #inter rotations
    vav = matrix(conformation['vav'])   #inter translations
    seq = open('conformation data/sequence.txt','r')
    sequence = str(seq.readline())
    nbbases = len(sequence)  #number of bases
    w = zeros(12*nbbases-6).reshape(4*nbbases-2,3)
    for i in range(nbbases):
        w[4*i] = etav[i]     #intra rotations
        w[4*i+1] = wav[i]    #intra translations
    for i in range(nbbases-1):
        w[4*i+2] = uav[i]    #inter rotations
        w[4*i+3] = vav[i]    #inter translations
    return w, sequence
    
def calculPoints(w): #COMPUTES COORDINATES IN THE LAB FRAME
    from scipy.linalg import sqrtm
    from bases import cayley
    w,sequence = w
    nbbases = (w.shape[0]+2)/4
    e = array([[1,0,0],[0,1,0],[0,0,1]]) #lab frame, arbitrairy
    r = zeros((2*nbbases,3))   #base positions #2n bases,3 components
    q = zeros((nbbases,3))     #basepairs positions
    d = zeros((2*nbbases,3,3)) #base frames #2n bases,3 vectors,3 components
    g = zeros((nbbases,3,3))   #intra base frames
    h = zeros((nbbases-1,3,3)) #inter base frames
    
    zeta0 = w[3]  #arbitrarily chosen
    Theta0 = w[2] #arbitrarily chosen
    for j in range(3):                              #INITIALISATION q0,g0
        q[0] += zeta0[j]*e[j]
    g[0] = matrix(dot(e,cayley(Theta0)))
    for a in range(nbbases-1):                      #GENERATES FRAMES g
        L = cayley(w[4*a+2])
        g[a+1] = dot(g[a],L)
    for a in range(nbbases-1):                      #GENERATES FRAMES h
        Theta = w[4*a+2]
        Phi = 2*arctan(sqrt(dot(Theta,Theta))/2.)
        Thnorm = 2*tan(Phi/4.)*(1./sqrt(dot(Theta,Theta)))*Theta
        h[a] = dot(g[a],cayley(Thnorm))
    for a in range(nbbases-1):                      #GENERATES POINTS q
        zeta = w[4*a+3]
        q[a+1] = q[a] + dot(h[a],zeta)
    for a in range(nbbases):                        #GENERATES POINTS r
        xi = w[4*a+1]
        r[a]         = q[a] - 0.5*dot(g[a],xi)
        r[a+nbbases] = q[a] + 0.5*dot(g[a],xi) #complementary
    for a in range(nbbases):                        #GENERATES FRAMES d
        theta = w[4*a]
        Lambda = matrix(cayley(theta))
        #CI-DESSOUS A VERIFIER!!!!!
        d[a+nbbases] = matrix(dot(g[a].T, sqrtm(Lambda).T)) #transpose=inverse
        d[a] = dot(d[a+nbbases],Lambda) 
    return r, q, d, g, h, sequence

def plot3d(cf, vp1, vp2, precision, m1,m2,centerline,showbases,seq,mult,ghost): #DRAWS THE MOLECULE(S)
    print "Loading mlab....."
    from anim import stiffness, ACP
    from draw import data3d
    from enthought.mayavi import mlab
    r, q, d, g, h, sequence = cf
    nbbases = q.shape[0]
    #mlab.options.backend = 'envisage' #to use the full Mayavi interface
    fig = mlab.figure(1, bgcolor=(0, 0, 0), size=(1024,768)) #creates a new window
    mlab.clf() #clears the window
    
    if centerline >0 :
        #DRAWS CENTRAL AXIS
        #axis_nodes = mlab.points3d(q[:,0],q[:,1],q[:,2],scale_factor=0.4,color=(1,0,1))
        ghost_axis = mlab.plot3d(q[:,0],q[:,1],q[:,2],tube_radius=0.05, color=(1,1,1)) #will never move
        axis = mlab.plot3d(q[:,0],q[:,1],q[:,2],tube_radius=0.05, color=(1,0,1))
    
    #DRAWS LAB FRAME e
    ee = mlab.quiver3d([0,0,0],[0,0,0],[0,0,0],[1,0,0],[0,1,0],[0,0,1], \
                  mode='arrow',scale_factor=2, color=(1,1,1))

    if showbases >0 :
        #DRAWS BASES AND BASE FRAMES
        points, triangles, centres, frames, i,j, repere, s = data3d(r,d,sequence) #see draw.py
        if ghost >0 :
            ghost = mlab.triangular_mesh(points[:,0],points[:,1],points[:,2],triangles,color=(1,1,1),opacity=0.4)
        k = array(repere)[nbbases,0] #index i : points
        l = array(repere)[nbbases,1] #index j : triangles
        main_strand = mlab.triangular_mesh(points[:k,0],points[:k,1],points[:k,2],triangles[:l],\
                      scalars=s[:k], colormap="autumn") #main strand
        comp_strand = mlab.triangular_mesh(points[k:i,0],points[k:i,1],points[k:i,2],triangles[l:j]-triangles[l,0],\
                      scalars=s[k:i], colormap="autumn") #complementary
        baseframes  = mlab.quiver3d(centres[:,0],centres[:,1],centres[:,2],frames[:,0],frames[:,1],frames[:,2], \
                      mode='arrow',scale_factor=2, colormap="Oranges")

    #ANIMATION
    if vp1 >0 :
        if centerline >0 :
            #axis_nodes_src = axis_nodes.mlab_source
            axis_src = axis.mlab_source
        if showbases >0 :
            main_strand_src = main_strand.mlab_source
            comp_strand_src = comp_strand.mlab_source
            baseframes_src  = baseframes.mlab_source
        vecp = stiffness(m1,1,vp1-1) #chosen eigenvector 1
        #mult = 2.0 #multiple of the added eigenvector
        fig.scene.anti_aliasing_frames = 0 #antialiasing, default=8, off=0
        for p in range(precision):
            wmod = ACP(readParameters(),vecp,mult*(p+1)/precision) #modified configuration
            cfmod = calculPoints(wmod) #modified coordinates
            r_p, q_p, d_p, g_p, h_p, sequence = cfmod
            fig.scene.disable_render = True #accelerates rendering
            if showbases >0 : 
                points_p, triangles_p, centres_p, frames_p, i,j, repere, s = data3d(r_p,d_p,sequence) #modified graphics
                main_strand_src.set(x=points_p[:k,0],y=points_p[:k,1],z=points_p[:k,2])
                comp_strand_src.set(x=points_p[k:i,0],y=points_p[k:i,1],z=points_p[k:i,2])
                baseframes_src.reset(x=centres_p[:,0],y=centres_p[:,1],z=centres_p[:,2],\
                                   u=frames_p[:,0],v=frames_p[:,1],w=frames_p[:,2])
            if centerline >0 :
                #axis_nodes_src.set(x=q_p[:,0],y=q_p[:,1],z=q_p[:,2],scale_factor=0.7,color=(1,0,1))
                axis_src.set(x=q_p[:,0],y=q_p[:,1],z=q_p[:,2],tube_radius=0.05, color=(1,0,1))
            fig.scene.disable_render = False
    if vp2 > 0:
        if showbases >0 :
            #new instance of the entire molecule
            new_main_strand = mlab.triangular_mesh(points[:k,0],points[:k,1],points[:k,2],triangles[:l],\
                      scalars=s[:k], colormap="winter") #main strand
            new_comp_strand = mlab.triangular_mesh(points[k:i,0],points[k:i,1],points[k:i,2],triangles[l:j]-triangles[l,0],\
                      scalars=s[k:i], colormap="winter") #complementary
            new_baseframes  = mlab.quiver3d(centres[:,0],centres[:,1],centres[:,2],frames[:,0],frames[:,1],frames[:,2], \
                      mode='arrow',scale_factor=2, colormap="Blues")
            new_main_strand_src = new_main_strand.mlab_source
            new_comp_strand_src = new_comp_strand.mlab_source
            new_baseframes_src  = new_baseframes.mlab_source
        if centerline >0 :
            new_axis = mlab.plot3d(q[:,0],q[:,1],q[:,2],tube_radius=0.05, color=(1,0,1))
            new_axis_src        = new_axis.mlab_source
        vecp = stiffness(m2,1,vp2-1) #chosen eigenvector 2
        #mult = 2.0 #multiple of the added eigenvector
        
        for p in range(precision):
            wmod = ACP(readParameters(),vecp,mult*(p+1)/precision) #modified configuration
            cfmod = calculPoints(wmod) #modified coordinates
            r_p, q_p, d_p, g_p, h_p, sequence = cfmod
            if showbases >0 :
                points_p, triangles_p, centres_p, frames_p, i,j, repere, s = data3d(r_p,d_p,sequence) #modified graphics
                fig.scene.disable_render = True #accelerates rendering
                new_main_strand_src.set(x=points_p[:k,0],y=points_p[:k,1],z=points_p[:k,2])
                new_comp_strand_src.set(x=points_p[k:i,0],y=points_p[k:i,1],z=points_p[k:i,2])
                new_baseframes_src.reset(x=centres_p[:,0],y=centres_p[:,1],z=centres_p[:,2],\
                                   u=frames_p[:,0],v=frames_p[:,1],w=frames_p[:,2])
            if centerline >0 :
                #new_axis_nodes_src.set(x=q_p[:,0],y=q_p[:,1],z=q_p[:,2],scale_factor=0.7,color=(1,0,1))
                new_axis_src.set(x=q_p[:,0],y=q_p[:,1],z=q_p[:,2],tube_radius=0.05, color=(1,0,1))
            fig.scene.disable_render = False
    mlab.orientation_axes()
    mlab.show()


def run(): 
    menu(plot3d,calculPoints(readParameters()))
    
run()

