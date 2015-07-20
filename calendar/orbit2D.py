#!/usr/bin/env python

from    numpy                 import    *
import  numpy                 as        np
import  matplotlib.pyplot     as        mpl
import  math
import  pylab
#import datetime,time
#from time import sleep
#import sys
from mpl_toolkits.mplot3d import axes3d
import planets

pi=math.pi
M=1000
AU=planets.AU

#############################################
# matrice matrix rotation : random around vector [a,b,c]
def rotationaxe(x,y,z,angle,a,b,c):

    R1=(1-cos(angle))*array([[a*a,a*b,a*c],[b*a,b*b,b*c],[c*a,c*b,c*c]])
    R2=cos(angle)*array([[1,0,0],[0,1,0],[0,0,1]])
    R3=sin(angle)*array([[0,-c,b],[c,0,-a],[-b,a,0]])
    R=R1+R2+R3
    x1=R[0,0]*x+R[0,1]*y+R[0,2]*z
    y1=R[1,0]*x+R[1,1]*y+R[1,2]*z
    z1=R[2,0]*x+R[2,1]*y+R[2,2]*z
    return x1,y1,z1

# matrice matrix rotation : around vector x,y or z
def rotation(x,y,z,axis,angle):
    i=angle
    if axis==1:
           print 'Rx, with angle: ',i*180./pi, ' deg'
           R=array([[1,0,0],[0,cos(i),-sin(i)],[0,sin(i),cos(i)]])
    elif axis==2:
           print 'Ry, with angle: ',i*180./pi, ' deg'
           R=array([[cos(i),0,sin(i)],[0,1,0],[-sin(i),0,cos(i)]])
    elif axis==3:
           print 'Rz, with angle: ',i*180./pi, ' deg'
           R=array([[cos(i),-sin(i),0],[sin(i),cos(i),0],[0,0,1]])
    else:
           print 'problem: choice of rotation axis out of bounds :'
           print '1: X, 2: Y, 3: Z'
    x1=R[0,0]*x+R[0,1]*y+R[0,2]*z
    y1=R[1,0]*x+R[1,1]*y+R[1,2]*z
    z1=R[2,0]*x+R[2,1]*y+R[2,2]*z
    return x1,y1,z1

###########################################

def getplan(dim):
        # plan horizontal ecliptique
        X, Y = np.meshgrid([-dim, dim], [-dim, dim])
        Z = np.zeros((2,2))
        return X,Y,Z

def getplanellipse(dim,incl,ascend):
        
        # plan de orbite --en cours
        X2, Y2 = np.meshgrid([-dim, dim], [-dim, dim])
        X3, Y3 = np.meshgrid([-dim, dim], [-dim, dim])

        angle=ascend*pi/180.
        Z2 = X2 * incl
        Z3 = X3 * incl
        #X3=X3*cos(angle)
        #Y3=Y3*sin(angle)
        print 'yeah'
        print X3,Y3,Z3

	# plan ellipse divise en moitie pour effet couleurs
        # moitie 1
#        X2, Y2 = np.meshgrid([-dim, dim], [0, dim])
#        Z2 = X2 * incl
        # moitie 2
#        X3, Y3 = np.meshgrid([-dim, dim], [-dim, 0])
#        Z3 = X3 * incl

        return X2,X3,Y2,Y3,Z2,Z3      
          
###########################################

def getellipse(ecc,a,ascend,periaps):

        b=a*pow(1-ecc*ecc,0.5)
        c=ecc*a
        apohelie=c+a
        perihelie=2*a-apohelie
        print "apohelie,perihelie (AU) =",apohelie,' / ',perihelie
             
        # resolution ellipses
        th = np.linspace(0, 2 * pi, M)
        x0= a * np.cos(th)+c 
        y0= b * np.sin(th) 
        z0= 0
        # rotation ascending node autour de Z0
        x1,y1,z1=rotationaxe(x0,y0,z0,ascend,0,0,1)
	# rotation angle omega periapsis autour de vx,vy,vz
        x2,y2,z2=rotationaxe(x1,y1,z1,periaps,0,0,1)

        return x2,y2

###########################################

def plotellipse(x,y,color,ax):

        lw=1  # linewidth : epaisseur du trait
        linestyle='-'
        color=color
        zorder=0
        #ax.plot(x[y < 0.], y[y < 0.], z[y < 0.], lw=lw, linestyle=linestyle, color=color, zorder=zorder)
        ax.plot(x, y, lw=lw, linestyle=linestyle, color=color, zorder=zorder)
        
###########################################

def plotplanellipse(X2,X3,Y2,Y3,Z2,Z3,ax):

        # plan ellipse moitie 1
        color='blue'
        alpha=.1  #transparence
        lw=0.  #contour du plan ( = 0 : pas de contour)
        zorder=-1 #property to tell matplotlib in what order the surfaces should be painted
        ax.plot_surface(X2, Y3 ,Z3, color=color, alpha=alpha, linewidth=lw, zorder=zorder)

        # plan ellipse moitie 2
        color='blue'
        alpha=.1  #transparence
        lw=0.  #contour du plan ( = 0 : pas de contour)
        zorder=3 #property to tell matplotlib in what order the surfaces should be painted
        ax.plot_surface(X2, Y2 ,Z2, color=color, alpha=alpha, linewidth=lw, zorder=zorder)
        
###########################################

def plotplaneclipt(X,Y,Z,ax):

        # plan ecliptic
        color='red'
        alpha=.1  #transparence
        lw=0.  #contour du plan ( = 0 : pas de contour)
        zorder=1 #property to tell matplotlib in what order the surfaces should be painted
        ax.plot_surface(X, Y ,Z, color=color, alpha=alpha, linewidth=lw, zorder=zorder)
        
###########################################

def plotstar(ax):
	# Star position 
        marker='o'
        color='k'
        ax.plot([0], [0], marker=marker,color=color)

###########################################

def plotplanet(x,y,ind,color,ax):
	# planet position 
        marker='o'
        ax.plot([x[ind]],[y[ind]], marker=marker,color=color)

###########################################

def getplanet(myplanet,color):

    #parametres ellipse
    ecc=myplanet.eccentricity
    a=myplanet.rsm/AU  #sma
    #incl=myplanet.incl
    ascend=myplanet.ascend #110.3 # ascending angle (deg)
    periaps=myplanet.omeg # periapsis angle (deg)

    dim=2*a
    #incl=incl*pi/180.
    ascend=ascend*pi/180.
    periaps=periaps*pi/180.
 
    x,y=getellipse(ecc,a,ascend,periaps)
    #X2,X3,Y2,Y3,Z2,Z3=getplanellipse(dim,incl,ascend)
    
    #plotplanet(x,y,ind,color)
    #plotplanellipse(X2,X3,Y2,Y3,Z2,Z3)
    return x,y,color,dim
###########################################

def main(nbplanet,planetes,col):
    fig = mpl.figure(figsize=(10,10))
    ax=mpl.gca()
    plotstar(ax)
    Ls=0
    ind=Ls*M/360.
    for i in range(nbplanet):
        x,y,color,dim=getplanet(planetes[i],col[i])
        plotplanet(x,y,ind,color,ax)
        plotellipse(x,y,color,ax)

    ax.set_xlim(-dim, dim)
    ax.set_ylim(-dim, dim)
    ax.set_xlabel('X Label (AU)')
    ax.set_ylabel('Y Label (AU)')

    mpl.axis('on')
    mpl.grid()
    mpl.show()


#col=['black','red','orange','green','blue','magenta','cyan','yellow']
#planetes=[planets.Venus,planets.Earth,planets.Mars]
#nbplanet=3
#main(nbplanet,planetes,col)













