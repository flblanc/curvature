
#!/usr/bin/env python
#-*- coding: utf-8-*-

# Version  1.0
# florian.blanc@unistra.fr

#  Copyright 2015 Florian Blanc <flblanc@vador>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

# 24.05.2015 00:00:37 CEST

# Linux vador 3.13.0-49-generic #83-Ubuntu SMP Fri Apr 10 20:11:33 UTC 2015 x86_64 x86_64 x86_64 GNU/Linux

# Importing the required modules...

# Authorship information in Python style:

__author__ = "Florian Blanc"
__copyright__="Copyright 2015, Florian Blanc"
__credits__="Florian Blanc"

__license__="GPL"
__email__="florian.blanc@unistra.fr"
__date__= "24.05.2015 00:00:37 CEST"

import numpy as np
import numpy.linalg

import matplotlib.pyplot as pl

# 3D plotting related stuff
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

def TargetFunction(p,x,y):
    u"""
    The general form of the quadratic surface. 
    p = array of parameters
    """
    
    return p[0] + p[1]*x + p[2]*(x**2) + p[3]*y + p[4]*(y**2) + p[5]*x*y
    
#p_real = np.ones((6)) # The TRUE parameters of the quadratic surface
p_real = 2.*np.random.random(6) - 1.0
p = p_real
p[0] = 0.0

#~ p[1] = 1.
#~ p[2] = 1.
#~ p[3] = 1.
#~ p[4] = 2.
#~ p[5] = 1.

print p 


# Point set generation; first I build a grid:

N = 30 # number of generated points

xx = np.linspace(-10.,10.,N)
yy = np.linspace(-10.,10.,N)

X, Y = np.meshgrid(xx,yy)

# Then I compute the actual surface

Z = TargetFunction(p, X,Y)

#~ Z -= (Z.min() - 1.) # offset

pl.ion()
#fig = pl.figure(figsize=pl.figaspect(0.5)*1.5)
fig = pl.figure()
ax = fig.gca(projection='3d')
#~ ax.set_aspect('equal')

#~ surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.bwr,#color='cyan',#cmap=cm.coolwarm,
        #~ linewidth=0.0, antialiased=False,alpha=0.6)
        
        
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, color='cyan',#cmap=cm.coolwarm,
        linewidth=0.0, antialiased=False,alpha=0.6)


cset = ax.contour(X, Y, Z, zdir='z', offset=-100., cmap=cm.coolwarm,linewidth=6.)
#~ cset = ax.contour(X, Y, Z, zdir='z', zs=0., colors='black',linewidth=6.)

# Computing and plotting principal axes
hessian_matrix = np.array([[ 2.*p[2], p[5]], [p[5], 2.*p[4]]])

eival, eivec = np.linalg.eigh(hessian_matrix)

sorted_indices = np.argsort(eival)[::-1]
eival = eival[sorted_indices]
eivec = eivec[:,sorted_indices]
l1 = eival[0]
l2 = eival[1]

v1 = eivec[:,0]
v2 = eivec[:,1]

xx = np.linspace(-10.,10.,10*N)

y1 = (v1[1]/v1[0])*xx
xx1 = xx[np.abs(y1[:]) <= 10.]
y11 = y1[np.abs(y1[:]) <= 10.]

y2 = (v2[1]/v2[0])*xx
xx2 = xx[np.abs(y2[:]) <= 10.]
y22 = y2[np.abs(y2[:]) <= 10.]

z1 = TargetFunction(p, xx1,y11) - p[0]
z2 = TargetFunction(p, xx2,y22) - p[0]
 
 
principal_axis1 = ax.plot(xx1,y11,z1, zs=0,zdir='z',color='midnightblue',linewidth=2.)
principal_axis2 = ax.plot(xx2,y22,z2, zs=0,zdir='z',color='maroon',linewidth=2.)

#~ y1 = (v1[1]/v1[0])*xx[N/2:3*N/2]
#~ y2 = (v2[1]/v2[0])*xx[N/2:3*N/2]
#~ 
#~ z1 = TargetFunction(p, xx[N/2:3*N/2],y1) - p[0]
#~ z2 = TargetFunction(p, xx[N/2:3*N/2],y2) - p[0]
 #~ 
 #~ 
#~ principal_axis1 = ax.plot(xx[N/2:3*N/2],y1,z1, zs=0,zdir='z',color='black',linewidth=1.)
#~ principal_axis2 = ax.plot(xx[N/2:3*N/2],y2,z2, zs=0,zdir='z',color='black',linewidth=1.)
#print eival




ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])

pl.show()

# Curvature elements
mu = p[2] + p[4]
gamma = 4*p[2]*p[4] - p[5]**2
k1 = mu + np.sqrt(mu**2 - gamma)
k2 = mu - np.sqrt(mu**2 - gamma)

print "k1 = %f" %k1
print "k2 = %f" %k2
print "m  = %f" %mu
print "g  = %f" %gamma
