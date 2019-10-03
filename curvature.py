#!/usr/bin/env python
# -*- coding: utf-8 -*-
u"""
Python/SciPy prototype of a function to compute the curvature of 
a set of points by fitting a quadratic surface on it.



"""
#
#  Copyright 2015 Florian Blanc <flblanc@cecchini5>
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





import numpy as np
import scipy.optimize

import matplotlib.pyplot as pl



def TargetFunction(p,x,y):
    u"""
    The general form of the quadratic surface. 
    p = array of parameters
    """
    
    return p[0] + p[1]*x + p[2]*(x**2) + p[3]*y + p[4]*(y**2) + p[5]*x*y
    

def ErrorFunction(p,x,y,z):
    u"""
    The error function with respect ot TargetFunction."""
    
    return TargetFunction(p,x,y) - z
    
def JacobianFunction(p,x,y,z):
    u"""
    Jacobian matrix explicitly computed for the quadratic surface.
    /!\ Here x and y are flattened arrays.
    z is here for argument consistency and is not used by the jacobian."""
    
    n = len(x)
    
    J = np.array([ np.ones((n)),x,x**2,y,y**2,x*y ])
    
    return J
    
# ----------------------------
#   TEST ON GENERATED DATA
#-----------------------------


p_real = np.ones((6)) # The TRUE parameters of the quadratic surface

# Point set generation; first I build a grid:

N = 10 # number of generated points

xx = np.linspace(-10.,10.,N)
yy = np.linspace(-10.,10.,N)

X, Y = np.meshgrid(xx,yy)

# Then I compute the actual surface

Z = TargetFunction(p_real, X,Y)


#~ # 2D contour line plotting
#~ fig,ax  = pl.subplots(nrows=1,ncols=1)
#~ ax.contourf(X,Y,Z)
#~ ax.contour(X,Y,Z,20,colors='black')
#~ fig.savefig("test1.png")

# I now add gaussian noise to the surface to model experimental error

sigma=1.
Z += sigma*np.random.randn(N,N)

#~ np.savetxt("X-series.dat",X.flatten(),fmt="%10.5f",header="%i" %N,comments='')
#~ np.savetxt("Y-series.dat",Y.flatten(),fmt="%10.5f",header="%i" %N,comments='')
#~ np.savetxt("Z-series.dat",Z.flatten(),fmt="%10.5f",header="%i" %N,comments='')
np.savetxt("series.dat",np.array([X.flatten(),Y.flatten(),Z.flatten()]).T,fmt="%10.5f",header="%i" %(N*N),comments='')

# 2D contour line plotting
#~ fig,ax  = pl.subplots(nrows=1,ncols=1)
#~ ax.contourf(X,Y,Z)
#~ ax.contour(X,Y,Z,20,colors='black')
#~ ax.set_title(u"""$\sigma = $ %f""" %sigma)
#~ fig.savefig("test2.png")

# Initialization of the problem
p_guess = np.random.random(6)  # randomly drawn initial guess for the parameters
print p_guess
R_init = np.sum(ErrorFunction(p_guess,X,Y,Z)**2)
print "Residual error: %f" %R_init

# Now we can commence the fitting procedure:

#p_fit = scipy.optimize.leastsq(ErrorFunction, p_guess,args=(X.flatten(),Y.flatten(),Z.flatten()))[0]
p_fit = scipy.optimize.leastsq(ErrorFunction, p_guess,args=(X.flatten(),Y.flatten(),Z.flatten()), Dfun=JacobianFunction,col_deriv=1)[0]
print "Fit done."
print p_fit

# Computation of residual error

R_min = np.sum(ErrorFunction(p_fit,X,Y,Z)**2)
print "Residual error after optimization: %f" %R_min


