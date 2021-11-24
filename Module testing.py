# -*- coding: utf-8 -*-
"""
@author: user
"""
from Vectorfunc import Vectorfunc as vec
from raytracer import Ray
import numpy as np
from OpticalElement import OpticalElement
from OpticalElement import SphericalRefraction
#%% Test vector func module

Norm1 = vec.normalisation([1,1,1]) # = 1/ sqrt(1/3)
#vec.check3dlistarray([1,1]) # raise Exception
#vec.check3dlistarray([1,1,1]) # OK
vec.check3dlistarray("hello") # raise TypeError as expected
#%% test ray initiation 

#testray0 = Ray("hello") # Rejects invalid inputs
#testray1 = Ray ([0,0,0],[0,1]) #Rejects incorrect dimensions
#testray2 = Ray ([0,0,0],[0,0,1]) # OK 
testray3 = Ray (np.array([0,0,1]), np.array([0,0,1])) # also OK

# test p() and k()
testray3.p() 
testray3.k()

#test append, all was OK!
testray3.append([0,0,5],[0,1,0])
testray3.p()
testray3.k()
testray3.directions()
testray3.vertices()

#%% Test Optical Element Module 
#Test Intercept function (Task 4)

#Ray Initialtion tested
ray1 = Ray([0,0,0],[0,0,1])
ray2 = Ray([0,0,0],[0,0,-1])
ray3 = Ray([0,0,4],[0,0,1])
ray4 = Ray([0,0,0],[0,1,1]) 
ray5 = Ray ([0,2,0],[0,0,1])

#SphericalRefraction Initial Tested for different cases including plane
Len1 = SphericalRefraction(1,0.2,1,1.5,3)   
Len2 = SphericalRefraction(-1,0.2,1,1,3)
Len3 = SphericalRefraction(1,0,1,1.5,3)  #Flat plane
Len4 = SphericalRefraction(1,0.2,1.5,1,3)
Len5 = SphericalRefraction(1,-0.2,1.5,1,3)

#Len1.intercept(ray1)  #Intersect at [0,0,1] which of course is z0
#Len1.intercept(ray2)  #Does not intersect, return None
#Len1.intercept(ray3)  #Does not intersect because of Aperture size, return None
#Len1.intercept(ray4)  #Correct! Calculated analytical solution: l = (6-sqrt(14))/2 
#Len2.intercept(ray2) #Intersect at [0,0,-1] which of course is z0
#Len3.intercept(ray1) #Intersect at [0,0,1] which of course is z0 for a plane
#Len3.intercept(ray5)  #Intersect at [0,2,1], correct!

#Test Snells function (Task 5)

#Len1.Snells(ray1) # As expected the direction is not altered for a beam in the centre and parrralel to optical axis
#Len1.Snells(ray2) #Return None, since there are no valid intercept
#Len1.Snells(ray4) #Ray seems refracted, checked by numerical calculation

#Now check total internal reflection is detected 

#Len4.intercept(ray4) # Ray does intercept with Len4
#Len4.Snells(ray4) #return None, because of total internal reflection

#Check code also has no bugs for negative curvature 

#Len5.Snells(ray1) #Direction is correct for a negative curvature, implying correctly chosen normal vector
#Len5.Snells(ray2) #return None, No intersect

#Check numerically also for flat plane 
#Len3.Snells(ray1) #No refraction as expected, parralel to optical axis
#Len3.Snells(ray4) # Solution checked numerically

#Since Snells and intercept has been tested for bugs, only need to see
#whether propagate_ray append the points and direction correctly (Task 6)

Len3.propagate_ray(ray4) # Ray is valid
print (ray4.p()) # Both are correct
print (ray4.k())

