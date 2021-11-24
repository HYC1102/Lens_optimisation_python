
"""
@author: user
""" 
from raytracer import Ray
from OpticalElement import OpticalElement
from OpticalElement import SphericalRefraction
from OpticalElement import Outputplane
from InvestigationTools import Tools as tl
import matplotlib.pyplot as pl
import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
from scipy.optimize import minimize
#%% Test Case 1 Parralle rays (Task 9 & Task 10)
SphericalLens1 = SphericalRefraction(100,0.03,1,1.5,50)
oop = Outputplane(250)   
#oop = Outputplane(200)      #Outputplane at approximate paraxial focal point
Test1 = [SphericalLens1, oop]
ray1 = Ray([0,0,0],[0,0,1])
ray2 = Ray([0.1,0,0],[0,0,1])
ray3 = Ray([0.2,0,0],[0,0,1])
ray4 = Ray([0.3,0,0],[0,0,1])
rays1 = [ray1, ray2, ray3, ray4]
Approximatefocus = tl.focus([SphericalLens1])  # Find focus 
# Focal Length by the lens maker formula
FocalLength = SphericalLens1._RadiusCurvature/(SphericalLens1._n2-1) 
FocalPoint1 = SphericalLens1._centre[2] + FocalLength
tl.xzsimpleplot(Test1,rays1, "Paraxial rays")

print ("Approximate focus is at", Approximatefocus,"Theorectial focus is at", 
       FocalPoint1, ".")


#%% Test Case 2 Rays from the focal point with different directions, refracted 
# to parralel (Task 11)

Test2 = [SphericalLens1, oop]
FocalPoint2 = SphericalLens1._z0 - FocalLength  # Focal Point on the other side of lens
ray5 = Ray([0,0,FocalPoint2],[0.001,0,0.5])
ray6 = Ray([0,0,FocalPoint2],[0.002,0,0.5])
ray7 = Ray([0,0,FocalPoint2],[0.003,0,0.5]) 
ray8 = Ray([0,0,FocalPoint2],[-0.001,0,0.5])
ray9 = Ray([0,0,FocalPoint2],[-0.002,0,0.5])
ray10 = Ray([0,0,FocalPoint2],[-0.003,0,0.5])
ray11 = Ray([0,0,FocalPoint2],[0.,0,0.5]) 
rays2 = [ray5, ray6, ray7, ray8, ray9, ray10, ray11]
tl.xzsimpleplot(Test2, rays2)

#%% Test Case 3 Demonstration of spherical abberation
Test3 = [SphericalLens1, oop]
ray12 = Ray([0,0,0],[0.,0,1])
ray13 = Ray([10,0,0],[0.,0,1])
ray14 = Ray([20,0,0],[0.,0,1])
ray15 = Ray([40,0,0],[0.,0,1])
ray16 = Ray([-10,0,0],[0.,0,1])
ray17 = Ray([-20,0,0],[0.,0,1])
ray18 = Ray([-40,0,0],[0.,0,1])
rays3 = [ray12, ray13, ray14, ray15, ray16, ray17, ray18]
tl.xzsimpleplot(Test3, rays3)

#%% Test Case 4 Ray Bundle and RMS calculation (Task 12 & 13)
oop2 = Outputplane (200)
Test4 = [SphericalLens1, oop2]
rays4 = tl.UniformCollimatedBeam(2.5,0.5,2,[0,0,1], 100)
rays5 = tl.UniformCollimatedBeam(2.5,0.5,2,[0.04,0,1], 100)  

#Diffraction limit
lamb = 650e-6 # in mm
FocalLength = SphericalLens1._RadiusCurvature/(SphericalLens1._n2-1) 
diffraction_limit = 1.22 * lamb * FocalLength / 5
#print ("diffraction limit is =" ,diffraction_limit )

#Spot diagram
#tl.xysimpleplot(Test4, rays4)  #Rays parralel to optical axis 
#tl.xysimpleplot(Test4, rays5)   #Rays with some x-component in rays direction

#Ray trajectories
rays6 = tl.UniformCollimatedBeam(2.5,0.5,1,[0,0,1], 100) # For clear graphs
rays7 = tl.UniformCollimatedBeam(2.5,0.5,1,[0.04,0,1], 100) 
pl.xlim(50,200)
pl.ylim(-4.5,3)
tl.xzsimpleplot(Test4, rays6, colour = "b") 
tl.xzsimpleplot(Test4, rays7, colour = "r")
#%% Modelling a plano - convex lens (Task 15), plotting performance graphs

def performancegraph (elements, end, num):
    """Compute the points of Root Mean Squared radius against the diameter of 
    the Uniform Collimated beam.
        
        Keyword arguments:
            
        elements -- List of Class OpticalElement object from the OpticalElement 
                    file 
        end -- Final value of diameter to compute
        
        num -- Number of data points
        
        Returns: 
            output: diameter: A list of all radius that has been computed
                    rms : A list of rms value corresponded to the diameter of 
                          beam
    """
    rays = []  #List of list of rays
    rms = []
    diameter = np.linspace(0.01, end, num)      #Conver radius to diameter
    z0 = tl.focus(elements)
    elements.append(Outputplane(z0))
    for i in diameter: 
        rays.append((tl.UniformCollimatedBeam(i/2, i/10, 3, [0,0,1], 
                                              elements[0]._z0)))
    
    for ray in rays: 
        rmsvalue, Npoints = tl.xysimpleplot(elements, ray, RMS = True, RMSn = -1)
        rms.append(rmsvalue)
    
    return diameter, rms

#Convex-Plano
SphericalLens2 = SphericalRefraction(100,0.02,1,1.5168,50)
Planosurface1 = SphericalRefraction(105,0,1.5168,1,50)
Convex_Plano = [SphericalLens2, Planosurface1]

#Plano-Convex
Planosurface2 = SphericalRefraction(100,0,1,1.5168,50)
SphericalLens3 = SphericalRefraction(105,-0.02,1.5168,1,50)
Plano_Convex = [Planosurface2, SphericalLens3]

#Performance Graph plots
diameter1, rms1 = performancegraph(Convex_Plano, 10, 50)
diameter2, rms2 = performancegraph(Plano_Convex, 10, 50)  
#%% Performance graph plotting and fitting 
def exponential(x, A, c, h):
    return A* np.exp(x + c) + h
popt, err = sp.optimize.curve_fit(exponential, diameter1, rms1, p0=[0.5, -1, 0])
#pl.plot(diameter1, rms1, '.')
#pl.plot (diameter1, exponential(diameter1, *popt)) # Exponential fit it poor

#Plot and fit quadratic
par1, residuals1, _, _, _= np.polyfit(diameter1, rms1, 2, full = True)
par1, cov1= np.polyfit(diameter1, rms1, 2, cov = True)
par2, residuals2, _, _, _= np.polyfit(diameter2, rms2, 2, full = True)
par2, cov2= np.polyfit(diameter2, rms2, 2, cov = True)

chisq1 = residuals1 # degrees of freedom = N - 4 (3 parameters to fit) 
chisq2 =  residuals2
pl.xlabel("Initial beam diameter / mm")
pl.ylabel("Root mean square radius / mm")
pl.plot(diameter1, rms1, '.', label = "Convex-Plano")
pl.plot(diameter1, diameter1**2 * par1[0] + diameter1 * par1[1] + par1[2])
pl.plot(diameter2, rms2, 'r.', label = "Plano-Convex")
pl.plot(diameter2, diameter2**2 * par2[0] + diameter2 * par2[1] + par2[2])
pl.legend()

#Root mean square comparision at 10mm
#print (rms2[-1] / rms1[-1], "times smaller rms for the convex-plano at 10mm diameter")

#Fitted parameters to the quadratic
#print ("Convex-Plano: Chi Square value is", chisq1[0], "Parameter a =", 
#       round(par1[0],6), "+/-", round(np.sqrt(cov1[0][0]), 6), "Parameter b =", 
#       round(par1[1],6), "+/-", round(np.sqrt(cov1[1][1]), 6), "Parameter c =",
#       round(par1[2],6), "+/-", round(np.sqrt(cov1[2][2]), 6))

#print ("Plano-Convex: Chi Square value is", chisq2[0], "Parameter a =", 
#       round(par2[0],6), "+/-", round(np.sqrt(cov2[0][0]), 6), "Parameter b =", 
#       round(par2[1],6), "+/-", round(np.sqrt(cov2[1][1]), 6), "Parameter c =",
#       round(par2[2],6), "+/-", round(np.sqrt(cov2[2][2]), 6))

    #%% Lens Optimisation
focal_point = 200  #Set focal point of the lens (location of outputplane)

def rmsfunc (curvature): 
    """Compute the rms of a singlet lens as a variable of 2 curvatures for the
        2 surfaces to be separated 5mm apart. The ray used is a 10mm diameter
        uniform circular beam.
        
        Keyword arguments:
            
        curvature -- Size 2 list with the 2 values of curvature
        
        Returns: 
            output: rms: float 
                    rms: None if rays are invalid
        Note: Focal point/focal length can be set by changing the focal_point
        variable defined globally.  
    """
    curvature_1 = curvature[0]
    curvature_2 = curvature[1]
    count = 0
    Sumsquared = 0
    
    SphericalLens1 = SphericalRefraction(100,curvature_1,1,1.5168,50)
    SphericalLens2 = SphericalRefraction(105,curvature_2,1.5168,1,50)
    lens = [SphericalLens1, SphericalLens2]
    rays = tl.UniformCollimatedBeam(5,1,2,[0,0,1])
    # Set the focal point of the optimised lens
    lens.append(Outputplane(focal_point))    
    
    for elem in lens:
        for ray in rays:
            elem.propagate_ray(ray)
            
    for ray in rays: 
        if ray._Rayvalid:
            count += 1
            vertices = np.array(ray.vertices())
            Sumsquared += (vertices[-1][0]) ** 2 + (vertices[-1][1]) ** 2 
    if count != 0:
        rms = np.sqrt(Sumsquared / count) 
        
    else:
        rms = None
    return rms

initial_guess = [0, 0]
#initial_guess = [ 0.01511615, -0.00480703] 

#Perform Optimisation of lens
sp.optimize.minimize(rmsfunc, initial_guess,method= "L-BFGS-B", bounds = 
                     ((0,np.inf),(-np.inf,0))) #return success false, change initial guess
#%% Optimised values, lowest rms is chosen

print (rmsfunc([ 0.01511615, -0.00480703])) #using L-BFGS-B
print (rmsfunc([ 0.01172145, -0.00834972])) #using TNC
print (rmsfunc([ 0.01418388, -0.00578591])) #using SLSQP

rays6 = tl.UniformCollimatedBeam(2.5,0.5,1,[0,0,1], 100)
OptimisedElement = [SphericalRefraction(100, 0.01511615, 1, 1.5168, 50), 
          SphericalRefraction(105, -0.00480703, 1.5168, 1, 50 )]
focus_optimised = tl.focus(OptimisedElement)
OptimisedElement.append(Outputplane(focus_optimised))
tl.xzsimpleplot(OptimisedElement, rays6)


