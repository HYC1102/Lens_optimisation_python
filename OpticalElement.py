# -*- coding: utf-8 -*-
"""

@author: user
"""
import numpy as np
from Vectorfunc import Vectorfunc as vec

class OpticalElement:
    
    def propagate_ray(self, ray):
        "propagate a ray through the optical element"
        raise NotImplementedError()
        
class SphericalRefraction(OpticalElement):
    
    def __init__(self, z0, Curvature, n1, n2, ApertureRadius):
        """Parameters initilisation of a spherical or plane refraction lens 
        object centred at the z axis.
        
        Keyword arguments:
            
        z0 -- the intercept of the surface with the z-axis
        
        Curvature -- the curvature of the surface, positive curvature 
        
        corresponds to a centre at z>z0
        
        n1,n2 -- the refractive indices either side of the surface
        
        ApertureRadius -- the maximum extent of the surface from the optical axis
        
        Note: RadiusCurvature is 1/Curvature. For a plane, the RadiusCurvature 
        is defined as None.
        centre -- the centre of the refractive object
        """
        self._z0 = float(z0)
        self._Curvature = float(Curvature)
        self._n1 = float(n1)
        self._n2 = float(n2)
        self._ApertureRadius = float(ApertureRadius)
        if self._Curvature != 0:
            self._RadiusCurvature = 1 / self._Curvature
            self._centre = np.array([0,0,self._z0 + self._RadiusCurvature])
        else: 
            self._RadiusCurvature = None  #RadiusCurvature is None for a plane
            self._centre = 0              #centre is not important for a plane
            
    def intercept (self, ray):
        """Find intercept of a ray with a Spherical Refraction Element or a 
        Plane in the special case.
        
        Keyword arguments:
        Ray -- Class Ray from the raytracer.py file 
        
        Returns: 
            output: 3-d ndarray - returns the point of intercept of the ray.
                    None - if ray does not intercect with the element. 
        """
        k_hat = vec.normalisation(ray.k())
        if self._RadiusCurvature == None:   #tackles the case when curvature is 0
           l = -(ray.p()[2] - self._z0) / k_hat[2]
           IntersectionPoint = ray.p() + l * k_hat
           #Checks ray is not going backwards and intersects with the plane
           if l > 0 and abs(IntersectionPoint[0]) < self._ApertureRadius and abs(IntersectionPoint[1]) < self._ApertureRadius: 
               return IntersectionPoint
           else: 
               return None
        else:             #Spherical Refraction with +ve/-ve curvature
           r = ray.p() - self._centre
           dotrk_hat = np.dot(k_hat, r)
           discriminant = dotrk_hat**2 - np.dot(r,r) + self._RadiusCurvature**2
           if discriminant > 0:   #ensures that l will be real (i.e intersect)
               l_long = -dotrk_hat + np.sqrt(discriminant)
               l_short = -dotrk_hat - np.sqrt(discriminant)
               
               #selects the correct l based on 4 different cases (wavedirection and curvature) 
#                   l=l_short
#               if k_hat[2]>0 and self._Curvature > 0:    
#               elif k_hat[2]>0 and self._Curvature <0: 
#                   l=l_long
#               elif k_hat[2]<0 and self._Curvature >0: 
#                   l=l_long
#               elif k_hat[2]<0 and self._Curvature <0 and l_short > 0: 
#                   l=l_short              
               #Simplification of the above conditions.
               if k_hat[2] * self._Curvature > 0: 
                   l = l_short
               elif k_hat[2] * self._Curvature < 0:
                   l = l_long
               IntersectionPoint = ray.p() + l*k_hat 
               if l > 0 and abs(IntersectionPoint[0]) < self._ApertureRadius and abs(IntersectionPoint[1]) < self._ApertureRadius: #checks for no backward travelling rays and IntersectionPoint is within the ApertureRadius         
                   return IntersectionPoint
               else:
                   return None
           else: 
               return None
         
    def Snells(self,ray):
        """Find the change in the direction of ray due to refraction
        when it intersects with a spherical / plane element basesd on the 
        vector form of Snell's law.
        
        Keyword arguments:
            
        Ray -- Class Ray from the raytracer.py file 
        
        Returns: 
            output: 3-d ndarray - returns the new direction of the ray after 
                                  refraction.
                    None - if ray does not intercect with the element or total  
                    internal reflection occurs.
        """
        #Checks for no interecpt case: If ray did not interect, intercept 
        #fucntion will return None (not ndarray).
        if not isinstance(self.intercept(ray),np.ndarray):    
            return None
        
        else:
            if self._RadiusCurvature == None:   # Checks for flat plane
                if ray.k()[2] > 0 :                    
                    n_hat = np.array([0,0,-1]) #Define normal vector to be opposite to the ray
                elif ray.k()[2] < 0:
                    n_hat = np.array([0,0,1])   
            else:
                #Normal vector defined to be towards rays direction for all cases 
                n_hat = vec.normalisation(self.intercept(ray) - self._centre) 
                #Considers different scenario of curvature and ray direction
                if self._Curvature * ray.k()[2] > 0: 
                    n_hat = n_hat
                elif self._Curvature * ray.k()[2] < 0:
                    n_hat = -n_hat
            
            k_hat = vec.normalisation(ray.k())
            r = self._n1 / self._n2
            cos1 = -np.dot(n_hat, k_hat)
            determinant = 1 - (r**2) * (1 - cos1**2)   
            if determinant <= 0:      #Checks for total internal reflection
                return None 
            else:
                cos2 = np.sqrt(determinant)
                RefractedRay = r * k_hat + n_hat * (r * cos1 - cos2)
                return RefractedRay
            
    def propagate_ray(self,ray):
        """Propagate the ray through the Spherical Refraction Element. The 
        function appends the intercept of the ray and new direction due to 
        refraction into Ray's point and direction list respectively. 
        
        Keyword arguments:
            
        Ray -- Class Ray from the raytracer.py file 
        
        Returns: 
            output: str - returns the status of the ray and whether ray has 
            been terminated. 
                
        Note: ray._Rayvalid will be set false if there is no valid intercept 
        or total internal reflection occurs. 
        """
        if self.intercept(ray) is None or self.Snells(ray) is None:
            ray._Rayvalid = False
            return f"Ray valid: {ray._Rayvalid}, Ray has been terminated"
        else:
            ray.append(self.intercept(ray),self.Snells(ray))
            return f"Ray valid: {ray._Rayvalid}"
        
class Outputplane(OpticalElement): 

    def __init__(self, z0):
        """Parameters initilisation of a infinite (10e20) detection plane on 
        the z-axis. 
        
        Keyword arguments:
            
        z0 -- the intercept of the surface with the z-axis
        
        Note: 
        Curvature is set to 0. 
        """
        self._z0 = float(z0)
        self._Curvature = 0            #Setting the optical element as a plane.
        
    def Intercept(self, ray):
        """Finds the intercept of the ray with the detection plane. 
        
        Keyword arguments:
            
        ray -- Class Ray from the raytracer.py file 
        
        Returns:
            output: 3-d ndarray - returns the point of intercept of the ray.
                    None - if ray does not intercect with the element. 
        Note: 
            The function is composed of function from SphericalRefraction class,
        Element made local to prohibit access to this object.
        
        ApertureRadius is set to np.inf. 
        """
        #Stored as a local variable to avoid access to other SphericalRefraction functions.
        Element = SphericalRefraction(self._z0, self._Curvature, 1, 1, np.inf) 
        return Element.intercept(ray)
    
    def propagate_ray (self, ray):
        """Propagate the ray through the Outputplane Element. The 
        function appends the intercept of the ray and the current ray direction
        into Ray's point and direction list respectively. 
        
        Keyword arguments:
            
        ray -- Class Ray from the raytracer.py file 
        
        Returns:
            output: str - returns the status of the ray and whether ray has 
            been terminated
        Note: 
        
        ApertureRadius is set to np.inf. 
        
        ray._Rayvalid will be set false if there is no valid intercept 
        or total internal reflection occurs. 
        """
        Element = SphericalRefraction(self._z0, self._Curvature, 1, 1, np.inf) 
        if Element.intercept(ray) is None:
            ray._Rayvalid = False
            return f"Ray valid: {ray._Rayvalid}, Ray has been terminated"
        else:
            ray.append(Element.intercept(ray),ray.k())
            return f"Ray valid: {ray._Rayvalid}"  


