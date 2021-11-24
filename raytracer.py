"""


@author: user
"""
import numpy as np
from Vectorfunc import Vectorfunc as vec
class Ray:
    
    def __init__(self, initial_point= [0,0,0], initial_direction = [0,0,0]):
        """Parameters initilisation of a ray.
        
        Keyword arguments:
            
        Intial_point -- size 3 ndarray
        
        initial_direction -- size 3 ndarray
        
        Note: Rayvalid is defined as True initially, when ray is terminated, 
        it will set false. 
        """
        self._p = [np.array(vec.check3dlistarray(initial_point))]
        self._k = [np.array(vec.check3dlistarray(initial_direction))]
        self._Rayvalid = True
        
    def p(self):
        """Return the current point of the ray.
        
        Returns:
            output - size 3 ndarray, current point of the ray
        """
        return self._p[-1]
    
    def k(self):
        """Return the current direction of the ray.
        
        Returns:
            output - size 3 ndarray, current direction of the ray
        """
        return self._k[-1]
    
    def append(self, new_point, new_direction):
        """Append a point and a direction to the ray. 
        
        Keyword arguments:
            
        new_point -- size 3 ndarray: Point to append into Ray.
        
        new_direction -- size 3 ndarray: Direction to append into Ray.
        """
        self._p.append(np.array(vec.check3dlistarray(new_point)))
        self._k.append(np.array(vec.check3dlistarray(new_direction)))

    def vertices (self):
        """Return the vertices of the ray.
        
        Returns:
            output: list - vertices of the ray.
        """
        return self._p
    
    def directions (self):
        """Return the directions of the ray.
        
        Returns:
            output - list: directions of the ray.
        """
        return self._k

    def rayvalid(self):
        """Return the status of the ray
        
        Returns:
            output: bool - True if ray is valid. If ray terminated, returns 
                           False.
        """
        return self._rayvalid 
    
    def __str__(self):
        return f"{self._p,self.k()}"
    

