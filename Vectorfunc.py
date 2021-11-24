# -*- coding: utf-8 -*-
"""

@author: user
"""
import numpy as np

class Vectorfunc:
    
    def normalisation (vector):
        """Normalise a vector to a normal vector.
        
        Keyword arguments:
            
            vector -- Array-like object
        
        Returns:
            
            output : ndarray - returns the normalised vector
                     None - if magniture of vector is 0
        """   
        if np.linalg.norm(vector) != 0:        #ensures the magnitude is not is not 0
            norm = vector / np.linalg.norm(vector)
            return norm
        else: 
            raise ZeroDivisionError ("Magnitude of list or array is 0, ray is not propagating" )
 
    def check3dlistarray (obj):       
        """Raise TypeErrors or Exceptions when input object is not a 3d list or 
        ndarray.
        
        Keyword arguments:
            
            vector -- Array-like object
        
        Returns:
            
            output : obj - returns the 3d list or ndarray
    
        Raise:
            TypeError: When object is not a list or ndarray
            Exception: When the list or array is not size 3
        """   
        if not isinstance (obj,(list,np.ndarray)):
            raise TypeError ("Input parameter must be a size 3 ndarray or list")
        
        elif len(obj) != 3:
            raise Exception ("Input parameter must be a size 3 ndarray or list")
        
        else:
            return obj