# -*- coding: utf-8 -*-
"""


@author: user
"""
import matplotlib.pyplot as pl
import numpy as np
from raytracer import Ray 


class Tools: 
    
    def xzsimpleplot (elements, rays, label = None, colour = "b"):
        """Propagate a bundle of rays through the elements in the order given
        by the element list and produce a x-z plot for tracing ray trajectories.
        
        Keyword arguments:
            
        elements -- List of Class OpticalElement object from the OpticalElement 
                    file 
        
        rays -- List of Class Ray object from the raytracer.py file 
        
        label -- string: label of the line and produce legend on plot

        Returns: 
            output: str - returns number of points plotted on the graph
            
                    plot: x-z axis plot showing ray trajectories
        """
        count = 0
        for elem in elements:
            for ray in rays:
                elem.propagate_ray(ray)
    
        for ray in rays: 
            if ray._Rayvalid:
                count += 1
                vertices = np.array(ray.vertices())
                pl.xlabel("z / mm")
                pl.ylabel("x / mm")
                pl.plot(vertices [:,2], vertices [:,0], colour)

                
            
        return f"Number of valid rays plotted: {count}"
            
    def xysimpleplot (elements, rays, RMS = True, RMSn = -1):
        """Propagate a bundle of rays through the elements in the order given
        by the element list and produce 2 x-y plot (Input and Output
        spot diagrams).
        
        Keyword arguments:
            
        elements -- List of Class OpticalElement object from the OpticalElement 
                    file 
        
        rays -- List of Class Ray object from the raytracer.py file 
        
        RMS -- bool : If True returns the value of the root mean squared 
        
        devations from centre spot. Default is True. 
        
        RMSn -- int: Specify the Plane to calculate Root Mean Square.
        Default is -1 (Outputplane).

        Returns: 
            output: str - returns number of points plotted on the graph
            
                    plot: x-z axis plot showing ray trajectories
        """
        count = 0
        Sumsquared = 0
        for elem in elements:
            for ray in rays:
                elem.propagate_ray(ray)
        #Used to calculate rms when ray does not focus at x = 0 and y = 0
        centreray_xpos = rays[0].vertices()[RMSn][0] 
        centreray_ypos = rays[0].vertices()[RMSn][1] 
        #Graph plotting
        fig, (ax1, ax2) = pl.subplots(1, 2)
        pl.subplots_adjust(wspace = 0.25)
        ax1.set_xlabel('x / mm')
        ax1.set_ylabel('y / mm')
        ax1.set_title ("z = 0")
        ax2.set_title ("Focus")
        ax2.set_xlabel('x / mm')
        pl.tight_layout()
        
        for ray in rays: 
            if ray._Rayvalid:
                vertices = np.array(ray.vertices())
                
                ax1.set_aspect(1./ax1.get_data_ratio())
                ax1.plot(vertices [0][0], vertices [0][1], "b.")
        
                ax2.set_aspect(1./ax2.get_data_ratio())
                ax2.plot(vertices [-1][0],vertices [-1][1],"b.")
                Sumsquared += ((vertices[RMSn][0] - centreray_xpos) ** 2 + 
                               (vertices[RMSn][1] - centreray_ypos) ** 2 )
                count += 1
        if RMS and count != 0:
            rms = np.sqrt(Sumsquared / count)                
        else:
            rms = None
        return rms, f"Number of valid points plotted: {count}"

    def UniformCollimatedBeam (radius, Radial_distance, Points_per_unit_circumference, 
                               ray_direction = [0,0,1], alignzpos = 0):
        """Generate a uniform circular bundle of collimated rays. The function 
        generate bundle of rays through circular layers.
        
        Keyword arguments:
        radius -- Radius of the circular bundle of rays.
        
        Radial_distance -- Separation between each circular layer of points. 
        
        Points_per_unit_circumference -- Specify the number of points per 1 
        unit of circumference in each layer.
        
        ray_direction -- The direction of the colliated ray. Default is [0,0,1]
        
        alignzpos -- Specify the z position to centre the bundle at. Default is 0.

        Returns: 
            output: str - returns number of points plotted on the graph
            
                    plot: x-z axis plot showing ray trajectories
        """
        if not isinstance (ray_direction,(list,np.ndarray)):
            raise TypeError ("Input parameter must be a size 3 ndarray or list")
        elif len(ray_direction) != 3:
            raise Exception ("Input parameter must be a size 3 ndarray or list")
        else: 
            Rays = []
            radial = np.arange(0, radius + Radial_distance, Radial_distance)
            k = alignzpos / ray_direction [2]   
#Calculate transaltion requrired to bring ray to x = 0 and y = 0            
            xtrans = -ray_direction[0] * k       
            ytrans = -ray_direction[1] * k
            for i in range (len(radial)): 
                if radial[i] == 0:      #Produces the centre spot
                    Rays.append(Ray([0 + xtrans,0 + ytrans,0],ray_direction))
                else:
                    angle = np.linspace(0, 2*np.pi, int(2 * np.pi * radial[i] 
                    * Points_per_unit_circumference))
                    for j in range (len (angle)):    
                        Rays.append(Ray([radial[i]*np.cos(angle[j]) + xtrans, 
                        radial[i] * np.sin(angle[j]) + ytrans, 0], ray_direction))
        return Rays

    
    def focus (elements):
        """Finds the approximate focus of a lens system assuming that it is a 
        converging lens. The system must not contain an outputplane. 
        
        Keyword arguments:
            
        elements -- List of Class OpticalElement object from the OpticalElement 
                    file 
        
        Returns: 
            output: float: the approximate paraxial focus.  
        """
        #Ray close to centre to find the paraxial focus.
        ray = Ray([0.1,0,0],[0,0,1])  
        for elem in elements:
              elem.propagate_ray(ray)
        t = -ray.p() [0] / ray.k()[0]
        focus = ray.p()[2] + t * ray.k()[2]
        return focus
        
  
              