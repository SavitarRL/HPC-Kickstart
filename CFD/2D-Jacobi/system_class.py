import time
import numpy as np
from numba import njit, prange

## only two 
class FluidPorts:
    def __init__(self, inlet_location, port_width, oulet_location):

        self.inlet_location = inlet_location
        self.port_width = port_width 
        self.oulet_location = oulet_location

class Box:
    def __init__(self, length, width, fluidport_instance):
        self.length = length
        self.width = width
        
        self.grid_rows = length 
        self.grid_cols = length 

        self.psi = self.initiate(fluidport_instance) #2D matrix of state function psi

    def _set_inlet(self, psi, inlet_location, oulet_location):
        ## set inlet boundary conditions
        for i in range(inlet_location + 1, inlet_location + oulet_location):
            psi[0][i] = float(i - inlet_location)
        return psi
        
    def _set_outlet(self, psi, inlet_location, oulet_location, port_width):
        
        ## set boundary along horizontal axis
        for i in range(inlet_location + oulet_location, self.grid_cols + 1):
            psi[0][i] = float(oulet_location)
       
        ## set boundary along vertical axis
        for j in range(1, port_width + 1):
            psi[j][self.grid_rows + 1] = float(oulet_location)
        
        ## set diagonal boundaries
        for j in range(port_width + 1, port_width + oulet_location):
            psi[j][self.grid_rows + 1] = float(oulet_location - j + port_width)
        
        return psi

    def initiate(self, fluidport_instance):

        inlet_location = fluidport_instance.inlet_location 
        oulet_location = fluidport_instance.port_width
        port_width = fluidport_instance.oulet_location

        psi = np.zeros((self.grid_rows + 2, self.grid_rows + 2))
        psi = self._set_inlet(psi, inlet_location, oulet_location)
        psi = self._set_outlet(psi, inlet_location, port_width, oulet_location)
        return psi

        