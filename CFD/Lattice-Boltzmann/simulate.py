from typing import Union
from components import Grid, Fluid
import numpy as np
import matplotlib.pyplot as plt
import os

class Simulate:
    def __init__(self, fluid_instance: Fluid, num_timesteps: int, timescale: float, wall_boundary: bool = True) -> None:
        self.num_timesteps = num_timesteps
        self.timescale = timescale
        self.wall_boundary = wall_boundary

        self.fluid_instance = fluid_instance
        self.cylinder = self.fluid_instance.set_cylinder()

        self.details = self._simulation_details()        

    def _simulation_details(self) -> dict:
        details_dict = {
            "size": "{}x{}".format(self.fluid_instance.N_x, self.fluid_instance.N_y),
            "timesteps": "{}".format(self.num_timesteps),
            "timescale": "{}".format(self.timescale),
            "wall boundary": "{}".format(self.wall_boundary)
        }
        return details_dict


    def _file_detail_suffix(self)-> str:
        suffix = ""
        for _, values in self.details.items():
            suffix = suffix + "_" + values 
        return suffix

    def apply_wall_boundary(self, F: np.ndarray) -> np.ndarray:
        F[:, -1, [6,7,8]] = F[:, -2, [6,7,8]]
        F[:, 0, [2,3,4]] = F[:, 1, [2,3,4]]
        return F
    
    def apply_obstacle_boundary(self, F: np.ndarray) -> Union[np.ndarray, np.ndarray, np.ndarray]:

        """
        calculate collisions with cylinder:
        getting all points where velocities in cylinder
        at each point --> invert velocities
        """
        
        boundary_F = F[self.cylinder, :]
        boundary_F = boundary_F[:,[0,5,6,7,8,1,2,3,4]]

        # rho = Fluid.rho(F = F)
        vel_x = self.fluid_instance.vel_x(F = F)
        vel_y = self.fluid_instance.vel_y(F = F)

        F[self.cylinder, :] = boundary_F
        ## no fluid movement in cylinder
        vel_x[self.cylinder] = 0
        vel_y[self.cylinder] = 0

        return (F, vel_x, vel_y)

    
    def streaming(self, F:np.ndarray)->np.ndarray:
        """
        streaming step:
        take every single nodal velocity and move to its corresponding neighbour
        in the direction of its discrete velocity
        """
        for i, v_x, v_y in zip(range(self.fluid_instance.num_velocity), 
                               self.fluid_instance._discrete_vx, 
                               self.fluid_instance._discrete_vy):
            F[:,:,i] = np.roll(F[:, :,i], v_x,axis = 1)
            F[:,:,i] = np.roll(F[:, :,i], v_y,axis = 0)
        return F
    
    # BGK approximation
    def apply_BGK_update(self, F: np.ndarray, rho: np.ndarray, vel_x: np.ndarray, vel_y: np.ndarray)->np.ndarray:

        F_eq = np.zeros(F.shape)

        ## streaming (L.H.S)
        for i, v_x, v_y, w in zip(range(self.fluid_instance.num_velocity), 
                                  self.fluid_instance._discrete_vx, 
                                  self.fluid_instance._discrete_vy, 
                                  self.fluid_instance._weights):
            
            F_eq[:,:,i] =  rho * w * ( 1 + 3*(v_x*vel_x + v_y*vel_y)  
                                    + 9*(v_x*vel_x + v_y*vel_y)**2/2 
                                    - 3*(vel_x**2+vel_y**2)/2 )
            
        ## collision (R.H.S)
        F += -(1.0/self.timescale) * (F - F_eq)

        return F
    
    
    def simulate(self, filename:str, mode: str = "speed", pause_step: float = 0.001)->None:
        """
        Simulate fluid dynamics and visualize the results.

        Parameters:
            filename (str): Filename to save the final simulation image.
            mode (str): Visualization mode. Options: "speed" (default), "vortices".
            pause_step (float): Pause duration between steps in the visualization.

        Raises:
            ValueError: If an invalid mode is provided.
        """
        if mode != "speed" and mode != "vortices":
            raise ValueError("Invalid mode {}. Choose either 'speed' or 'vortices'.".format(mode))
        
        dist_func = self.fluid_instance.initialise()
        
        # initially flowing to right --> assign 3rd value of every node to be non-zero
        ## can be different
        dist_func[:,:,3] = 2.3

        for dt in range(self.num_timesteps + 1):
            
            # applying absorbing wall conditions 
            if self.wall_boundary:
                dist_func = self.apply_wall_boundary(dist_func)
            
            # streaming step
            dist_func = self.streaming(dist_func)

            # collision with cylinder, returning distribution function and velocities
            dist_func, vel_x, vel_y = self.apply_obstacle_boundary(dist_func)

            # apply BGK approx
            dist_func = self.apply_BGK_update(dist_func, 
                                              rho = self.fluid_instance.rho(F=dist_func), 
                                              vel_x = vel_x, 
                                              vel_y = vel_y)

            if mode == "speed" and dt % 5 == 0:
                im = plt.imshow(np.sqrt(vel_x**2 + vel_y**2))
                cbar = plt.colorbar(im, shrink=0.5)  # Add color bar and set its size
                cbar.set_label('Speeds')
                plt.title("Time {}".format(dt))
                plt.pause(pause_step)
                
                if dt < self.num_timesteps:
                    cbar.remove()
                    plt.cla()
                    
            
            elif mode == "vortices" and dt % 5 == 0:
                # curl, as by mathematical definition
                dfy_dx = vel_x[2:,1:-1] - vel_x[0:-2, 1:-1]
                dfx_dy = vel_y[1:-1,2:] - vel_y[1:-1, 0:-2]

                curl_F = dfy_dx - dfx_dy
                im = plt.imshow(curl_F, cmap="bwr")
                cbar = plt.colorbar(im, shrink=0.5)  # Add color bar and set its size
                cbar.set_label('Curl')
                plt.title("Time {}".format(dt))
                plt.pause(pause_step)
                
                if dt < self.num_timesteps:
                    cbar.remove()
                    plt.cla() 

                    
                
        file_details = filename + "_" + self._file_detail_suffix() + "_" + mode 

        # Constructing the image path
        img_name = file_details + ".png"
        img_path = os.path.join(os.getcwd(),"Lattice-Boltzmann", "images", img_name)
        plt.savefig(img_path ,dpi=240)
        
            
        