from simulate import Simulate
from components import Fluid

"""
With credits to Philip Mocz (2020) Princeton Univeristy, as from tutorial:
https://www.youtube.com/watch?v=JFWqCQHg-Hs

Simulating flow past cylinder for an isothermal fluid
"""


def main():
	
    # Create an instance of the Fluid class
	fluid_instance = Fluid(N_x=400, N_y=100,perturbations=True)
    
    
    # Create a simulation instance
	num_timesteps = 1000
	timescale = 0.6
	wall_boundary = True
	simulation = Simulate(fluid_instance, num_timesteps, timescale, wall_boundary)

	filename = "LatticeBoltzmann"

    # Run the simulation with default settings (mode='speed', pause_step=0.001)
	simulation.simulate(filename, mode="speed", pause_step=0.001)

    # # Run the simulation with 'vortices' mode and a longer pause step
	simulation.simulate(filename, mode='vortices', pause_step=0.001)

if __name__ == "__main__":
    main()

