import time
import numpy as np
from numba import njit, prange
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import cm


class Simulate:
    def __init__(self, num_iter, psi):
        self.num_iter = num_iter
        self.num_rows = len(psi) 
        self.nums_cols = len(psi[0]) 

    @staticmethod
    @njit(parallel=True)
    def _jacobi(num_iter, psi, rows, cols):
        
        temp_psi = np.zeros((rows + 2, cols + 2))

        for _ in range(num_iter):
            # Compute Jacobi iteration in parallel
            for i in prange(1, rows - 1):
                for j in prange(1, cols - 1):
                    temp_psi[i][j] = 0.25 * (psi[i + 1][j] + psi[i - 1][j] + psi[i][j + 1] + psi[i][j - 1])

            # Update psi in parallel
            for i in prange(1, rows - 1):
                for j in prange(1, cols - 1):
                    psi[i][j] = temp_psi[i][j]
        return psi

    def simulate(self, psi, filename = "jacobi_flow.dat", scaling = 0.1):
        
        print("Simulation start\n")
        tstart = time.time()
        psi = self._jacobi(self.num_iter, psi, self.num_rows, self.nums_cols)
        tend = time.time()
        print("\nCalculation took {0:.5f}s".format(tend - tstart))
        self._write_data(self.num_rows - 2, self.nums_cols - 2, psi, ".\datafiles\\"+ filename, scaling )

    def _write_data(self, m, n, psi, outfile, scaling):
        out = open(outfile, "w")
        ## dimension of box
        out.write("{0} {1}\n".format(m, n))
    
        for i in range(1, m+1):
            for j in range(1, n+1):

                vel_x = (psi[i][j+1] - psi[i][j-1])/2.0
                vel_y = (psi[i-1][j] - psi[i+1][j])/2.0
                speed = (vel_x + vel_y)**2
                    
                out.write("{0:5d} {1:5d} {2} {3} {4}\n".format(i-1, j-1, vel_x, vel_y, speed**scaling))
        out.close()
        print("\nData writing complete")

class Visualise:
    def __init__(self,file_path):
        self.file_path = file_path
        

    def _read_file(self):
        # Open the input file
        input_file = open(self.file_path, "r")

        # Read dimensions   
        line = input_file.readline()
        line = line.rstrip()
        tokens = line.split()
        num_rows = int(tokens[0]) + 2
        num_cols = int(tokens[1]) + 2

        # Define and zero the numpy arrays
        speed = np.zeros((num_rows, num_cols))
        vel_x = np.zeros((num_rows, num_cols))
        vel_y = np.zeros((num_rows, num_cols))

        # Loop over the grid reading the data into the arrays
        for _ in range(1, num_rows-1):
            for _ in range(1, num_cols-1):
                line = input_file.readline()
                line = line.rstrip()
                tokens = line.split()
                
                i1 = int(tokens[0])
                j1 = int(tokens[1])

                vel_x[i1, j1] = float(tokens[2])
                vel_y[i1, j1] = float(tokens[3])
                speed[i1, j1] = float(tokens[4])

        input_file.close()

        return num_rows, num_cols, speed, vel_x, vel_y

    def draw_flow(self, savefig_name = "jacobi_flow.png", colorbar = True):
        matplotlib.rcParams['font.size'] = 8
        matplotlib.rcParams['axes.facecolor'] = 'white'
        matplotlib.use("Agg")

        ## params from file
        num_rows, num_cols, speed, vel_x, vel_y = self._read_file()

        fig, ax = plt.subplots()

        # Regular grids
        x = np.linspace(0, num_rows-1, num_rows)
        y = np.linspace(0, num_cols-1, num_cols)

        # Line widths scaled by the modulus of velocity
        lw = 3 * speed / speed.max()

        # Create the stream lines denoting the velocities
        streamplot = ax.streamplot(x, y, vel_x, vel_y, color='k', density=1, linewidth=lw)

        # Create the heatmap denoting the modulus of the velocities
        heatmap = ax.imshow(speed, interpolation='nearest', cmap=cm.jet)

        if colorbar:
            cbar = fig.colorbar(heatmap, ax=ax, orientation='vertical')
            cbar.set_label('Speed (scaled)')

        # Save the figure to the output PNG file
        fig.savefig(".\images\\"+ savefig_name)
        print("Figure saved")