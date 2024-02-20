import numpy as np

class Grid:
    def __init__(self, N_x: int = 400, N_y: int = 100, num_velocity: int = 9) -> None:
        """
        Attributes
        """
        self.N_x = N_x
        self.N_y = N_y

        self.num_velocity = num_velocity

    def set_cylinder(self, pos_x: float = None, pos_y: float = None, radius: float = 20.0) -> np.ndarray:
        
        # default values for pos_x and pos_y
        pos_x = pos_x if pos_x is not None else self.N_x // 4
        pos_y = pos_y if pos_y is not None else self.N_y // 2

        def distance(x1,y1,x2,y2):
            return np.sqrt((x2-x1)**2 + (y2 - y1)**2)
        
        cylinder = np.full((self.N_y, self.N_x), False)
        for x in range(self.N_x):
            for y in range(self.N_y):    
                # define where cylinder is
                # if less than radius of cylinder --> set boundaries
                if(distance(pos_x, pos_y//2,x,y)) < radius:
                    cylinder[y][x] = True
        return cylinder

class Fluid(Grid):
    def __init__(self, N_x: int = 400, N_y: int = 100, perturbations: bool = True) -> None:
        super().__init__(N_x, N_y)
        self._discrete_vx = np.array([0, 0, 1, 1, 1, 0, -1, -1, -1])
        self._discrete_vy = np.array([0, 1, 1, 0, -1, -1, -1, 0, 1])
        self._weights = np.array([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36])
        self.dist_func = self.initialise(perturbations)
        self.num_velocity = len(self._weights)

    def initialise(self, perturbations: bool = True) -> np.ndarray:
        perturb = 0
        if perturbations:
            perturb = 0.01 * np.random.randn(self.N_y, self.N_x, self.num_velocity)
        return np.ones((self.N_y, self.N_x, self.num_velocity)) + perturb

    ## fluid variables
    def rho(self, F: np.ndarray) -> np.ndarray:
        return np.sum(F, 2)
    
    def vel_x(self, F: np.ndarray) -> np.ndarray:
        return np.sum(F * self._discrete_vx, 2) / self.rho(F)
    
    def vel_y(self, F: np.ndarray) -> np.ndarray:
        return np.sum(F * self._discrete_vy, 2) / self.rho(F)
   
    