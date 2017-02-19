'''This is a demonstration component model (Hypomodel) which can 
be incorporated with other models in a composite modelling (Hypermodelling) 
framework. It is being developed as part of the CHIC project.

This is the Tumour component. It takes in spatial oxygen concentration
and returns spatial density of cells in each stage of the cell cycle.
'''

from scipy.integrate import odeint
from scipy import ndimage
import sys, os
import numpy as np
import time
from watchdog.observers import Observer
import warnings

import ChicCommon

################## Functions ##################
def cellCycleSystem(y, t, k):
    
    '''Return the operator for the cell cycle ODE system.
    
    The system is specified in the form:
        dy/dt = f(y, t) return f.
        
    Here k is a vector of cell cycle constants and y is a
    vector of fractions of cells in each stage of the cycle.
    '''
    
    deriv = np.array([
                      2.0 * k[3] * y[3] + k[5] * y[4] - k[0] * y[0],
                      k[0] * y[0] - k[1] * y[1],
                      k[1] * y[1] - k[2] * y[2],
                      k[2] * y[2] - k[3] * y[3] - k[4] * y[3],
                      k[4] * y[3] - k[5] * y[4]
                      ])
    return deriv

################## Classes ##################
class CellSimulation(ChicCommon.Simulation):
    
    '''Class for managing simulations
    '''   
    
    def add_grid(self):
        self.grid = CellGrid(self.extents, self.spacing)
        
    def apply_cell_constants(self, parameters):
        cellConstants = {"k1s": float(parameters["k1s"]),
                         "ks2": float(parameters["ks2"]),
                         "k2m": float(parameters["k2m"]),
                         "km1": float(parameters["km1"]),
                         "km0Base": float(parameters["km0Base"]),
                         "k01Base": float(parameters["k01Base"]),
                         "thresholdOx": float(parameters["thresholdOx"]),
                         "maxCellNumber": float(parameters["maxCellNumber"])}

        for eachCell in self.grid.cells:
            eachCell.set_constants(cellConstants)  
            
    def get_input(self, fileName, inputList):
        inputList = self.read_input(fileName, inputList)     
    
        for idx, eachCell in enumerate(self.grid.cells):
            if inputList["Number of Cells"]:
                eachCell.self.cellNumber = inputList["Number of Cells"][idx]
            if inputList["G1 Cells"]:
                eachCell.cellPopulation[0] = inputList["G1 Cells"][idx]
            if inputList["S Cells"]:
                eachCell.cellPopulation[1] = inputList["S Cells"][idx]
            if inputList["G2 Cells"]:
                eachCell.cellPopulation[2] = inputList["G2 Cells"][idx]
            if inputList["M Cells"]:
                eachCell.cellPopulation[3] = inputList["M Cells"][idx]
            if inputList["G0 Cells"]:
                eachCell.cellPopulation[4] = inputList["G0 Cells"][idx]
            if inputList["Oxygen Concentration"]:
                eachCell.concOx = inputList["Oxygen Concentration"][idx]
            if inputList["Interface Distance"]:
                eachCell.distance = inputList["Interface Distance"][idx] 
                  
    def run(self):
        tol = 0.01
        self.grid.update_cell_distances(self.grid.cells)
        self.write_output(fileName = self.name + "0")
        
        for idx in range(self.numSteps):
            
            # Cycle the cells
            for eachCell in self.grid.cells:
                eachCell.update(self.timeIncrement)
               
            # Redistribute cells
            for eachEntry in self.grid.distanceList:
                thisCell = self.grid.cells[eachEntry]
                if thisCell.cellNumber > thisCell.maxCellNumber + tol:
                    thisCell.redistribute_cells()  
            self.grid.update_cell_distances(self.grid.cells)
            
            # Write output
            if idx%self.outputFrequency == 0:
                self.write_output(fileName = self.name + str(idx + 1))  
                
    def step(self, parentInput = None, fileName = None):    
        tol = 0.01
        
        if parentInput:
            result = [] 
        
        if fileName:
            self.read_input(fileName) 
            
        self.grid.update_cell_distances(self.grid.cells) 
        
        # Cycle the cells
        for eachCell in self.grid.cells:
            eachCell.update(self.timeIncrement)

        # Redistribute cells
        for eachEntry in self.grid.distanceList:
            thisCell = self.grid.cells[eachEntry]
            if thisCell.cellNumber > thisCell.maxCellNumber + tol:
                thisCell.redistribute_cells()        
            
        if parentInput:
            for eachCell in self.grid.cells:
                result.append([eachCell.location, 
                               eachCell.cellNumber])
        # Write output
        if self.increment%self.outputFrequency == 0:
            self.write_output(fileName = self.name + str(self.increment + 1))       
        
        self.increment = self.increment + 1
        
        if parentInput:
            return result
        
class CellGrid(ChicCommon.Grid):
    
    '''Class for domain level storage and operations.
    
    This class stores the definition of the computational
    grid and a list of geometric cell objects at each grid
    point. Methods involving inter-cell operations are
    included in this class.
    '''
    
    def __init__(self, extents = (10, 10, 10), spacing = 1.0, 
                 neighbourhood = "Neumann"):
        
        ''' On initialization a grid is created and filled with
        cells.
        '''
        
        self.extents = extents
        self.extentsR = extents[::-1]
        self.spacing = spacing
        self.distanceMap = []
        self.distanceList = []
        self.outerCells = None
        self.cells = []
        
        nx, ny, nz = self.extents
        self.numCells = nx * ny * nz
        idx, idy, idz = np.unravel_index(np.arange(self.numCells), 
                                 dims = self.extents, order = "F")
       
        for count in range(len(idx)):
            point = (idx[count] * self.spacing, 
                     idy[count] * self.spacing,
                     idz[count] * self.spacing)
            cell = CellCell(location = point, index = count)
            self.cells.append(cell)
        
        self.update_cell_neighbours(self.cells, neighbourhood)
        
    def update_cell_neighbours(self, cells, neighbourhood = "Neumann"):
        
        self.update_neighbour_list(neighbourhood)
        for eachCell in cells:
            eachCell.neighbours = []
            for eachNeighbour in self.neighbours[eachCell.index]:
                if eachNeighbour != eachCell.index:
                    eachCell.neighbours.append(self.cells[eachNeighbour])   
            
    def update_neighbour_list(self, neighbourhood = "Neumann"):
        
        '''Find the neighbours of each cell using the Neumann or Moore
        neighbourhoods.
        '''
   
        numCells = len(self.cells)       
                  
        idx, idy, idz = np.unravel_index(np.arange(numCells), 
                                         dims = self.extents, order = "F")
        
        # Use a von Neumann or Moore Neighbourhood
        if neighbourhood == "Neumann":
            neigh_idx = np.vstack((idx-1, idx+1, idx, idx, idx, idx))
            neigh_idy = np.vstack((idy, idy, idy-1, idy+1, idy, idy))
            neigh_idz = np.vstack((idz, idz, idz, idz, idz -1, idz +1))
        else:
            neigh_idx = np.vstack((idx-1, idx, idx+1, idx-1, idx, idx+1, 
                                   idx-1, idx, idx+1, idx-1, idx, idx+1,
                                   idx-1, idx+1, idx-1, idx, idx+1, 
                                   idx-1, idx, idx+1, idx-1, idx, idx+1,
                                   idx-1, idx, idx+1))
            neigh_idy = np.vstack((idy-1, idy-1, idy-1, idy, idy, idy,   
                                   idy+1, idy+1, idy-1, idy-1, idy-1, idy-1,
                                   idy, idy, idy+1, idy+1, idy+1, 
                                   idy-1, idy-1, idy-1, idy, idy, idy,
                                   idy+1, idy+1, idy+1))
            neigh_idz = np.vstack((idz-1, idz-1, idz-1, idz-1, idz-1, idz-1, 
                                   idz-1, idz-1, idz-1, idz, idz, idz,  
                                   idz, idz, idz, idz, idz, 
                                   idz+1, idz+1, idz+1, idz+1, idz+1, idz+1,
                                   idz+1, idz+1, idz+1))           
        
        # Get neighbour indices in list form                
        self.neighbours = np.ravel_multi_index((neigh_idx, neigh_idy, neigh_idz), 
                                    dims = self.extents, mode='clip', order = "F").T
                              
    def update_cell_distances(self, cells):
        self.update_distance_list()
        for eachCell in cells:
            eachCell.distance = self.distanceMap[eachCell.index]
            
    def update_distance_list(self):
        self.update_distance_map()
        x = np.array(self.distanceMap)
        reverse = np.argsort(x)
        self.distanceList = reverse[::-1]
        
    def update_distance_map(self):
        self.outerCells = []
        for eachCell in self.cells:
            if eachCell.cellNumber == 0:
                self.outerCells.append(eachCell.index)
        if len(self.outerCells) == 0:
            raise SystemExit("No free space found for cell propagation.")
            
        interfaceMap = np.ones(shape = self.extentsR)
        interfaceMap[np.unravel_index(self.outerCells, dims = self.extentsR)] = 0.0
            
        # Get distance transform
        distanceMap = ndimage.distance_transform_edt(interfaceMap, 
                                                          sampling = self.spacing)
        
        # Update distance list
        indices = np.arange(self.numCells)
        self.distanceMap = distanceMap[np.unravel_index(indices, dims = self.extentsR)]
        
class CellCell(ChicCommon.GeometricCell):
        
    def further_initialization(self):
        self.neighbours = []
        self.distance = -1.0
        self.time = 0.0
        self.updateTime = 0.0
        self.set_cell_population()
        self.concOx = 0.0
                  
    def set_cell_population(self, cellPopulation = None):
        if cellPopulation is None:
            self.cellPopulation = np.array([0.0, 0.0, 0.0, 0.0, 0.0]) 
        else:
            self.cellPopulation = cellPopulation
        self.cellNumber = np.sum(self.cellPopulation)
        
    def set_constants(self, constants):
        self.k1s = constants["k1s"]
        self.ks2 = constants["ks2"]
        self.k2m = constants["k2m"]
        self.km1 = constants["km1"]
        self.km0_base = constants["km0Base"]
        self.k01_base = constants["k01Base"]
        self.thresholdOx = constants["thresholdOx"]
        self.maxCellNumber = constants["maxCellNumber"]
        
    def update(self, dt):
        if self.cellNumber > 0:
            k = [0] * 6
            k[0] = self.k1s
            k[1] = self.ks2
            k[2] = self.k2m
            k[3] = self.km1
            
            if self.concOx >= self.thresholdOx:
                k[5] = self.k01_base
            else:
                k[4] = self.km0_base
    
            timePoints = [self.time, self.time + dt]
            solution = odeint(cellCycleSystem, self.cellPopulation, timePoints, args = (k,))
            self.cellPopulation = solution[1]
            self.cellNumber = np.sum(self.cellPopulation)

        self.time = self.time + dt
        
    def redistribute_cells(self):
        
        ##t todo: account for cell removal using type 'remove'. Only call
        # on non-interface cells and if a neighbouring interface cell is
        # empty the present cell can begin to empty.
        
        # Get fraction of excess cells
        frac = self.cellNumber/self.maxCellNumber
        
        # Get number of each population
        transfers = [0] * len(self.cellPopulation)
        for idx, eachPopulation in enumerate(self.cellPopulation):
            transfers[idx] = eachPopulation * (frac - 1.0) / frac
            self.cellPopulation[idx] = eachPopulation - transfers[idx]
            
        # Identify neighbours that can accept cells
        updatableNeighbours = []
        for eachNeighbour in self.neighbours:       
            if eachNeighbour.distance < self.distance:
                updatableNeighbours.append(eachNeighbour)
                
        numNeighbours = len(updatableNeighbours) 
    
        if numNeighbours > 0: 
            for eachNeighbour in updatableNeighbours:
                for idx, eachPopulation in enumerate(eachNeighbour.cellPopulation):
                    eachNeighbour.cellPopulation[idx] = eachPopulation + \
                        transfers[idx] / numNeighbours 
                eachNeighbour.cellNumber = np.sum(eachNeighbour.cellPopulation)

            self.cellNumber = np.sum(self.cellPopulation)
            self.updateTime = self.time
            
        else:   
            message = "Cell" + str(self.index) + " failed to redistribute its contents at time: " +str(self.time) + ". The tumour has reached the domain boundary." 
            warnings.warn(message)  
            sys.exit(1)  
            
    
    def get_output(self):
        
        output = {"Number of Cells": self.cellNumber,
                "G1 Cells": self.cellPopulation[0],
                "S Cells": self.cellPopulation[1],
                "G2 Cells": self.cellPopulation[2],
                "M Cells": self.cellPopulation[3],
                "G0 Cells": self.cellPopulation[4],
                "P Cells": self.cellNumber - self.cellPopulation[4],
                "Interface Distance": self.distance,
                "Oxygen Concentration": self.concOx}
        
        return output
    
if __name__ == "__main__":
    
    # Parse the simulation inputs
    inputFile = ChicCommon.parse_command_line_args(sys.argv[1:], "CellComponent")
    
    parameters = {"name": "CellSimulation",
                  "numSteps": 1,
                  "timeIncrement": 1.0,
                  "outputFrequency": 1,
                  "spacing": 1.0,
                  "numX": 1,
                  "numY": 1,
                  "numZ": 1,
                  "k1s": 0.1,
                  "ks2": 0.1,
                  "k2m": 0.1,
                  "km1": 0.1,
                  "km0Base": 0.1,
                  "k01Base": 0.03,
                  "thresholdOx": 0.4,
                  "maxCellNumber": 1.0,
                  "initialConditionFile": None,
                  "coexecutionFile": None,
                  "coexecutionWait": 0.0}
    
    parameters = ChicCommon.parse_input_file(inputFile, parameters)
    
    # Set up the simulation
    outputDirectory = os.getcwd()
    simulation = CellSimulation(outputDirectory, parameters, inputFile)
    
    # If there are initial conditions apply them
    if simulation.initialConditionFile:
        inputList = {"Factor":[],
                    "Vessel Density":[]}
        
        inputList = {"G1 Cells": [],
                     "S Cells": [],
                     "G2 Cells": [],
                     "M Cells": [],
                     "G0 Cells": [],
                     "Oxygen Concentration": []}
        
        simulation.get_input(simulation.initialConditionFile, inputList)
    
    # If the simulation is a single run perform the run
    if not simulation.coexecutionFile:
        simulation.run()
    else:
        # watch the diffusible component file for updates and step through 
        # the simulation as they appear
        event_handler = ChicCommon.SimulationHandler(patterns = [simulation.coexecutionFile])
        event_handler.add_simulation(simulation, simulation.coexecutionFile)
        
        observer = Observer()
        observer.schedule(event_handler, path='.', recursive=False)
        observer.start()
    
        try:
            counter = 0.0
            while counter <= simulation.coexecutionWait:
                time.sleep(1)
                counter = counter + 1.0
            observer.stop()
        except KeyboardInterrupt:
            observer.stop()
        observer.join()