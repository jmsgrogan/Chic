#! /usr/bin/env python
'''This is a demonstration component model (Hypomodel) which can 
be incorporated with other models in a composite modelling (Hypermodelling) 
framework. It is being developed as part of the CHIC project.

This is the Vasculature component. It takes in spatial factor concentration
and returns spatial density of vessels.
'''

import sys, os
import time
from watchdog.observers import Observer

import ChicCommon

class VesselSimulation(ChicCommon.Simulation):  
    
    def add_grid(self):
        self.grid = VesselGrid(self.extents, self.spacing)
        
    def get_input(self, fileName, inputList):
        inputList = self.read_input(fileName, inputList)
        
        for idx, eachCell in enumerate(self.grid.cells):
            if inputList["Factor"]:
                eachCell.factor = inputList["Factor"][idx]
            if inputList["Vessel Density"]:
                eachCell.cellPopulation = inputList["Vessel Density"][idx]
        
    def apply_cell_constants(self, parameters):
        cellConstants = {"rMax": float(parameters["rMax"]),
                         "rEq": float(parameters["rEq"]),
                         "vMax": float(parameters["vMax"]),
                         "vEq": float(parameters["vEq"])}

        for eachCell in self.grid.cells:
            eachCell.set_constants(cellConstants)    
    
class VesselGrid(ChicCommon.Grid):
    
    def add_cells(self):
        for count in range(len(self.idx)):
            point = (self.idx[count] * self.spacing, 
                     self.idy[count] * self.spacing,
                     self.idz[count] * self.spacing)
            cell = VesselCell(location = point, index = count)
            self.cells.append(cell)    
        
class VesselCell(ChicCommon.GeometricCell):
    
    def further_initialization(self):
        self.factor = 0.0
        self.cellPopulation = 0.5
        
    def set_constants(self, constants): 
        self.r0 = constants["rMax"]
        self.r1 = constants["rEq"]
        self.Veq = constants["vEq"]
        self.Vmax = constants["vMax"]
        
    def update(self, dt):
        Vn = self.cellPopulation
        if self.factor < 0.5: 
            alpha = 0.0
        else:
            alpha = 1.0
        delta = (self.r0 * alpha * (self.Vmax - Vn) - self.r1 * (Vn - self.Veq) )
        self.cellPopulation = Vn  + dt * delta
            
    def get_output(self):
        output = {"Vessel Density": self.cellPopulation,
                "Factor": self.factor}
        return output 
    
if __name__ == "__main__":
    
    # Parse the simulation inputs
    inputFile = ChicCommon.parse_command_line_args(sys.argv[1:], "VesselComponent")
    
    parameters = {"name": "VesselSimulation",
                  "numSteps": 1,
                  "timeIncrement": 1.0,
                  "outputFrequency": 1,
                  "spacing": 1.0,
                  "numX": 1,
                  "numY": 1,
                  "numZ": 1,
                  "vMax": 1.0,
                  "vEq": 0.5,
                  "rMax": 0.1,
                  "rEq": 0.1,
                  "initialConditionFile": None,
                  "coexecutionFile": None,
                  "coexecutionWait": 0.0}
    
    parameters = ChicCommon.parse_input_file(inputFile, parameters)
    
    # Set up the simulation
    outputDirectory = os.getcwd()
    simulation = VesselSimulation(outputDirectory, parameters, inputFile)
    
    # If there are initial conditions apply them
    if simulation.initialConditionFile:
        inputList = {"Factor":[],
                    "Vessel Density":[]}
        
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
