'''This is a demonstration composite model (Hypermodel) developed as part of
the CHIC project. 

This script manages the execution of the Hypermodel. The Hypermodel has the
following components (Hypomodels):

* The Diffusible component
This component takes in spatial cell and vessel density distributions and
returns Oxygen and Growth Factor distributions. 

* The Tumour component
This component takes in spatial oxygen distributions and returns spatial
densities of cells at particular stages of the cell cycle or cells that
have undergone apoptosis.

* The Vasculature component
This component takes in spatial growth factor distributions and returns
spatial vessel densities.

'''
import numpy as np

import os
import CellComponent
import VesselComponent
import ChemicalComponent

################## Utility Functions ##################

def isInCircle(centre, radius, point, dim = 3):
    
    '''Return true if the input point is in or on the shape
    '''
    radiusSq = radius * radius
    distSq = sum((p[0] - p[1])**2 for p in zip(centre[:dim], point[:dim]))
    
    return distSq <= radiusSq

################## Main ##################

if __name__ == "__main__":
    
    outputDirectory = os.getcwd()
    numSteps = 80
    timeIncrement = 1.0
    outputFrequency = 5
    numX = 20.0
    numY = 20.0
    numZ = 20.0
    spacing = 1.0
    
    # Create the simulation instances
    vesselParameters = {"name": "VesselSimulation",
                  "numSteps": numSteps,
                  "timeIncrement": timeIncrement,
                  "outputFrequency": outputFrequency,
                  "spacing": spacing,
                  "numX": numX,
                  "numY": numY,
                  "numZ": numZ,
                  "vMax": 1.0,
                  "vEq": 0.5,
                  "rMax": 0.2,
                  "rEq": 0.1,
                  "initialConditionFile": None,
                  "coexecutionFile": None,
                  "coexecutionWait": 0.0}   
    
    vesselSimulation = VesselComponent.VesselSimulation(outputDirectory, 
                                                    vesselParameters)
    
    cellParameters = {"name": "CellSimulation",
                  "numSteps": numSteps,
                  "timeIncrement": timeIncrement,
                  "outputFrequency": outputFrequency,
                  "spacing": spacing,
                  "numX": numX,
                  "numY": numY,
                  "numZ": numZ,
                  "k1s": 0.2,
                  "ks2": 0.2,
                  "k2m": 0.2,
                  "km1": 0.2,
                  "km0Base": 0.2,
                  "k01Base": 0.06,
                  "thresholdOx": 0.4,
                  "maxCellNumber": 1.0,
                  "initialConditionFile": None,
                  "coexecutionFile": None,
                  "coexecutionWait": 0.0}  
    
    cellSimulation = CellComponent.CellSimulation(outputDirectory, 
                                                      cellParameters)
    
    chemicalParameters = {"name": "ChemicalSimulation",
                  "spacing": spacing,
                  "numX": numX,
                  "numY": numY,
                  "numZ": numZ,
                  "Dc": 1.0,
                  "Dv": 1.0,
                  "permeability": 0.1,
                  "vesselConc": 1.0,
                  "consumptionRate": 0.3,
                  "factorSensitivity": 0.1,
                  "decayRate": 0.1,
                  "initialConditionFile": None,
                  "coexecutionFile": None,
                  "coexecutionWait": 0.0}
    
    chemicalSimulation = ChemicalComponent.ChemicalSimulation(outputDirectory, 
                                                              chemicalParameters)
    
    # Specify initial conditions for the vessel simulation
    for eachCell in vesselSimulation.grid.cells:
        eachCell.set_cell_population(0.1)
        
    # Specify initial conditions for the cell simulation
    G1 = 0.25 
    S = 0.25 
    G2 = 0.25
    M = 0.25 
    G0 = 0.0
    
    # Uncomment to add and extra sphere of cells to the domain
#     radius = 10.0
#     centre = (numX * spacing / 3.0, 
#               numY * spacing / 3.0, 
#               numZ * spacing / 3.0)
#     
#     for eachCell in cellSimulation.grid.cells:
#         eachCell.concOx = 0.8
#         if isInCircle(centre, radius, eachCell.location, dim = 3):
#             eachCell.set_cell_population(np.array([G1, S, G2 ,M, G0]))
            
    radius = 4.0
    centre = (numX * spacing / 2.0, 
              numY * spacing / 2.0, 
              numZ * spacing / 6.0)
    
    for eachCell in cellSimulation.grid.cells:
        eachCell.concOx = 0.8
        if isInCircle(centre, radius, eachCell.location, dim = 3):
            eachCell.set_cell_population(np.array([G1, S, G2 ,M, G0]))
      
      
    ####### Standalone simulations #####
    # Uncomment the next two lines and comment the Coexecution block to run standalone models
    #cellSimulation.run()
    #vesselSimulation.run()
    ###### End Standalone module executions #####
    
    ####### Coexecution simulations #####
    # Uncomment the next ten lines and comment the Standalone block to run coexecution models
    for idx in range(numSteps):
        cellData = cellSimulation.step(parentInput = True)
        vesselData = vesselSimulation.step(parentInput = True)
        chemicalSimulation.sources = [cellData, vesselData]
        oxygenData = chemicalSimulation.run("Oxygen", idx)
        factorData = chemicalSimulation.run("Factor", idx)
        for jdx, eachCell in enumerate(cellSimulation.grid.cells):
            eachCell.concOx = oxygenData[jdx]
        for jdx, eachCell in enumerate(vesselSimulation.grid.cells):
            eachCell.factor = factorData[jdx]
    ###### End Coexecution simulations #####
    
