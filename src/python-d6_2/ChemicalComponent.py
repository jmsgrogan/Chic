'''This is a demonstration component model (Hypomodel) which can 
be incorporated with other models in a composite modelling (Hypermodelling) 
framework. It is being developed as part of the CHIC project.

This is the Diffusible component. It takes in spatial densities of
cells in various stages of the cell cycle and vessel density and
returns Oxygen and Growth factor density.
'''

import dolfin as df
import sys, os
import time
from watchdog.observers import Observer
import vtk
import numpy as np

import ChicCommon

class ChemicalSimulation():
    
    def __init__(self, outputDirectory, parameters, inputFile = None):      
        self.name = parameters["name"]
        self.outputDirectory = outputDirectory
        self.coexecutionFile = parameters["coexecutionFile"]
        self.coexecutionWait = float(parameters["coexecutionWait"])
        self.initialConditionFile = parameters["initialConditionFile"]
        
        self.extents = (int(parameters["numX"]), 
                        int(parameters["numY"]), 
                        int(parameters["numZ"]))
        self.spacing = float(parameters["spacing"])        
    
        if self.extents[2] > 1:
            self.mesh = df.BoxMesh(0, 0, 0, 
                                   self.extents[0] - 1 * self.spacing, 
                                   self.extents[1] - 1 * self.spacing, 
                                   self.extents[2] - 1 * self.spacing, 
                                   self.extents[0] - 1, 
                                   self.extents[1] - 1, 
                                   self.extents[2] - 1)
        else:
            self.mesh = df.RectangleMesh(0, 0, 
                                     self.extents[0] - 1 * self.spacing, 
                                     self.extents[1] - 1 * self.spacing,
                                     self.extents[0] - 1,
                                     self.extents[1] - 1)
        
        # Oxygen parameters
        self.Dc = float(parameters["Dc"])
        self.permeability = float(parameters["permeability"])
        self.vesselConc = float(parameters["vesselConc"])
        self.consumptionRate = float(parameters["consumptionRate"])
        
        # Factor Parameters
        self.Dv = float(parameters["Dv"])
        self.factorSensitvity = float(parameters["factorSensitivity"]) 
        self.decayRate = float(parameters["decayRate"]) 
        self.increment = 0
        self.sources = [[], []]
    
    def run(self, species = "Oxygen", increment = 0):
            
        # Specify the test/trial function space
        
        df.set_log_active(False)
        V = df.FunctionSpace(self.mesh, "Lagrange", 1)
        
         
        # Specify the boundary conditions, Dirichlet on all domain faces
        def u0_boundary(x, on_boundary):
            return on_boundary
         
        # Define the problem
        u = df.TrialFunction(V)
        v = df.TestFunction(V)
        if species == "Oxygen":
            bc = df.DirichletBC(V, df.Constant(1), u0_boundary)
            f = df.Constant(0.0)
            a = self.Dc * df.inner(df.nabla_grad(u), df.nabla_grad(v))*df.dx
            L = f*v*df.dx
        elif species == "Factor":
            bc = df.DirichletBC(V, df.Constant(0), u0_boundary)
            f = df.Constant(0)
            a = (self.decayRate * u * v + self.Dv * df.inner(df.nabla_grad(u), df.nabla_grad(v)))*df.dx
            L = f*v*df.dx
         
        # Assemble the system
        A, b = df.assemble_system(a, L, bc)
         
        if species == "Oxygen":
        
            # Add vessel source terms
            vesselSources = self.sources[1]
            for eachSource in vesselSources:
                location = [point - self.spacing/2.0 for point in eachSource[0]]
                if self.extents[2] > 1:
                    delta = df.PointSource(V, df.Point(location[0], 
                                                       location[1], 
                                                       location[2]), 
                                           self.permeability * eachSource[1])
                else:
                    delta = df.PointSource(V, df.Point(location[0], 
                                                       location[1]), 
                                           self.permeability * eachSource[1])                  
                try:
                    delta.apply(b)
                except:
                    pass
                
            # Add cell sink terms
            cellSources = self.sources[0]
            for eachSource in cellSources:
                location = [point - self.spacing/2.0 for point in eachSource[0]]
                if self.extents[2] > 1:
                    delta = df.PointSource(V, df.Point(location[0], 
                                                       location[1], 
                                                       location[2]), 
                                           -self.consumptionRate * eachSource[1])
                else:
                    delta = df.PointSource(V, df.Point(location[0], 
                                                       location[1]), 
                                           -self.consumptionRate * eachSource[1])                   
                try:
                    delta.apply(b)
                except:
                    pass
        
        elif species == "Factor":     
            # Add cell source terms
            cellSources = self.sources[0]
            for eachSource in cellSources:
                location = [point - self.spacing/2.0 for point in eachSource[0]]
                if self.extents[2] > 1:
                    delta = df.PointSource(V, df.Point(location[0], 
                                                       location[1], 
                                                       location[2]), 
                                           self.factorSensitvity * eachSource[1])
                else:
                    delta = df.PointSource(V, df.Point(location[0], 
                                                       location[1]), 
                                           self.factorSensitvity * eachSource[1])               
                try:
                    delta.apply(b)   
                except:
                    pass          
                                      
        # Set up solution vector
        u = df.Function(V)
        U = u.vector()
        
        # Set up and run solver
        solver = df.KrylovSolver("cg", "ilu")
        solver.solve(A, U, b)  
        
        self.result = []
        for eachEntry in self.sources[0]:
            location = [point - self.spacing/2.0 for point in eachEntry[0]]
            if self.extents[2] <= 1:
                location = location[:2]           
            try:
                result = u(location)
            except:
                if species == "Oxygen":
                    result = 1.0
                else:
                    result = 0.0
            self.result.append(result)
            
        self.write_output(species, increment)
        return self.result
    
    def write_output(self, species, increment):
        
        data = vtk.vtkImageData()
        data.SetDimensions(self.extents) 
        data.SetOrigin(0, 0, 0)
        spacing = self.spacing
        data.SetSpacing(spacing, spacing, spacing)

        result = vtk.vtkDoubleArray()
        result.SetNumberOfComponents(1)
        result.SetNumberOfTuples(len(self.result))
        if species == "Oxygen":
            result.SetName("Oxygen Concentration")
        elif species == "Factor":
            result.SetName("Factor")
            
        # VTK uses [X, Y, Z] ordering so need to reverse the list
        # of cells.
        numCells = len(self.result)
        for idx in range(numCells):
            eachValue = self.result[idx]
            result.SetValue(idx, eachValue)
            
        data.GetPointData().AddArray(result)
        
        writer = vtk.vtkXMLImageDataWriter()    
        writer.SetFileName(self.outputDirectory + "/" + str(self.name) + species + 
                          str(increment) + ".vti")
        if vtk.VTK_MAJOR_VERSION <= 5:
            writer.SetInput(data)
        else:
            writer.SetInputData(data)
        writer.Write()   
        
    def step(self, fileName):
        
        inputList = {"CellNumber": [],
             "VesselNumber": []}
        
        self.read_input(fileName, inputList)
        
        self.sources[0] = inputList["CellNumber"]
        self.sources[1] = inputList["VesselNumber"]
        self.increment = self.increment + 1
        
    def read_input(self, fileName, inputList):
        reader = vtk.vtkXMLImageDataReader()
        reader.SetFileName(fileName)
        reader.Update()
        data = reader.GetOutput()  
        
        pointData = data.GetPointData()
        
        for key in inputList:
            eachList = []
        
            for idx in range(pointData.GetNumberOfArrays()):
                if pointData.GetArrayName(idx) == key:
                    for jdx in range(pointData.GetArray(idx).GetSize()):
                        entry = pointData.GetArray(idx).GetTuple(jdx)[0]
                        eachList.append[entry]
                    
            inputList[key] = eachList
        return inputList 
        

if __name__ == "__main__":
    
    # Parse the simulation inputs
    inputFile = ChicCommon.parse_command_line_args(sys.argv[1:], "ChemicalComponent")
    
    parameters = {"name": "ChemicalSimulation",
                  "spacing": 1.0,
                  "numX": 1,
                  "numY": 1,
                  "numZ": 1,
                  "Dc": 1.0,
                  "Dv": 1.0,
                  "permeability": 1.0,
                  "vesselConc": 1.0,
                  "consumptionRate": 1.0,
                  "factorSensitivity": 0.1,
                  "decayRate": 0.1,
                  "initialConditionFile": None,
                  "coexecutionFile": None,
                  "coexecutionWait": 0.0}
    
    parameters = ChicCommon.parse_input_file(inputFile, parameters)
    
    # Set up the simulation
    outputDirectory = os.getcwd()
    simulation = ChemicalSimulation(outputDirectory, parameters, inputFile)
    
    # If there are initial conditions apply them
#     if simulation.initialConditionFile:
#         inputList = {"Factor":[],
#                     "Vessel Density":[]}
#         
#         simulation.get_input(simulation.initialConditionFile, inputList)
    
    # If the simulation is a single run perform the run
    if not simulation.coexecutionFile:
        simulation.run("Oxygen")
        simulation.run("Factor")
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
        
    
    
    
