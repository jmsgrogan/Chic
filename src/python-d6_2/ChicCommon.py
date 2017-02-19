'''Common classes and functions for CHIC component Hypermodels
'''

import sys, getopt
import numpy as np
import vtk
from xml.dom import minidom
from watchdog.events import PatternMatchingEventHandler

################# Functions #################
def parse_command_line_args(args, componentName = None):
    
    inputfile = ""
    
    try:
        opts, args = getopt.getopt(args,"hi:")
    except getopt.GetoptError:
        print 'Usage: ModelName.py -i <inputfile>'
        sys.exit(2)
        
    for opt, arg in opts:
        if opt == '-h':
            print 'Usage: ModelName.py -i <inputfile>'
            sys.exit()
        elif opt in ("-i"):
            inputfile = arg
            
    if inputfile == "" and componentName is not None:
        inputfile = componentName + ".xml"
        
    return inputfile

def parse_input_file(inputFile, parameters):
    
    try:
        xmldoc = minidom.parse(inputFile)
    except:
        print 'Error parsing input file.'
        sys.exit()    
    
    documentElement
    
    for key in parameters:
        value = None
        tag = xmldoc.getElementsByTagName(key) 
        if tag[0].childNodes:
            value = tag[0].childNodes[0].nodeValue
        parameters[key] = value
    
    return parameters

################## Classes ##################
class SimulationHandler(PatternMatchingEventHandler):
    
    def on_modified(self, event):
        self.simulation.step(self.fileName)

    def add_simulation(self, simulation, fileName = None):
        self.simulation = simulation
        self.fileName = fileName
        
class Simulation():
    
    '''Class for managing simulations
    '''   
    
    def __init__(self, outputDirectory, parameters, child = False):
        
        self.name = parameters["name"]
        self.outputDirectory = outputDirectory
        self.numSteps = int(parameters["numSteps"])
        self.timeIncrement = float(parameters["timeIncrement"])
        self.outputFrequency = int(parameters["outputFrequency"])
        self.coexecutionFile = parameters["coexecutionFile"]
        self.coexecutionWait = float(parameters["coexecutionWait"])
        self.initialConditionFile = parameters["initialConditionFile"]
        
        self.extents = (int(parameters["numX"]), 
                        int(parameters["numY"]), 
                        int(parameters["numZ"]))
        self.spacing = float(parameters["spacing"])
        self.add_grid()
        self.apply_cell_constants(parameters)  
        self.increment = 0
        
    def add_grid(self):
        self.grid = Grid(self.extents, self.spacing)
        
    def apply_cell_constants(self, parameters):
        cellConstants = {}

        for eachCell in self.grid.cells:
            eachCell.set_constants(cellConstants)
            
    def get_inputs(self, fileName):
        pass
        
    def run(self):
        
        self.write_output(fileName = self.name + "0")
        for idx in range(self.numSteps):
            
            # Cycle the cells
            for eachCell in self.grid.cells:
                eachCell.update(self.timeIncrement)
            
            # Write output
            if idx%self.outputFrequency == 0:
                self.write_output(fileName = self.name + str(idx + 1)) 
        
    def step(self, parentInput = None, fileName = None):
        
        if parentInput:
            result = [] 
        
        if fileName:
            self.read_input(fileName)  
        
        # Cycle the cells
        for eachCell in self.grid.cells:
            eachCell.update(self.timeIncrement)
            
            if parentInput:
                result.append([eachCell.location, 
                               eachCell.cellPopulation])
        # Write output
        if self.increment%self.outputFrequency == 0:
            self.write_output(fileName = self.name + str(self.increment + 1))       
        
        self.increment = self.increment + 1
        
        if parentInput:
            return result     
        
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
            
    def write_output(self, fileName):
        
        '''Write the cell data to a vtk file. 
        '''
        data = vtk.vtkImageData()
        data.SetDimensions(self.grid.extents) 
        data.SetOrigin(0, 0, 0)
        spacing = self.grid.spacing
        data.SetSpacing(spacing, spacing, spacing)

        dataArrays = []
        for eachOutput in self.grid.cells[0].get_output(): 
            array = vtk.vtkDoubleArray()
            array.SetNumberOfComponents(1)
            array.SetNumberOfTuples(self.grid.numCells)
            array.SetName(eachOutput)
            dataArrays.append(array)
            
        # VTK uses [X, Y, Z] ordering so need to reverse the list
        # of cells.
        numCells = len(self.grid.cells)
        for idx in range(numCells):
            eachCell = self.grid.cells[idx]
            output = eachCell.get_output()
            for eachArray in dataArrays:
                eachArray.SetValue(idx, output[eachArray.GetName()])
            
        for eachArray in dataArrays: 
            data.GetPointData().AddArray(eachArray)
        
        writer = vtk.vtkXMLImageDataWriter()    
        writer.SetFileName(self.outputDirectory + "/" + fileName + ".vti")
        if vtk.VTK_MAJOR_VERSION <= 5:
            writer.SetInput(data)
        else:
            writer.SetInputData(data)
        writer.Write()   
        
        
class Grid():
    
    '''Class for domain level storage and operations.
    
    This class stores the definition of the computational
    grid and a list of geometric cell objects at each grid
    point. Methods involving inter-cell operations are
    included in this class.
    '''
    
    def __init__(self, extents = (10, 10, 10), spacing = 1.0):
        
        ''' On initialization a grid is created and filled with
        cells.
        '''
        
        self.extents = extents
        self.spacing = spacing
        self.cells = []
        
        nx, ny, nz = self.extents
        self.numCells = nx * ny * nz
        idx, idy, idz = np.unravel_index(np.arange(self.numCells), 
                                 dims = self.extents, order = "F")
        
        self.idx = idx
        self.idy = idy
        self.idz = idz
        self.add_cells()
        
    def add_cells(self):
        
        for count in range(len(self.idx)):
            point = (self.idx[count] * self.spacing, 
                     self.idy[count] * self.spacing,
                     self.idz[count] * self.spacing)
            cell = GeometricCell(location = point, index = count)
            self.cells.append(cell)      
        
class GeometricCell():
    
    def __init__(self, location, index = -1):
        self.location = location
        self.index = index 
        self.further_initialization()
        
    def further_initilization(self):
        pass
                  
    def set_cell_population(self, cellPopulation = None):
        self.cellPopulation = cellPopulation
        
    def set_constants(self, constants):
        self.constant = constants
        
    def update(self, dt):
        pass
            
    def get_output(self):    
        pass