"""Tests for the vasculature component of the lung cancer hypermodel
"""

import ChicCommon
import VesselComponent
import math
import matplotlib.pyplot as plt
import numpy as np

def analytical_solution(time, parameters):
    
    Vmax = parameters[0]
    Veq = parameters[1]
    V0 = parameters[2]
    
    req = parameters[3]
    rs = parameters[4]
    
    feq = parameters[5]
    f = parameters[6]
    
    if f >= feq:
        rmax = rs
    else:
        rmax = 0.0
        
    V = (rmax * (Vmax - (Vmax - V0) * math.exp(-(rmax + req) * time)) + 
         req * (Veq - (Veq - V0) * math.exp(-(rmax + req) * time))) / (rmax + req)

    return V

def numerical_solution(timeIncrement, parameters):
    
    Vmax = parameters[0]
    Veq = parameters[1]
    V0 = parameters[2]
    
    req = parameters[3]
    rs = parameters[4]
    
    feq = parameters[5]
    f = parameters[6]
    
    if f >= feq:
        rmax = rs
    else:
        rmax = 0.0
        
    delta = (rmax * (Vmax - V0) - req * (V0 - Veq) )
       
    V = V0 + delta * timeIncrement 
    return V

Vmax = 1.0
Veq = 0.5
V0 = 0.05
req = 0.1
rs = 0.1
feq = 0.5
f = 0.6
endTime = 20
timeIncrement = 0.1

parameters = [Vmax, Veq, V0, req, rs, feq, f]

# Get the analytical solution
timeAnalytical = np.arange(0, endTime + timeIncrement, timeIncrement)
vascAnalytical = []
for eachTime in timeAnalytical:
    vascAnalytical.append(analytical_solution(eachTime, parameters))
    
# Get the numerical solution
timeNumerical = np.arange(0, endTime + timeIncrement, timeIncrement)
vascNumerical = [V0]
for eachIncrement in range(len(timeNumerical) - 1):
    newVasc = numerical_solution(timeIncrement, parameters)
    parameters[2] = newVasc
    vascNumerical.append(newVasc)
    
# Get the module solution
vesselParameters = {"name": "VesselSimulation",
              "numSteps": 100,
              "timeIncrement": timeIncrement,
              "outputFrequency": 1,
              "spacing": 0.1,
              "numX": 10,
              "numY": 10,
              "numZ": 10,
              "vMax": 1.0,
              "vEq": 0.5,
              "rMax": 0.1,
              "rEq": 0.1,
              "initialConditionFile": None,
              "coexecutionFile": None,
              "coexecutionWait": 0.0} 
  
outputDirectory = os.getcwd()
vesselSimulation = VesselComponent.VesselSimulation(outputDirectory, 
                                                vesselParameters)
for eachCell in vesselSimulation.grid.cells:
    eachCell.set_cell_population(0.05)
    
for jdx, eachCell in enumerate(vesselSimulation.grid.cells):
    eachCell.factor = 0.6
    
timeModel = [0.0]
vascModel = [0.05]
diff = [0]
for idx in range(100):
    vesselData = vesselSimulation.step(parentInput = True)
    timeModel.append((idx + 1) * timeIncrement)
    
    rsum = 0.0
    for eachEntry in vesselData:
        rsum = rsum + eachEntry[1]
    vascModel.append(rsum / len(vesselData))
    diffResult = abs((vascModel[idx + 1] - vascAnalytical[idx + 1]) / vascAnalytical[idx + 1])
    diff.append(diffResult)
    

# Get the standalone solution
#plt.plot(timeAnalytical, vascAnalytical, label='Analytical')
#plt.plot(timeNumerical, vascNumerical, label='Numerical')
plt.plot(timeModel, diff, marker='o',label='Component Model')
plt.ylabel('Abs Relative Error')
plt.xlabel('Time [T]')
legend = plt.legend(loc='upper center', shadow=True)
#plt.axis([0, endTime , 0, 0.0015])
plt.show()


    