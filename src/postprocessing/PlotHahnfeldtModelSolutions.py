#! /usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import os


# plot species concentrations over time 
#######################################

# read in file contatining results for a normal cell

filepath = os.environ.get('CHASTE_TEST_OUTPUT') + '/TestHahnfeldtSimulation/solution.dat'

file = open(filepath,'r')

# extract data from file

variableNames = []
values = {}

for index,line in enumerate(file.readlines()):
    if index == 0:
        variableNames = line.split(' ')
        variableNames[-1] = variableNames[-1].replace('\n','')
        for i in range(len(variableNames)):
            variableNames[i] = variableNames[i].replace('_',' ')
            values[variableNames[i]] = []
    else:
        valarray = line.split(' ')
        vals = []
        
        for index2, val in enumerate(valarray):
            if 'e' in val:
                vals.append(float(val))
                
        for i in range(len(vals) - 1):
            values[variableNames[i]].append(vals[i])
        values[variableNames[-1]].append(vals[-1])
        

file.close()

# plot results

fig,ax = plt.subplots(figsize=(15,10))

ax.hold(True)

mark = []
mark.append('x')
mark.append('o')

varNames = []

for i in range(len(variableNames) - 1):   
    multiTime = []
    multiVars = []
    for ind,t in enumerate(values['Time(days)']):
       multiTime.append(t)
       multiVars.append(values[variableNames[i+1]][ind])
    varNames.append(variableNames[i+1].replace('(mm3)','').replace('V','V - Tumour').replace('K','K - Carrying capacity'))
    
    ax.plot(multiTime ,multiVars,marker=mark[i],markevery=50,markerfacecolor='none',markersize=10)


ax.hold(False)

plt.xlabel('Time [ $days$ ]',fontsize=20)
plt.ylabel('Volume [ $mm^3$ ]',fontsize=20)
ax.tick_params(axis='both',which='major',labelsize=20)
plt.legend(varNames,loc='best')

plt.show()

