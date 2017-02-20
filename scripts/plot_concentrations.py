import numpy as np
import matplotlib.pyplot as plt
import math
import csv
import matplotlib

def get_cons(positions, theta_1, theta_2):
    
    A = (theta_2*theta_2/(theta_1*theta_1))
    top = np.sinh(theta_1*positions)
    bottom = positions*np.sinh(theta_1)

    return (1.0-A)*(top/bottom) + A

def read_csv_vals(file):
    with open(file, 'r') as csvFile:
        csvReader = csv.reader(csvFile)
        locations = []
        concentrations = []
        for index, row in enumerate(csvReader):
            if index>0:
                locations.append((float(row[10])-120.0)/100.0)
                concentrations.append(float(row[8]))
    return locations, concentrations

radius = 100.0 # mm
diffusivity = 0.396 # mm2/hr
num_cells = 1.e9
cn = 1.0
vasc_efficiency = 0.25 #/hr

positions = np.linspace(0.001, 1, 110)
consumption_rate = 7.6e-10 # per cell/per hr
theta_2 = math.sqrt((radius*radius/diffusivity)*vasc_efficiency)
theta_1 = math.sqrt((radius*radius/diffusivity)*(num_cells*consumption_rate+vasc_efficiency))
print theta_1, theta_2
cr1 = get_cons(positions, theta_1, theta_2)

consumption_rate = 3.8e-10 # per cell/per hr
theta_2 = math.sqrt((radius*radius/diffusivity)*vasc_efficiency)
theta_1 = math.sqrt((radius*radius/diffusivity)*(num_cells*consumption_rate+vasc_efficiency))
print theta_1, theta_2
cr2 = get_cons(positions, theta_1, theta_2)

consumption_rate = 1.9e-10 # per cell/per hr
theta_2 = math.sqrt((radius*radius/diffusivity)*vasc_efficiency)
theta_1 = math.sqrt((radius*radius/diffusivity)*(num_cells*consumption_rate+vasc_efficiency))
print theta_1, theta_2
cr3 = get_cons(positions, theta_1, theta_2)

floc1, fconc1 = read_csv_vals("/home/grogan/nutrient_1.csv")
floc2, fconc2 = read_csv_vals("/home/grogan/nutrient_2.csv")
floc3, fconc3 = read_csv_vals("/home/grogan/nutrient_3.csv")

font = {'family' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)

fig, ax = plt.subplots()
ax.set_ylim([0,1])
ax.set_xlim([0,1.1])
ax.plot(positions, cr1, color="red", lw=1.5, )
ax.plot(positions, cr2, color="red", lw=1.5, )
ax.plot(positions, cr3, color="red", lw=1.5, )
ax.plot(floc1, fconc1, color="blue", linestyle="--")
ax.plot(floc2, fconc2, color="blue", linestyle="--")
ax.plot(floc3, fconc3, color="blue", linestyle="--")
#plt.axis('equal')
plt.show()
