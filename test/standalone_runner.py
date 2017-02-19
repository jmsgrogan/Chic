from subprocess import call

output_dir = "output/standalone/hypermodel"

cell_binary_location = "../bin/CellSimulator"
vessel_binary_location = "../bin/VesselSimulator"
metabolic_binary_location = "../bin/MetabolicSimulator"

dt = 1.0
num_increments = 10

for idx in range(num_increments):
    cell_exe_string = cell_binary_location + " -input " + output_dir + "_metabolic_t_" + str(idx-1) + ".vti" + " -output " + output_dir
    cell_exe_string += " -current_time " + str(idx) + " -time_increment " + str(dt)
    
    vessel_exe_string = vessel_binary_location + " -input " + output_dir + "_cell_t_" + str(idx) + ".vti" + " -output " + output_dir
    vessel_exe_string += " -current_time " + str(idx) + " -time_increment " + str(dt)
    
    metabolic_exe_string = metabolic_binary_location + " -input " + output_dir + "_vessel_t_" + str(idx) + ".vti" + " -output " + output_dir
    metabolic_exe_string += " -current_time " + str(idx)
    
    retval = call(cell_exe_string, shell=True)  
    retval = call(vessel_exe_string, shell=True)  
    retval = call(metabolic_exe_string, shell=True)  