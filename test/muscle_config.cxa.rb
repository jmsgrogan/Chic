# Parameters declared here as $env can be queried using Muscle library in from instances

######## Parameters for all components ##########################
$env['preparation_steps'] = 1
$env['GC_size_x'] = 20 # none
$env['GC_size_y'] = 20 # none
$env['GC_size_z'] = 20 # none
$env['GC_spacing'] = 5 # mm
$env['GC_origin_x'] = 0 # mm
$env['GC_origin_y'] = 0 # mm
$env['GC_origin_z'] = 0 # mm
$env['end_time'] = 20 # hour, (range: 1-)
$env['max_timesteps'] = 20 # none, (range: 1-)
$env['default_dt'] = 1 # hour, (range: 1-)
$env['oncosimulator_vasculature_time_interval'] = 1 # hour, (range: 1-)
$env['output_frequency'] = 1 #none, (range: 1-)

# Check Muscle environment is loaded
abort "Run 'source [MUSCLE_HOME]/etc/muscle.profile' before this script" if not ENV.has_key?('MUSCLE_HOME')

# Declare kernels (model binaries) that will be called from Muscle
cell_simulator = NativeInstance.new('CellSimulator', "../bin/CellSimulator", args: "-standalone 0 -output output/muscle/hypermodel")
vessel_simulator = NativeInstance.new('VesselSimulator', "../bin/VesselSimulator", args: "-standalone 0")
metabolic_simulator = NativeInstance.new('MetabolicSimulator', "../bin/MetabolicSimulator", args: "-standalone 0")

# configure connection scheme
cell_simulator.couple(vessel_simulator, 'proliferating_out' => 'proliferating_in')
cell_simulator.couple(vessel_simulator, 'quiescent_out' => 'quiescent_in')
cell_simulator.couple(vessel_simulator, 'apoptotic_out' => 'apoptotic_in')
cell_simulator.couple(vessel_simulator, 'necrotic_out' => 'necrotic_in')
cell_simulator.couple(vessel_simulator, 'differentiated_out' => 'differentiated_in')
cell_simulator.couple(vessel_simulator, 'tumour_out' => 'tumour_in')
vessel_simulator.couple(metabolic_simulator, 'Nutrient_out' => 'Nutrient_in')
metabolic_simulator.couple(cell_simulator, 'proliferation_rate_factor_out' => 'proliferation_rate_factor_in')