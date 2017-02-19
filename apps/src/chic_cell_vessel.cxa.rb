# Parameters declared here as $env can be queried using Muscle library in from instances

######## Parameters for all components ##########################
$env['preparation_steps'] = 1
$env['input_file'] = './data/clinical_image_3d.vti'
$env['output_file'] = '/home/grogan/Chic/hypermodel_output_3d'

$env['GC_size_x'] = 50 # none
$env['GC_size_y'] = 50 # none
$env['GC_size_z'] = 50 # none
$env['GC_spacing'] = 2 # mm
$env['GC_origin_x'] = 0 # mm
$env['GC_origin_y'] = 0 # mm
$env['GC_origin_z'] = 0 # mm
$env['end_time'] = 20 # hour, (range: 1-)
$env['max_timesteps'] = 20 # none, (range: 1-)
$env['default_dt'] = 1 # hour, (range: 1-)
$env['oncosimulator_vasculature_time_interval'] = 1 # hour, (range: 1-)
$env['output_frequency'] = 1 #none, (range: 1-)

######## Parameters just for the vessel component ##########################
$env['run_standalone_vessel'] = 0 # none (bool: 0, 1) 
$env['nutrient_diffusivity'] = 7.2 # mm^2/hr (range: 1.e-6-)
$env['stimulus_diffusivity'] = 0.36 # mm^2/hr (range: 1.e-6-)
$env['stimulus_decay_rate'] = 0.36 # [S]^-1 hr^-1 (range: any)
$env['stimulus_release_rate'] = 1.48 # [S]^-1 hr^-1 (range: any)
$env['nutrient_consumption_rate'] = 0.001 # [N]^-1 hr^-1 (range: any)
$env['nutrient_delivery_rate'] = 0.0 # [N]^-1 hr^-1 (range: any)
$env['vessel_nutrient_concentration'] = 40.0 # [N] (range: 0-)
$env['stimulus_concentration_healthy'] = 0.0 # [S] (range: 0-)
$env['nutrient_concentration_healthy'] = 40.0 # [N] (range: 0-)
$env['initial_vessel_fraction'] = 0.25 # none (range: 0-1)
$env['max_vessel_fraction'] = 0.5 # none (range: 0-1)
$env['eq_vessel_fraction'] = 0.25 # none (range: 0-1)
$env['rate_of_vessel_growth'] = 0.1 # hr-1 (range: 0-)
$env['rate_of_vessel_regression'] = 0.01 # hr-1 (range: 0-)
$env['vessel_growth_timestep'] = 1.0 # hour (range: 1.e-6-)

######## Parameters just for the OXFORD cell component ##########################
$env['run_standalone_cell'] = 0 # none (bool: 0, 1) 
$env['proliferation_rate'] = 0.04 # [proliferating]^-1 hr^-1 (range: 0-)
$env['quiescence_threshold'] = 10.0 # [N] (range: 0-)
$env['pq_rate'] = 0.01 # [proliferating]^-1 hr^-1 (range: 0-)
$env['apoptotic_threshold'] = 0.02 # [N] (range: 0-)
$env['qa_rate'] = 0.01 # [quiescent]^-1 hr^-1 (range: 0-)
$env['recovery_threshold'] = 20.0 # [N] (range: 0-)
$env['qp_rate'] = 0.005 # [quiescent]^-1 hr^-1 (range: 0-)
$env['initial_nutrient'] = 0.005 # [N] (range: 0-)

# Check Muscle environment is loaded
abort "Run 'source [MUSCLE_HOME]/etc/muscle.profile' before this script" if not ENV.has_key?('MUSCLE_HOME')

# Declare kernels (model binaries) that will be called from Muscle
cell_simulator = NativeInstance.new('CellSimulator', "./CellSimulator" )
vessel_simulator = NativeInstance.new('VesselSimulator', "./VesselSimulator" )

# configure connection scheme
vessel_simulator.couple(cell_simulator, 'Nutrient_out' => 'Nutrient_in')
cell_simulator.couple(vessel_simulator, 'proliferating_out' => 'proliferating_in')
cell_simulator.couple(vessel_simulator, 'quiescent_out' => 'quiescent_in')
cell_simulator.couple(vessel_simulator, 'apoptotic_out' => 'apoptotic_in')
cell_simulator.couple(vessel_simulator, 'necrotic_out' => 'necrotic_in')
cell_simulator.couple(vessel_simulator, 'tumour_out' => 'tumour_in')
