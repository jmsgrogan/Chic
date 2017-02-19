# Parameters declared here as $env can be queried using Muscle library in from instances

######## Parameters for all components ##########################
$env['preparation_steps'] = 1
$env['input_file'] = '/home/grogan/Chic/clinical_image_3d.vti'
$env['output_file'] = '/home/grogan/Chic/cell_standalone'

$env['GC_size_x'] = 25 # none
$env['GC_size_y'] = 25 # none
$env['GC_size_z'] = 25 # none
$env['GC_spacing'] = 4 # mm
$env['GC_origin_x'] = 0 # mm
$env['GC_origin_y'] = 0 # mm
$env['GC_origin_z'] = 0 # mm
$env['end_time'] = 20 # hour, (range: 1-)
$env['max_timesteps'] = 20 # none, (range: 1-)
$env['default_dt'] = 1 # hour, (range: 1-)
$env['oncosimulator_vasculature_time_interval'] = 1 # hour, (range: 1-)
$env['output_frequency'] = 1 #none, (range: 1-)

######## Parameters just for the OXFORD cell component ##########################
$env['run_standalone_cell'] = 1 # none (bool: 0, 1) 
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