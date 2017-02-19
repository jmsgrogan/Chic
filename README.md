# CHIC Hypermodel - Demo Cell and Vessel Components

## Installation (Ubuntu 14.10)
This is only tested on Ubuntu 14.10. If you are using a different OS try running on a Ubuntu virtual machine.

First install [Chaste](https://chaste.cs.ox.ac.uk/trac/wiki/InstallGuides/UbuntuPackage). Do:

```bash
sudo gedit /etc/apt/sources.list.d/chaste.list
deb http://www.cs.ox.ac.uk/chaste/ubuntu utopic/
sudo apt-key adv --recv-keys --keyserver hkp://keyserver.ubuntu.com:80 422C4D99
sudo apt-get update
sudo apt-get install chaste-source
cd $INSTALL_DIR
tar -jxf /usr/src/chaste-source.tar.bz2
```

Next install [Muscle2](http://apps.man.poznan.pl/trac/muscle/wiki/Installation). Do:

```bash
cd $INSTALL_DIR
mkdir MUSCLE-Build
curl -#O http://apps.man.poznan.pl/trac/muscle/downloads/MUSCLE-2.1.0-sources.zip
unzip MUSCLE-2.1.0-sources.zip
cd MUSCLE-2.1.0-sources/build
sudo ./build.sh $INSTALL_DIR/MUSCLE-Build
```

Set up the MUSCLE_HOME env variable and tell Chaste where the MUSCLE headers and
libraries are. Add:

```bash
source $INSTALL_DIR/MUSCLE-Build/etc/muscle.profile
export LD_LIBRARY_PATH=$INSTALL_DIR/MUSCLE-Build/lib:$LD_LIBRARY_PATH
alias muscle2="$INSTALL_DIR/MUSCLE-Build/bin/muscle2"
```

to ~/.bashrc. Now add:

```python
    use_muscle = int(prefs.get('use-muscle', True))
    if use_muscle:
        other_includepaths.append('$INSTALL_DIR/MUSCLE-Build/include/')
        other_libpaths.append('$INSTALL_DIR/MUSCLE-Build/lib/')
        other_libraries.insert(0,'muscle2')  
```

to line 149 of the file `$INSTALL_DIR/Chaste/python/hostconfig/ubuntu.py`. Keep the 
indentation similar to the lines below.

You should also add `vtkImaging' to the `other_libraries' list just below where you have 
just added the above lines. So you will have something like:

```python
    if use_vtk:
        # Note: 10.10 uses VTK 5.4, 10.04 uses 5.2, and early use 5.0
        other_includepaths.extend(vtk_include_path)
        other_libraries.extend(['vtkImaging', 'vtkIO', 'vtkCommon', 'vtkGraphics', 'z'])
```

around lines 159-162.

Clone this repo and link the contents with Chaste:

```bash
cd $INSTALL_DIR
git clone https://jmsgrogan@bitbucket.org/jmsgrogan/chic_lc_demo.git
cd $INSTALL_DIR/Chaste/projects
ln -s ../../chic_lc_demo
```

Build the model executables:

```bash
cd $INSTALL_DIR/Chaste
scons projects/chic_lc_demo/apps cl=1 exe=1
```

cd to the executables:

```bash
cd $INSTALL_DIR/chic_lc_demo/apps/src
```

Open `chic_cell_vessel.cxa.rb' and change lines 4 and 5 to reflect where your
input file (`input_file') is and where you would like your output files to
be placed (`output_file'). For example, you could use something like:

```bash
$env['input_file'] = `$INSTALL_DIR/chic_lc_demo/apps/src/data/clinical_image_3d.vti'
$env['output_file'] = '/scratch/connor/ChicOutput/hypermodel_output_3d'
```

You will need to create the output file directory (if it does not already exist) before running the hypermodel.

You can then run an example hypermodel by doing:

```bash
muscle2 -mac chic_cell_vessel.cxa.rb
```

## Repo Overview

* src: source code directory, containing the simulation base class, cell and vessel simulations and ode models governing cell and vessel growth.
* test: test directory, tests for the ode models and standalone simulations.
* test/data input files for tests. clinical_image is an artifical vtk file with a binary mask of a tumour, it is the input for the cell component. 
vessel_standalone_image is a vtk file with cell populations. It is the input for the vessel standalone component.
* python-d6_2: old python code for the 6_2 deliverable
* apps/src: the main files for the binaries and the input file for the muscle co-execution.

## Running the models

### Standalone

There are two ways to run the standalone models. The first is using the unit testing framework. Do:

```bash
cd $INSTALL_DIR/Chaste
scons projects/chic_lc_demo/test/TestVesselSimulation.hpp cl=1 
```
to run the vessel standalone.

It can also be run using muscle....(TBC)

### Coupled

Do:

```bash
cd $INSTALL_DIR/chic_lc_demo/apps/src
muscle2 -mac chic_cell_vessel.cxa.rb
```

## Input and Output

The cell standalone takes a vtk image data file with and array 'Tumour' and values 1 in tumour regions and 0 in non-tumour regions.
It takes grid size, spacing and origin via the muscel config file. It outputs vtk image data files containing arrays 'P', 'Q', 'A', 'N'
 corresponding to cell populations at specified intervals. 

The vessel standalone takes a vtk image data file with arrays 'P', 'Q', 'N' and outputs the original arrays, along with 'Nutrient'. 
It uses the vtk image data to define gird size, spacing and locaiton.

For the coupled versions both components read grid size, spacing and location from the muscle config file. The cell component passes vector
doubles of 'P', 'Q', 'N' to the vessel component and receives the vector double 'Nutrient'. The vessel component does not have any written 
output in this case. 'Nutrient' is written out by the cell component.
