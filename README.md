# CHIC Hypermodel Collection
This is a collection of instructions for building 'Hypomodels' with Chaste for use in the [Computational Horizons in Cancer - Chic](http://chic-vph.eu/) project. 

## Repo Overview

* `src`: source code directory, containing the simulation base class, cell and vessel simulations and ode models governing cell and vessel growth.
* `test`: test directory, tests for the ode models and standalone simulations.
* `test/data` input files for tests. clinical_image is an artifical vtk file with a binary mask of a tumour, it is the input for the cell component. vessel_standalone_image is a vtk file with cell populations. It is the input for the vessel standalone component.
* `python-d6_2`: old python code for the 6_2 deliverable
* `apps/src`: the main files for the binaries and the input file for the muscle co-execution.

## Build From Source (Linux)

### Install Chaste

First install [Chaste](https://chaste.cs.ox.ac.uk/trac/wiki/InstallGuides/UbuntuPackage). Do:

```bash
sudo gedit /etc/apt/sources.list.d/chaste.list
deb http://www.cs.ox.ac.uk/chaste/ubuntu utopic/
sudo apt-key adv --recv-keys --keyserver hkp://keyserver.ubuntu.com:80 422C4D99
sudo apt-get update
sudo apt-get install chaste-source
cd $WORK_DIR
tar -jxf /usr/src/chaste-source.tar.bz2 Chaste
```

### Install Muscle 2

Next install [Muscle2](http://apps.man.poznan.pl/trac/muscle/wiki/Installation). Do:

```bash
cd $WORK_DIR
git clone https://github.com/psnc-apps/muscle2.git Muscle
mkdir Muscle-Build
cd Muscle/build
sudo ./build.sh $WORK_DIR/Muscle-Build
```

Add Muscle 2 to the environment. Add:

```bash
source $WORK_DIR/Muscle-Build/etc/muscle.profile
export LD_LIBRARY_PATH=$WORK_DIR/Muscle-Build/lib:$LD_LIBRARY_PATH
export PATH=$WORK_DIR/Muscle-Build/bin/:$PATH
```

to `~/.bashrc`

### Build Sample Hypermodels

Clone this repo and link the contents with Chaste:

```bash
cd $WORK_DIR
git clone https://github.com/jmsgrogan/Chic.git Chic
cd $WORK_DIR/Chaste/projects
ln -s ../../Chic
```

Build the model executables:

```bash
cd $WORK_DIR
mkdir Chaste_Build
cd Chaste_Build
cmake ../Chaste -DMUSCLE_DIR=$WORK_DIR
make project_Chic -j $NUM_CPUS
```

The exectuables will be built here:

```bash
cd $WORK_DIR/Chaste_Build/projects/Chic/apps/src
```

## Running Hypermodels with Muscle

Open `$WORK_DIR/Chic/apps/src/chic_cell_vessel.cxa.rb' and change lines 4 and 5 to reflect where your
input file (`input_file') is and where you would like your output files to
be placed (`output_file'). For example, you could use something like:

```bash
$env['input_file'] = `$WORK_DIR/Chic/apps/src/data/clinical_image_3d.vti'
$env['output_file'] = '$WORK_DIR/TestOutput/hypermodel_output_3d'
```

You will need to create the output file directory (if it does not already exist) before running the hypermodel.

You can then run an example hypermodel by doing:

```bash
muscle2 -mac chic_cell_vessel.cxa.rb
```

## Running Standalone Tests

There are two ways to run the standalone models. The first is using the unit testing framework. Do:

```bash
cd $WORK_DIR/Chaste-Build
ctest -L project_Chic -j $NUM_CPUS
```

## Input and Output

The cell standalone takes a vtk image data file with and array 'tumour' and values 1 in tumour regions and 0 in non-tumour regions. It takes grid size, spacing and origin via the muscel config file. It outputs vtk image data files containing arrays `'P`, `Q`, `A`, `N` corresponding to cell populations at specified intervals. 

The vessel standalone takes a vtk image data file with arrays `P`, `Q`, `N` and outputs the original arrays, along with 'Nutrient'. It uses the vtk image data to define gird size, spacing and locaiton.

For the coupled versions both components read grid size, spacing and location from the muscle config file. The cell component passes vector doubles of `P`, `Q`, `N` to the vessel component and receives the vector double `Nutrient`. The vessel component does not have any written  output in this case. `Nutrient` is written out by the cell component.

## Building Distributables (Developer Only)
The current development version of Muscle should be used while an old a GCC version as you are
supporting should be used for the build. The Muscle library and usual GCC related libraries will still be dynamically linked.

### Download and Build Static Boost
```bash
wget http://downloads.sourceforge.net/project/boost/boost/1.54.0/boost_1_54_0.tar.gz
tar -xvf boost_1_54_0.tar.gz
mv boost_1_54_0 Boost
rm boost_1_54_0.tar.gz
mkdir Boost-Install
cd Boost
./bootstrap.sh --prefix=$WORK_DIR/Boost-Install --with-libraries=system,filesystem,serialization,program_options
./b2 install link=static
```

### Download and Build Static VTK
```bash
git clone git://vtk.org/VTK.git VTK
cd VTK
git checkout tags/v5.10.1
cd ..
mkdir VTK-Build
cd VTK-Build
ccmake ../VTK
```
Turn off shared libs and set to release build. Line below sets legacy GLX.

```bash
cmake -DCMAKE_C_FLAGS=-DGLX_GLXEXT_LEGACY -DCMAKE_CXX_FLAGS=-DGLX_GLXEXT_LEGACY -Wno-dev ../VTK
make
make install
```

### Download and Build Static PETSc

```bash
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.7.5.tar.gz
mv petsc-3.7.5 PETSc
export PETSC_DIR=$PWD/PETSc
export PETSC_ARCH=linux-gnu-opt
./config/configure.py --download-f2cblaslapack --download-mpich --download-hdf5 --download-parmetis --download-metis --with-x=false --with-clanguage=cxx --with-gnu-compilers --with-debugging=0 --with-shared-libraries=0  --with-ssl=0 --with-pthread=0
make PETSC_DIR=$PETSC_DIR PETSC_ARCH=$PETSCH_ARCH all
```

### Build Static Exectables

```bash
cd $WORK_DIR
mkdir Chaste_Build
cd Chaste_Build
ccmake ../Chaste -DMUSCLE_DIR=$WORK_DIR
```

point to static rather than system versions for Boost, VTK, PETSc. Turn on static Chaste build.

```bash
make project_Chic -j $NUM_CPUS
```

The binaries in `$WORK_DIR/Chaste_Build/projects/Chic/apps/src` will now be statically linked, except for the Muscle library and the usual GCC related libraries. 
