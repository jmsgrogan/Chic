/*

 Copyright (c) 2005-2015, University of Oxford.
 All rights reserved.

 University of Oxford means the Chancellor, Masters and Scholars of the
 University of Oxford, having an administrative office at Wellington
 Square, Oxford OX1 2JD, UK.

 This file is part of Chaste.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
 contributors may be used to endorse or promote products derived from this
 software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#ifndef SIMULATION_HPP_
#define SIMULATION_HPP_

#include <vector>
#include <string>
#include <map>
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include "SmartPointers.hpp"
#include "UblasVectorInclude.hpp"

/**
 * Base simulation class with common functionality for vessel and cell simulation components.
 */
class Simulation
{

protected:

	/**
	 * The path to an input file
	 */
    std::string mInputFile;

	/**
	 * The path to an output file
	 */
    std::string mOutputFile;

	/**
	 * The maximum number of time increments
	 */
    unsigned mMaxIncrements;

	/**
	 * The simulation end time
	 */
    double mEndTime;

	/**
	 * The desired time increment
	 */
    double mTargetTimeIncrement;

	/**
	 * The current simulation time
	 */
    double mCurrentTime;

	/**
	 * The frequency of simulaion output
	 */
    unsigned mOutputFrequency;

	/**
	 * The number of grid points in each direction
	 */
    c_vector<unsigned, 3> mGridSize;

	/**
	 * The grid spacing
	 */
    double mGridSpacing;

	/**
	 * The grid origin
	 */
    c_vector<double, 3> mGridOrigin;

	/**
	 * The solution in vtk form for output
	 */
    vtkSmartPointer<vtkImageData> mpVtkSolution;

	/**
	 * The solutions keyed with a Field name
	 */
    std::map<std::string, std::vector<double> > mSolutionVectors;

	/**
	 * Should be run in standalone mode, or with MUSCLES
	 */
    bool mStandalone;

	/**
	 * Grid neighbours
	 */
    std::vector<std::vector<unsigned> > mNeighbours;

    /**
     * Spatial parameters to be read from file
     */
    std::vector<std::string> mFileInputSpatialParameters;

    /**
     * Parameters to be output to files
     */
    std::vector<std::string> mFileOutputSpatialParameters;

    /**
     * Parameters to be received in muscle
     */
    std::vector<std::string> mMuscleInputSpatialParameters;

    /**
     * Parameters to be sent through muscle
     */
    std::vector<std::string> mMuscleOutputSpatialParameters;

public:

    /**
     * Constructor.
     */
    Simulation();

    /**
     * Destructor
     */
    virtual ~Simulation();

    /**
     * Run the model
     */
    virtual void Run();

    /**
     * Set the names of spatial parameters to be read from file
     * @param parameters the spatial parameters to be read from files
     */
    void SetFileInputSpatialParameters(std::vector<std::string> parameters);

    /**
     * Set the names of spatial parameters to be written to file
     * @param parameters the spatial parameters to be written to files
     */
    void SetFileOutputSpatialParameters(std::vector<std::string> parameters);

    /**
     * Set the names of spatial parameters to be received from muscle
     * @param parameters the spatial parameters to be received from muscle
     */
    void SetMuscleInputSpatialParameters(std::vector<std::string> parameters);

    /**
     * Set the names of spatial parameters to be sent with muscle
     * @param parameters the spatial parameters to be sent with muscle
     */
    void SetMuscleOutputSpatialParameters(std::vector<std::string> parameters);

    /**
     * Set the current simulation time
     * @param time the current simulation time
     */
    void SetCurrentTime(double time);

    /**
     * Set model end time
     * @param endTime the model end time
     */
    void SetEndTime(double endTime);

    /**
     * Set the origin for the regular grid
     * @param x x coordinate
     * @param y z coordinate
     * @param y z coordinate
     */
    void SetGridOrigin(double x, double y = 0.0, double z = 0.0);

    /**
     * Set the grid spacing
     * @param spacing the grid spacing
     */
    void SetGridSpacing(double spacing);

    /**
     * Set the number of grid points in each direction
     * @param x number of grid points in x direction
     * @param y number of grid points in y direction
     * @param y number of grid points in z direction
     */
    void SetGridSize(unsigned x, unsigned y = 10, unsigned z = 10);

    /**
     * The path to an input file
     * @param rInputFile the input file path
     */
    void SetInputFile(const std::string& rInputFile);

    /**
     * Is the model to be run in standalone mode or through MUSCLE
     * @param standalone run the model in standalone modes
     */
    void SetIsStandalone(bool standalone);

    /**
     * The path to an output file
     * @param rOutputFile the output file path
     */
    void SetOutputFile(const std::string& rOutputFile);

    /**
     * Set the maximum number of time increments
     * @param maxIncrements the max number of increments
     */
    void SetMaxIncrements(unsigned maxIncrements);

    /**
     * Set target time increment
     * @param targetTimeIncrement the target time increment
     */
    void SetTargetTimeIncrement(double targetTimeIncrement);

    /**
     * Set the output frequency
     * @param outputFrequency the frequency for output
     */
    void SetOutputFrequency(unsigned outputFrequency);


protected:

    /**
     * Code to be run before calling Run()
     */
    virtual void Initialize();

    /**
     * Do a muscle send
     */
    void Send();

    /**
     * Do a muscle receive
     */
    void Receive();

    /**
     * Read a VTK file
     * @param rFilename the path to the file
     */
    vtkSmartPointer<vtkImageData> ReadVtk(const std::string& rFilename);

    /**
     * Write a VTK file
     * @param rFilename the path to the file
     * @param outputArrayNames the array names for vtk output
     */
    void WriteVtk(const std::string& rFilename);
};

#endif /*SIMULATION_HPP_*/
