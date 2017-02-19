/*

 Copyright (c) 2005-2017, University of Oxford.
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

#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkXMLImageDataReader.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <muscle2/cppmuscle.hpp>
#include "Exception.hpp"

#include "Simulation.hpp"

Simulation::Simulation()
    : mInputFile(),
      mOutputFile(),
      mMaxIncrements(100),
      mEndTime(100.0),
      mTargetTimeIncrement(1.0),
      mCurrentTime(0.0),
      mOutputFrequency(1),
      mGridSize(scalar_vector<unsigned>(3, 10)),
      mGridSpacing(1.0),
      mGridOrigin(zero_vector<double>(3)),
      mpVtkSolution(),
      mSolutionVectors(),
      mStandalone(true),
      mNeighbours(),
      mFileInputSpatialParameters(),
      mFileOutputSpatialParameters(),
      mMuscleInputSpatialParameters(),
      mMuscleOutputSpatialParameters()
{

}

Simulation::~Simulation()
{
}

void Simulation::SetCurrentTime(double time)
{
    mCurrentTime = time;
}

void Simulation::SetFileInputSpatialParameters(std::vector<std::string> parameters)
{
	mFileInputSpatialParameters = parameters;
}

void Simulation::SetFileOutputSpatialParameters(std::vector<std::string> parameters)
{
	mFileOutputSpatialParameters = parameters;
}

void Simulation::SetMuscleInputSpatialParameters(std::vector<std::string> parameters)
{
	mMuscleInputSpatialParameters = parameters;
}

void Simulation::SetMuscleOutputSpatialParameters(std::vector<std::string> parameters)
{
	mMuscleOutputSpatialParameters = parameters;
}

void Simulation::SetGridOrigin(double x, double y, double z)
{
    mGridOrigin[0] = x;
    mGridOrigin[1] = y;
    mGridOrigin[2] = z;
}

void Simulation::SetGridSpacing(double spacing)
{
    mGridSpacing = spacing;
}

void Simulation::SetGridSize(unsigned x, unsigned y, unsigned z)
{
    mGridSize[0] = x;
    mGridSize[1] = y;
    mGridSize[2] = z;
}

void Simulation::SetInputFile(const std::string& rInputFile)
{
    mInputFile = rInputFile;
}

void Simulation::SetOutputFile(const std::string& rOutputFile)
{
    mOutputFile = rOutputFile;
}

void Simulation::SetMaxIncrements(unsigned maxIncrements)
{
    mMaxIncrements = maxIncrements;
}

void Simulation::SetEndTime(double endTime)
{
    mEndTime = endTime;
}

void Simulation::SetTargetTimeIncrement(double targetTimeIncrement)
{
    mTargetTimeIncrement = targetTimeIncrement;
}

void Simulation::SetOutputFrequency(unsigned outputFrequency)
{
    mOutputFrequency = outputFrequency;
}

vtkSmartPointer<vtkImageData> Simulation::ReadVtk(const std::string& rFilename)
{
    if(rFilename.empty())
    {
        EXCEPTION("Input File Name Not Specified.");
    }
    else
    {
        vtkSmartPointer<vtkXMLImageDataReader> p_reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
        p_reader->SetFileName(rFilename.c_str());
        p_reader->Update();
        try
        {
            return p_reader->GetOutput();
        }
        catch(...)
        {
            EXCEPTION("Error reading input VTK file.");
        }
    }
}

void Simulation::Send()
{
    // Send data with muscle
    for(unsigned idx=0;idx<mMuscleOutputSpatialParameters.size();idx++)
    {
        muscle::env::sendDoubleVector(mMuscleOutputSpatialParameters[idx] + "_out",
                mSolutionVectors[mMuscleOutputSpatialParameters[idx]]);
    }
}

void Simulation::Receive()
{
    // Take in any data from muscle
    unsigned num_points = mGridSize[0] * mGridSize[1] *mGridSize[2];
    for(unsigned idx=0; idx<mMuscleInputSpatialParameters.size(); idx++)
    {
        std::vector<double> incoming_vector = muscle::env::receiveDoubleVector(mMuscleInputSpatialParameters[idx] + "_in");
        if(incoming_vector.size() != num_points)
        {
            EXCEPTION("Number of points in incoming vector does not match number of points in grid");
        }
        mSolutionVectors[mMuscleInputSpatialParameters[idx]]=incoming_vector;
    }
}

void Simulation::Initialize()
{
    vtkSmartPointer<vtkImageData> p_input_data;
    if(mStandalone)
    {
        // Read any spatial input data from file
        p_input_data = ReadVtk(mInputFile);

        // Set up the grid
        for(unsigned idx=0;idx<3;idx++)
        {
            mGridOrigin[idx] = p_input_data->GetOrigin()[idx];
            mGridSize[idx] = p_input_data->GetDimensions()[idx];
        }
        mGridSpacing = p_input_data->GetSpacing()[0];
    }

    // Set up the vtk solution data
    mpVtkSolution = vtkSmartPointer<vtkImageData>::New();
    mpVtkSolution->SetOrigin(mGridOrigin[0], mGridOrigin[1], mGridOrigin[2]);
    mpVtkSolution->SetSpacing(mGridSpacing, mGridSpacing, mGridSpacing);
    mpVtkSolution->SetDimensions(mGridSize[0], mGridSize[1], mGridSize[2]);

    // Populate the solution vectors
    unsigned num_points = mGridSize[0] * mGridSize[1] *mGridSize[2];
    for(unsigned idx=0; idx < mFileOutputSpatialParameters.size(); idx++)
    {
        vtkSmartPointer<vtkDoubleArray> p_point_data = vtkSmartPointer<vtkDoubleArray>::New();
        p_point_data->SetNumberOfComponents(1);
        p_point_data->SetNumberOfTuples(num_points);
        p_point_data->SetName(mFileOutputSpatialParameters[idx].c_str());
        mpVtkSolution->GetPointData()->AddArray(p_point_data);
        mSolutionVectors[mFileOutputSpatialParameters[idx]] = std::vector<double>(num_points, 0.0);
    }

    if(mStandalone)
    {
        // Set all required data based on VTK file values
        for(unsigned idx=0; idx<mFileInputSpatialParameters.size(); idx++)
        {
            if(p_input_data->GetPointData()->HasArray(mFileInputSpatialParameters[idx].c_str()))
            {
                std::vector<double> point_values(num_points);
                vtkSmartPointer<vtkDataArray> p_point_data =
                        p_input_data->GetPointData()->GetArray(mFileInputSpatialParameters[idx].c_str());
                for(unsigned jdx=0; jdx<num_points; jdx++)
                {
                    point_values[jdx] = p_point_data->GetTuple1(jdx);
                }
                mSolutionVectors[mFileInputSpatialParameters[idx]] = point_values;
            }
            else
            {
                std::vector<double> point_values(num_points, 0.0);
                mSolutionVectors[mFileInputSpatialParameters[idx]] = point_values;
            }
        }
    }
}

void Simulation::SetIsStandalone(bool standalone)
{
    mStandalone = standalone;
}

void Simulation::WriteVtk(const std::string& rFilename)
{
    if(rFilename.empty())
    {
        EXCEPTION("Output filename is empty");
    }

    double num_grid_points = mGridSize[0]*mGridSize[1]*mGridSize[2];

    // Update the vtk solution
    for(unsigned idx=0; idx< mFileOutputSpatialParameters.size(); idx++)
    {
        if(mSolutionVectors[mFileOutputSpatialParameters[idx]].size()!= num_grid_points)
        {
            EXCEPTION("Number of grid points differs from the size of the solution vector");
        }

        if(mpVtkSolution->GetPointData()->GetArray(mFileOutputSpatialParameters[idx].c_str())->GetNumberOfTuples()!=num_grid_points)
        {
            EXCEPTION("Number of grid points differs from the size of the vtk solution vector");
        }

        for(unsigned jdx=0; jdx<num_grid_points; jdx++)
        {
            double solution = mSolutionVectors[mFileOutputSpatialParameters[idx]][jdx];
            mpVtkSolution->GetPointData()->GetArray(mFileOutputSpatialParameters[idx].c_str())->SetTuple1(jdx, solution);
        }
    }

    vtkSmartPointer<vtkXMLImageDataWriter> p_image_data_writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    p_image_data_writer->SetFileName(rFilename.c_str());
    p_image_data_writer->SetInputData(mpVtkSolution);
    p_image_data_writer->Update();

    try
    {
        p_image_data_writer->Write();
    }
    catch(...)
    {
        EXCEPTION("Error writing to VTK file.");
    }
}

