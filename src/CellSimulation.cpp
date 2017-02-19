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

#include <math.h>
#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkProbeFilter.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <boost/lexical_cast.hpp>
#include <muscle2/cppmuscle.hpp>
#include "Exception.hpp"
#include "Debug.hpp"

#include "CellSimulation.hpp"

using namespace muscle;

CellSimulation::CellSimulation() : Simulation(),
        mProliferationRate(100.0),
        mInitialVolume(100000.0),
        mCurrentVolume(100000.0),
        mCentre(zero_vector<double>(3))
{
    mFileInputSpatialParameters.push_back("proliferation_rate_factor");

    mFileOutputSpatialParameters.push_back("proliferating");
    mFileOutputSpatialParameters.push_back("quiescent");
    mFileOutputSpatialParameters.push_back("differentiated");
    mFileOutputSpatialParameters.push_back("apoptotic");
    mFileOutputSpatialParameters.push_back("necrotic");
    mFileOutputSpatialParameters.push_back("tumour");
    mFileOutputSpatialParameters.push_back("nutrient");
    mFileOutputSpatialParameters.push_back("proliferation_rate_factor");

    mMuscleInputSpatialParameters.push_back("proliferation_rate_factor");

    mMuscleOutputSpatialParameters.push_back("proliferating");
    mMuscleOutputSpatialParameters.push_back("differentiated");
    mMuscleOutputSpatialParameters.push_back("quiescent");
    mMuscleOutputSpatialParameters.push_back("apoptotic");
    mMuscleOutputSpatialParameters.push_back("necrotic");
    mMuscleOutputSpatialParameters.push_back("proliferating");
    mMuscleOutputSpatialParameters.push_back("tumour");
}

CellSimulation::~CellSimulation()
{

}

void CellSimulation::SetParameters(double proliferationRate, double initialVolume, c_vector<double, 3> centre)
{
    mProliferationRate = proliferationRate;
    mInitialVolume = initialVolume;
    mCentre = centre;
}

void CellSimulation::Initialize()
{
    Simulation::Initialize();

    // Set up the initial cell populations
    unsigned num_points = mGridSize[0] * mGridSize[1] *mGridSize[2];
    for(unsigned idx=0; idx<mGridSize[2]; idx++)
    {
        for(unsigned jdx=0; jdx<mGridSize[1]; jdx++)
        {
            for(unsigned kdx=0; kdx<mGridSize[0]; kdx++)
            {
                unsigned index = kdx + jdx * mGridSize[0] + idx * mGridSize[0] * mGridSize[1];

                // If the point is in the tumour set the tumour flag
                c_vector<double,3> location;
                location[0] = kdx*mGridSpacing + mGridOrigin[0];
                location[1] = jdx*mGridSpacing + mGridOrigin[1];
                location[2] = idx*mGridSpacing + mGridOrigin[2];
                if(norm_2(location - mCentre) < cbrt(3.0*mInitialVolume/(4.0*M_PI)))
                {
                    mSolutionVectors["proliferating"][index] = 1.0;
                    mSolutionVectors["tumour"][index] = 1.0;
                }
            }
        }
    }

    // If there is an input file use it to calculate the proliferation rates
    if(mCurrentTime==0.0)
    {
        vtkSmartPointer<vtkImageData> p_input_data = vtkSmartPointer<vtkImageData>::New();
        p_input_data = ReadVtk(mInputFile);
        double average_prolif_rate_factor = 0.0;
        unsigned num_tumour = 0;
        for(unsigned jdx = 0; jdx<num_points; jdx++)
        {
            if(mSolutionVectors["tumour"][jdx]==1)
            {
                average_prolif_rate_factor += p_input_data->GetPointData()->GetArray("proliferation_rate_factor")->GetTuple1(jdx);
                num_tumour++;
            }
        }
        average_prolif_rate_factor /= double(num_tumour);
        mProliferationRate *=average_prolif_rate_factor;
    }
}

void CellSimulation::Run()
{
    // Simulation main loop
    MARK;
    Initialize();
    if(this->mOutputFile.empty())
    {
        EXCEPTION("Output file name required.");
    }
    MARK;
    std::ofstream output_file((this->mOutputFile + "_t_"+boost::lexical_cast<std::string>(mCurrentTime)+ ".dat").c_str());
    output_file << "time, volume \n";

    MARK;
    // Update the tumour radius
    unsigned counter = 0;
    while(counter <= this->mMaxIncrements)
    {
        // Communicate with the other simulators
        if(!mStandalone && counter>0)
        {
            Receive();
        }

        if (output_file.is_open())
        {
            output_file << double(counter) * this->mTargetTimeIncrement << ", " << mCurrentVolume<<" \n";
        }

        double current_radius = cbrt(3.0*mCurrentVolume/(4.0*M_PI));

        current_radius += mProliferationRate / (current_radius*current_radius) * this->mTargetTimeIncrement;

        mCurrentVolume = (4.0/3.0)*M_PI*current_radius*current_radius*current_radius;

        // Update the solution data
        for(unsigned idx=0; idx<mGridSize[2]; idx++)
        {
            for(unsigned jdx=0; jdx<mGridSize[1]; jdx++)
            {
                for(unsigned kdx=0; kdx<mGridSize[0]; kdx++)
                {
                    unsigned index = kdx + jdx * mGridSize[0] + idx * mGridSize[0] * mGridSize[1];

                    // If the point is in the tumour set the tumour flag
                    c_vector<double,3> location;
                    location[0] = kdx*mGridSpacing + mGridOrigin[0];
                    location[1] = jdx*mGridSpacing + mGridOrigin[1];
                    location[2] = idx*mGridSpacing + mGridOrigin[2];
                    if(norm_2(location - mCentre) < cbrt(3.0*mCurrentVolume/(4.0*M_PI)))
                    {
                        mSolutionVectors["proliferating"][index] = 1.e6;
                        mSolutionVectors["tumour"][index] = 1.0;
                    }
                }
            }
        }

        // Write the output at the specified frequency
        if(counter % mOutputFrequency == 0)
        {
            if(!mOutputFile.empty() && mStandalone)
            {
                std::string output_file = mOutputFile + "_cell_t_" + boost::lexical_cast<std::string>(mCurrentTime)+".vti";
                WriteVtk(output_file);
            }
            else if(!mOutputFile.empty())
            {
                std::string output_file = mOutputFile + "_cell_t_" + boost::lexical_cast<std::string>(counter)+".vti";
                WriteVtk(output_file);
            }
            else
            {
                EXCEPTION("Output file not specified.");
            }
        }
        if(!mStandalone)
        {
            Send();
        }
        counter ++;
    }
    MARK;
    output_file.close();
}
