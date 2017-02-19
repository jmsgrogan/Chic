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
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <boost/lexical_cast.hpp>
#include <muscle2/cppmuscle.hpp>
#include "Exception.hpp"
#include "LinearSystem.hpp"
#include "ReplicatableVector.hpp"
#include "OdeSolution.hpp"
#include "VesselGrowthOde.hpp"
#include "EulerIvpOdeSolver.hpp"

#include "VesselSimulation.hpp"

VesselSimulation::VesselSimulation() : Simulation(),
        mInitialVolumeFraction(0.25),
        mNutrientDiffusivity(1.e-6),
        mStimulusDiffusivity(1.e-6),
        mStimulusDecayRate(0.36),
        mStimulusReleaseRate(1.48),
        mNutrientConsumptionRate(1.48),
        mNutrientDeliveryRate(1.0),
        mVesselNutrientConcentration(40.0),
        mStimulusConcentrationInHealthy(0.0),
        mNutrientConcentrationInHealthy(40.0),
        mMaxVesselFraction(0.5),
        mEquilibriumVesselFraction(0.25),
        mRateOfVesselGrowth(0.1),
        mRateOfVesselRegression(0.01),
        mVesselGrowthTimstep(1.0)
{
      // Set default parameter array names
      this->mFileInputSpatialParameters.push_back("proliferating");
      this->mFileInputSpatialParameters.push_back("quiescent");
      this->mFileInputSpatialParameters.push_back("apoptotic");
      this->mFileInputSpatialParameters.push_back("necrotic");
      this->mFileInputSpatialParameters.push_back("differentiated");
      this->mFileInputSpatialParameters.push_back("tumour");

      this->mFileOutputSpatialParameters.push_back("proliferating");
      this->mFileOutputSpatialParameters.push_back("quiescent");
      this->mFileOutputSpatialParameters.push_back("differentiated");
      this->mFileOutputSpatialParameters.push_back("apoptotic");
      this->mFileOutputSpatialParameters.push_back("necrotic");
      this->mFileOutputSpatialParameters.push_back("tumour");
      this->mFileOutputSpatialParameters.push_back("vessel");
      this->mFileOutputSpatialParameters.push_back("stimulus");
      this->mFileOutputSpatialParameters.push_back("nutrient");

      this->mMuscleInputSpatialParameters.push_back("proliferating");
      this->mMuscleInputSpatialParameters.push_back("quiescent");
      this->mMuscleInputSpatialParameters.push_back("apoptotic");
      this->mMuscleInputSpatialParameters.push_back("necrotic");
      this->mMuscleInputSpatialParameters.push_back("differentiated");
      this->mMuscleInputSpatialParameters.push_back("tumour");

	  //mMuscleOutputSpatialParameters.push_back("nutrient");
}

VesselSimulation::~VesselSimulation()
{
}

void VesselSimulation::SetParameters(double initialVolumeFraction,
                                     double nutrientDiffusivity,
                                     double stimulusDiffusivity,
                                     double stimulusDecayRate,
                                     double stimulusReleaseRate,
                                     double nutrientConsumptionRate,
                                     double nutrientDeliveryRate,
                                     double vesselNutrientConcentration,
                                     double stimulusConcentrationHealthy,
                                     double nutrientConcentrationInHealthy,
                                     double maxVesselFraction,
                                     double equilibriumVesselFraction,
                                     double rateOfVesselGrowth,
                                     double rateOfVesselRegression,
                                     double vesselGrowthTimstep)
{
    mInitialVolumeFraction = initialVolumeFraction;
    mNutrientDiffusivity = nutrientDiffusivity;
    mStimulusDiffusivity = stimulusDiffusivity;
    mStimulusDecayRate = stimulusDecayRate;
    mStimulusReleaseRate = stimulusReleaseRate;
    mNutrientConsumptionRate = nutrientConsumptionRate;
    mNutrientDeliveryRate = nutrientDeliveryRate;
    mVesselNutrientConcentration = vesselNutrientConcentration;
    mStimulusConcentrationInHealthy = stimulusConcentrationHealthy;
    mNutrientConcentrationInHealthy = nutrientConcentrationInHealthy;
    mMaxVesselFraction = maxVesselFraction;
    mEquilibriumVesselFraction = equilibriumVesselFraction;
    mRateOfVesselGrowth = rateOfVesselGrowth;
    mRateOfVesselRegression = rateOfVesselRegression;
    mVesselGrowthTimstep = vesselGrowthTimstep;
}

void VesselSimulation::Initialize()
{
    // Do the base class initialization
    Simulation::Initialize();

    unsigned num_points = mGridSize[0] * mGridSize[1] * mGridSize[2];

    // Over-ride to set initial vessel volume fraction
    for (unsigned jdx = 0; jdx < num_points; jdx++)
    {
        this->mSolutionVectors["vessel"][jdx] = mInitialVolumeFraction;
    }
}

void VesselSimulation::Send()
{
    // Send the nutrient solution with muscle
    Simulation::Send();

    // Special case for nutrients
    muscle::env::sendDoubleVector("Nutrient_out", this->mSolutionVectors["nutrient"]);
}

void VesselSimulation::UpdateFields(unsigned speciesIndex)
{
	double diffusivity = 0.0;
    if(speciesIndex == 0)
    {
        diffusivity = mStimulusDiffusivity;
    }
    else
    {
        diffusivity = mNutrientDiffusivity;
    }

    unsigned number_of_points = mGridSize[0] * mGridSize[1] * mGridSize[2];

    // Set up the system
    LinearSystem linear_system(number_of_points, 7);
    unsigned grid_index;
    unsigned grid_index2;
    double diff_term = diffusivity / (mGridSpacing * mGridSpacing);

    for (unsigned i = 0; i < mGridSize[2]; i++) // Z
    {
        for (unsigned j = 0; j < mGridSize[1]; j++) // Y
        {
            for (unsigned k = 0; k < mGridSize[0]; k++) // X
            {
                grid_index = k + mGridSize[0] * j + mGridSize[0] * mGridSize[1] * i;

                if(speciesIndex == 0)
                {
                    linear_system.AddToMatrixElement(grid_index, grid_index,
                            -mStimulusDecayRate - 6.0 * diff_term);
                }
                else
                {
                    double cell_numbers = mSolutionVectors["proliferating"][grid_index] +
                            mSolutionVectors["quiescent"][grid_index] +
                            mSolutionVectors["differentiated"][grid_index];
                    double linear_term = -(mSolutionVectors["vessel"][grid_index] +
                            mNutrientConsumptionRate * (cell_numbers));
                    linear_system.AddToMatrixElement(grid_index, grid_index,
                            linear_term - 6.0 * diff_term);
                }

                // Apply boundary conditions
                // No flux at x bottom
                if (k > 0)
                {
                    grid_index2 = grid_index - 1;
                    linear_system.AddToMatrixElement(grid_index, grid_index2, diff_term);
                }
                else
                {
                    linear_system.AddToMatrixElement(grid_index, grid_index, diff_term);
                }

                // No flux at x top
                if (k < mGridSize[0] - 1)
                {
                    grid_index2 = grid_index + 1;
                    linear_system.AddToMatrixElement(grid_index, grid_index2, diff_term);
                }
                else
                {
                    linear_system.AddToMatrixElement(grid_index, grid_index, diff_term);
                }

                // No flux at y bottom
                if (j > 0)
                {
                    grid_index2 = grid_index - mGridSize[0];
                    linear_system.AddToMatrixElement(grid_index, grid_index2, diff_term);
                }
                else
                {
                    linear_system.AddToMatrixElement(grid_index, grid_index, diff_term);
                }

                // No flux at y top
                if (j < mGridSize[1] - 1)
                {
                    grid_index2 = grid_index + mGridSize[0];
                    linear_system.AddToMatrixElement(grid_index, grid_index2, diff_term);
                }
                else
                {
                    linear_system.AddToMatrixElement(grid_index, grid_index, diff_term);
                }

                // No flux at z bottom
                if (i > 0)
                {
                    grid_index2 = grid_index - mGridSize[0] * mGridSize[1];
                    linear_system.AddToMatrixElement(grid_index, grid_index2, diff_term);
                }
                else
                {
                    linear_system.AddToMatrixElement(grid_index, grid_index, diff_term);
                }

                // No flux at z top
                if (i < mGridSize[2] - 1)
                {
                    grid_index2 = grid_index + mGridSize[0] * mGridSize[1];
                    linear_system.AddToMatrixElement(grid_index, grid_index2, diff_term);
                }
                else
                {
                    linear_system.AddToMatrixElement(grid_index, grid_index, diff_term);
                }
                double constant_term;
                if(speciesIndex==0)
                {
                    constant_term = mStimulusReleaseRate * (mSolutionVectors["quiescent"][grid_index] +
                            mSolutionVectors["apoptotic"][grid_index]);
                }
                else
                {
                    constant_term = mVesselNutrientConcentration * mSolutionVectors["vessel"][grid_index];
                }

                linear_system.SetRhsVectorElement(grid_index, -constant_term);
            }
        }
    }

    // Dirichlet for non-tumour regions
    std::vector<unsigned> bc_indices;
    for (unsigned row = 0; row < number_of_points; row++)
    {
        if(mSolutionVectors["proliferating"][row] +
                mSolutionVectors["quiescent"][row] +
                mSolutionVectors["apoptotic"][row] +
                mSolutionVectors["differentiated"][row]<   1.e-3)
        {
            bc_indices.push_back(row);
            if(speciesIndex == 0)
            {
                linear_system.SetRhsVectorElement(row, mStimulusConcentrationInHealthy);
            }
            else
            {
                linear_system.SetRhsVectorElement(row, mNutrientConcentrationInHealthy);
            }
        }
    }

    linear_system.ZeroMatrixRowsWithValueOnDiagonal(bc_indices, 1.0);

    // Solve the linear system
    linear_system.AssembleFinalLinearSystem();
    ReplicatableVector soln_repl(linear_system.Solve());

    // Update the solution
    if(speciesIndex == 0)
    {
        for (unsigned row = 0; row < number_of_points; row++)
        {
            mSolutionVectors["stimulus"][row] = soln_repl[row];
        }
    }
    else
    {
        for (unsigned row = 0; row < number_of_points; row++)
        {
            mSolutionVectors["nutrient"][row] = soln_repl[row];
        }
    }
}

void VesselSimulation::Run()
{
    // Simulation main loop
    Initialize();

    // Set up the ODE
    unsigned num_points = mGridSize[0] * mGridSize[1] *mGridSize[2];
    VesselGrowthOde vessel_growth_ode;
    EulerIvpOdeSolver euler_solver;
    double r0;
    double total_time = 0.0;

    for(unsigned idx = 0; idx<mMaxIncrements; idx++)
    {
        // Communicate with the other simulators
        if(!mStandalone)
        {
            Receive();
        }

        // Update the nutrient and factor fields
        UpdateFields(0);
        UpdateFields(1);

        // Update the vessel volume fractions
        for(unsigned jdx = 0; jdx<num_points; jdx++)
        {
            double growth_stimulus =  this->mSolutionVectors["stimulus"][jdx];
            if(growth_stimulus > 0.5)
            {
                r0 = mRateOfVesselGrowth*mSolutionVectors["nutrient"][jdx];
            }
            else
            {
                r0 = 0.0;
            }
            vessel_growth_ode.SetParameterValues(mMaxVesselFraction,
            		mEquilibriumVesselFraction, r0, mRateOfVesselRegression);

            std::vector<double> initial_condition;
            initial_condition.push_back(mSolutionVectors["vessel"][jdx]);
            OdeSolution solutions = euler_solver.Solve(&vessel_growth_ode,
                                                       initial_condition,
                                                       0.0,
                                                       mTargetTimeIncrement,
                                                       mVesselGrowthTimstep,
                                                       mTargetTimeIncrement);

            unsigned num_steps = solutions.rGetTimes().size();
            mSolutionVectors["vessel"][jdx] = solutions.rGetSolutions()[num_steps-1][0];
        }

        // Write the output at the specified frequency
        if(idx % mOutputFrequency == 0 && mStandalone)
        {
            std::string output_file = mOutputFile + "_vessel_t_" + boost::lexical_cast<std::string>(mCurrentTime)+".vti";
            WriteVtk(output_file);
        }

        if(!mStandalone)
        {
            Send();
        }

        // Update the total time
        total_time += mTargetTimeIncrement;
        if(total_time >= mEndTime)
        {
            break;
        }
    }
}
