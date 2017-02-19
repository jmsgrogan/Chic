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

#include "HahnfeldtModelSimulation.hpp"
#include "OdeSolution.hpp"
#include "FileFinder.hpp"
#include "CvodeAdaptor.hpp"

HahnfeldtModelSimulation::HahnfeldtModelSimulation() :
mOutputFolder(),
mInitialCarryingCapacity(1.0),
mInitialTumourVolume(0.05),
mODESystem(new HahnfeldtOdeSystem),
mStartTime(0),
mEndTime(100),
mSampleTime(1),
mRelTol(1e-4),
mAbsTol(1e-6),
mMaxSteps(500)
{

}

void HahnfeldtModelSimulation::SetOutputFolder(std::string output_file)
{
    mOutputFolder = output_file;
}

void HahnfeldtModelSimulation::SetInitialCarryingCapacity(double initial_carrying_capacity)
{
    assert(initial_carrying_capacity > 0);
    mInitialCarryingCapacity = initial_carrying_capacity;
}

void HahnfeldtModelSimulation::SetInitialTumourVolume(double initial_tumour_volume)
{
    assert(initial_tumour_volume > 0);
    mInitialTumourVolume = initial_tumour_volume;
}

void HahnfeldtModelSimulation::SetODE(boost::shared_ptr<HahnfeldtOdeSystem> odeSystem)
{
    mODESystem = odeSystem;
}

void HahnfeldtModelSimulation::SetStartTime(double time)
{
    assert(time >= 0);
    mStartTime = time;
}

void HahnfeldtModelSimulation::SetEndTime(double time)
{
    assert(time > 0);
    mEndTime = time;
}

void HahnfeldtModelSimulation::SetSolutionSampleTime(double time)
{
    assert(time > 0);
    mSampleTime = time;
}

void HahnfeldtModelSimulation::SetRelTol(double tol)
{
    assert(tol > 0);
    mRelTol = tol;
}

void HahnfeldtModelSimulation::SetAbsTol(double tol)
{
    assert(tol > 0);
    mAbsTol = tol;
}

void HahnfeldtModelSimulation::SetMaxStepsBetweenSolutionSampleTimes(long int value)
{
    assert(value > 1);
    mMaxSteps = value;
}


void HahnfeldtModelSimulation::Run()
{

    CvodeAdaptor solver(mRelTol, mAbsTol);
    solver.SetMaxSteps(mMaxSteps);

    OdeSolution solutions;

    std::vector<double> initial_conditions = OdeSystemInformation<HahnfeldtOdeSystem>::Instance()->GetInitialConditions();

    initial_conditions[0] = mInitialTumourVolume;
    initial_conditions[1] = mInitialCarryingCapacity;

    solutions = solver.Solve(mODESystem.get(), initial_conditions, mStartTime, mEndTime, mSampleTime/2.0, mSampleTime);

    if (mOutputFolder.empty())
    {
        EXCEPTION("Output folder not specified.");
    }
    else
    {
        FileFinder file_finder(mOutputFolder, RelativeTo::ChasteTestOutput);
        if (!file_finder.Exists())
        {
            EXCEPTION("Specified output folder does not exist.");
        }
        if (!file_finder.IsDir())
        {
            EXCEPTION("Specified output folder is not a directory.");
        }

        solutions.WriteToFile(mOutputFolder, "solution", "days");
    }
}
