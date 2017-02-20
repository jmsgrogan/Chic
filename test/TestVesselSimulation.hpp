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
 * Redistributions in binary form must reproduce the abovea copyright notice,
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

#ifndef TESTVESSELSIMULATION_HPP_
#define TESTVESSELSIMULATION_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>
#include "VesselSimulation.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestVesselSimulation : public CxxTest::TestSuite
{

public:

    void xTestStandaloneVesselSimulation2d()
    {
        // Create a new simulation
        VesselSimulation simulation;

        FileFinder file_finder("projects/Chic/apps/src/data/clinical_image_2d.vti",
                RelativeTo::ChasteSourceRoot);
        simulation.SetInputFile(file_finder.GetAbsolutePath());

        OutputFileHandler output_file_handler("TestStandaloneVesselSimulation", false);
        simulation.SetOutputFile(output_file_handler.GetOutputDirectoryFullPath() + "/vessel_sim_output_2d");
        simulation.SetMaxIncrements(3);
        simulation.SetEndTime(3);
        simulation.SetTargetTimeIncrement(1);
        simulation.SetOutputFrequency(1);

        // Run the simulation
        simulation.Run();
    }

    void TestStandaloneVesselSimulation3d()
    {
//        FileFinder file_finder("projects/Chic/apps/src/data/clinical_image_3d.vti",
//                RelativeTo::ChasteSourceRoot);
        OutputFileHandler output_file_handler("TestStandaloneVesselSimulation", false);

        // Create a new simulation
        VesselSimulation simulation;
        //simulation.SetInputFile(file_finder.GetAbsolutePath());
        simulation.SetInputFile("/home/grogan/rad_100.vti");
        simulation.SetOutputFile(output_file_handler.GetOutputDirectoryFullPath() +
                "/vessel_sim_output_3d");
        simulation.SetMaxIncrements(1);
        simulation.SetEndTime(1);
        simulation.SetTargetTimeIncrement(1);
        simulation.SetOutputFrequency(1);

        // Run the simulation
        simulation.Run();
    }

    void XTestStandaloneVesselSimulation()
    {
        // Create a new simulation
        VesselSimulation simulation;

        FileFinder file_finder("projects/Chic/apps/src/data/clinical_image_full.vti",
                RelativeTo::ChasteSourceRoot);
        simulation.SetInputFile(file_finder.GetAbsolutePath());

        OutputFileHandler output_file_handler("TestStandaloneVesselSimulation", false);
        simulation.SetOutputFile(output_file_handler.GetOutputDirectoryFullPath() + "/vessel_sim_output_full");
        simulation.SetMaxIncrements(3);
        simulation.SetEndTime(3);
        simulation.SetTargetTimeIncrement(1);
        simulation.SetOutputFrequency(1);

        // Run the simulation
        simulation.Run();
    }
};

#endif /*TESTVESSELSIMULATION_HPP_*/
