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

#ifndef TESTCELLSIMULATION_HPP_
#define TESTCELLSIMULATION_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include "CellSimulation.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestCellSimulation : public CxxTest::TestSuite
{

public:

    void TestRunModel()
    {
        // Specify the time steps
        double increment_size = 1.0; //hr
        unsigned num_increments = 300;

        // Set up the cell simulation
        double initial_tumour_volume = 144000.0; //mm^3
        double proliferation_rate = 150.0; // mm^3/hr

        std::string output_file = "/home/grogan/test.dat";

        CellSimulation cell_simulation;
        cell_simulation.SetParameters(proliferation_rate, initial_tumour_volume);
        cell_simulation.SetTargetTimeIncrement(increment_size);
        cell_simulation.SetMaxIncrements(num_increments);
        cell_simulation.SetOutputFile(output_file);
        cell_simulation.Run();
    }

};

#endif /*TESTVESSELSIMULATION_HPP_*/
