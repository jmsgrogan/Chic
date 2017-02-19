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

#ifndef TESTHAHNFELDTSIMULATION_HPP_
#define TESTHAHNFELDTSIMULATION_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>
#include "HahnfeldtModelSimulation.hpp"
#include "HahnfeldtOdeSystem.hpp"
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp"
#include <boost/shared_ptr.hpp>

class TestHahnfeldtSimulation : public CxxTest::TestSuite
{

public:

    void TestSimulation()
    {
        // Create a new simulation
        HahnfeldtModelSimulation simulation;

        double lambda = 0.05; // day^(-1)
        double alpha = 1; // dimensionless
        double gamma = 0.008; // day^(-1)
        double b = 0.08; // mm(K)^3 day^(-1) mm(V)^(-3)
        double d = 0.01; // day^(-1) mm^(-2)

        boost::shared_ptr<HahnfeldtOdeSystem> odeSystem(new HahnfeldtOdeSystem);
        odeSystem->SetLambda(lambda);
        odeSystem->SetAlpha(alpha);
        odeSystem->SetGamma(gamma);
        odeSystem->SetB(b);
        odeSystem->SetD(d);

        simulation.SetODE(odeSystem);

        OutputFileHandler output_file_handler("TestHahnfeldtSimulation", false);

        simulation.SetOutputFolder(output_file_handler.GetRelativePath());
        simulation.SetInitialCarryingCapacity(0.5);
        simulation.SetInitialTumourVolume(0.05);
        simulation.SetStartTime(0);
        simulation.SetEndTime(500);
        simulation.SetSolutionSampleTime(1);

        try
        {

        // Run the simulation
        simulation.Run();

        }
        catch(Exception& e)
        {
            std::cout << e.GetMessage() << std::endl;
        }
    }

};

#endif /*TESTHAHNFELDTSIMULATION_HPP_*/
