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

#ifndef TESTVESSELGROWTHODE_HPP_
#define TESTVESSELGROWTHODE_HPP_

#include <cxxtest/TestSuite.h>
#include <vector>
#include "VesselGrowthOde.hpp"
#include "EulerIvpOdeSolver.hpp"
#include "OdeSolution.hpp"

class TestVesselGrowthOde : public CxxTest::TestSuite
{

public:

    void TestVesselGrowthOdeHighFactor()
    {
        double v_max = 1.0;
        double v_eq = 0.5;
        double r0 = 0.2;
        double r1 = 0.1;

        VesselGrowthOde vessel_growth_ode;
        EulerIvpOdeSolver euler_solver;
        vessel_growth_ode.SetParameterValues(v_max, v_eq, r0, r1);

        std::vector<double> initial_condition;
        initial_condition.push_back(0.5);
        OdeSolution solutions = euler_solver.Solve(&vessel_growth_ode, initial_condition, 0.0, 1.0, 0.1, 0.1);

        for(unsigned idx=0; idx < solutions.rGetTimes().size(); idx++)
        {
            std::cout << "idx:" << idx << " V:" << solutions.rGetSolutions()[idx][0] << std::endl;
        }
    }

    void TestVesselGrowthOdeLowFactor()
    {
        double v_max = 1.0;
        double v_eq = 0.5;
        double r0 = 0.0;
        double r1 = 0.1;

        VesselGrowthOde vessel_growth_ode;
        EulerIvpOdeSolver euler_solver;
        vessel_growth_ode.SetParameterValues(v_max, v_eq, r0, r1);

        std::vector<double> initial_condition;
        initial_condition.push_back(0.75);
        OdeSolution solutions = euler_solver.Solve(&vessel_growth_ode, initial_condition, 0.0, 1.0, 0.1, 0.1);

        for(unsigned idx=0; idx < solutions.rGetTimes().size(); idx++)
        {
            std::cout << "idx:" << idx << " V:" << solutions.rGetSolutions()[idx][0] << std::endl;
        }
    }
};

#endif /*TESTVESSELGROWTHODE_HPP_*/
