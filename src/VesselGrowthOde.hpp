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

#ifndef VESSELGROWTHODE_HPP_
#define VESSELGROWTHODE_HPP_

#include "AbstractOdeSystem.hpp"
#include "OdeSystemInformation.hpp"

/*
 * ODE system representing the change in vascular volume fraction
 * from an equilibrium value as a function of a stimulus.
 */

class VesselGrowthOde : public AbstractOdeSystem
{
    double mVmax;
    double mVeq;
    double mR0;
    double mR1;

public:
    VesselGrowthOde() : AbstractOdeSystem(1),
        mVmax(1.0),
        mVeq(0.5),
        mR0(0.01),
        mR1(0.01)
    {
        mpSystemInfo = OdeSystemInformation<VesselGrowthOde>::Instance();
    }

    void EvaluateYDerivatives(double time, const std::vector<double>& rY,
                              std::vector<double>& rDY)
    {
        rDY[0] = mR0 * (mVmax - rY[0]) - mR1 * (rY[0] - mVeq);
    }

    void SetParameterValues(double vMax = 1.0, double vEq = 0.5, double r0 = 0.01, double r1 = 0.01)
    {
        mVmax = vMax;
        mVeq = vEq;
        mR0 = r0;
        mR1 = r1;
    }
};

template<>
void OdeSystemInformation<VesselGrowthOde>::Initialise()
{
    this->mVariableNames.push_back("V");
    this->mVariableUnits.push_back("dimensionless");
    this->mInitialConditions.push_back(0.5);
    this->mInitialised = true;
}

#endif /*VESSELGROWTHODE_HPP_*/
