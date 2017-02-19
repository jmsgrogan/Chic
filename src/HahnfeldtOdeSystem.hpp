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

#ifndef HAHNFELDTODESYSTEM_HPP_
#define HAHNFELDTODESYSTEM_HPP_


#include "AbstractOdeSystem.hpp"
#include "OdeSystemInformation.hpp"

/*
 * System of ODEs representing Hahnfeldt model of vascular
 * tumour growth (in absence of treatment).
 *
 * Tumour volume: V(t)
 * Carrying capacity: K(t)
 *
 * dV/dt = lambda*V*(1 - (V/K)^alpha)
 * dK/dt = -gamma*K + b*V - d*K*V^(2/3)
 *
 */

class HahnfeldtOdeSystem : public AbstractOdeSystem
{

private:

    double mLambda; // day^(-1)
    double mAlpha; // dimensionless
    double mGamma; // day^(-1)
    double mB; // mm(K)^3 day^(-1) mm(V)^(-3)
    double mD; // day^(-1) mm^(-2)

public:


    HahnfeldtOdeSystem();

    void EvaluateYDerivatives(double time, const std::vector<double>& rY,
                              std::vector<double>& rDY);

    void SetLambda(double value);

    void SetAlpha(double value);

    void SetGamma(double value);

    void SetB(double value);

    void SetD(double value);


};

#endif /* HAHNFELDTODESYSTEM_HPP_ */
