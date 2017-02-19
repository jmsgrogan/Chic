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

#ifndef CELLSIMULATION_HPP_
#define CELLSIMULATION_HPP_

#include "Simulation.hpp"

/**
 * Simple non-mechanical cell spheroid growth simulation. See "Mathematical Models of
 * Avascular Tumor Growth, Roose, Chapman and Maini" for details.
 */
class CellSimulation : public Simulation
{
    /**
     * Cell proliferation rate with plentiful nutrients, mm^3/hr
     */
    double mProliferationRate;

    /**
     * Initial tumour volume, mm^3
     */
    double mInitialVolume;

    /**
     * The current volume
     */
    double mCurrentVolume;

    /**
     * Set the centre of the tumour
     */
    c_vector<double, 3> mCentre;

public:

    /**
     * Constructor.
     */
    CellSimulation();

    /**
     * Destructor
     */
    ~CellSimulation();

    /**
     * Overridden method. Call once before running
     */
    void Initialize();

    /**
     * Set the model parameters
     * @param proliferationRate cell proliferation rate with plentiful nutrients, mm^3/hr
     * @param initialVolume initial tumour volume, mm^3
     * @param coordinates of tumour centre
     */
    void SetParameters(double proliferationRate,
                       double initialVolume,
                       c_vector<double, 3> centre = zero_vector<double>(3));

    /**
     * Run the simulation
     */
    void Run();
};

#endif /*CELLSIMULATION_HPP_*/
