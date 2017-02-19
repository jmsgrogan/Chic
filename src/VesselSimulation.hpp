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

#ifndef VESSELSIMULATION_HPP_
#define VESSELSIMULATION_HPP_

#include "Simulation.hpp"

/**
 * Vessel component for Chic Updates nutrient and growth factor fields and
 * vessel volume fraction.
 */
class VesselSimulation : public Simulation
{
    /**
     * Initial vessel volume fraction
     */
    double mInitialVolumeFraction;

    /**
     * Nutrient diffusion coefficient
     */
    double mNutrientDiffusivity;

    /**
     * Stimulus diffusion coefficient
     */
    double mStimulusDiffusivity;

    /**
     * Stimulus decay rate
     */
    double mStimulusDecayRate;

    /**
     * Stimulus release rate
     */
    double mStimulusReleaseRate;

    /**
     * Nutrient consumption rate
     */
    double mNutrientConsumptionRate;

    /**
     * Nutrient delivery rate
     */
    double mNutrientDeliveryRate;

    /**
     * Nutrient concentration in vessels
     */
    double mVesselNutrientConcentration;

    /**
     * Stimulus concentration in normal tissue
     */
    double mStimulusConcentrationInHealthy;

    /**
     * Nutrient concentration in normal tissue
     */
    double mNutrientConcentrationInHealthy;

    /**
     * Max volume fraction of vessels
     */
    double mMaxVesselFraction;

    /**
     * Equilibrium volume fraction of vessels
     */
    double mEquilibriumVesselFraction;

    /**
     * Rate of vessel growth in nutrient
     */
    double mRateOfVesselGrowth;

    /**
     * Rate of vessel regression without nutrient and stimulus
     */
    double mRateOfVesselRegression;

    /**
     * Time step for vessel growth update
     */
    double mVesselGrowthTimstep;

public:

    /**
     * Constructor.
     */
    VesselSimulation();

    /**
     * Destructor
     */
    ~VesselSimulation();

    /**
     * Set the model parameters
     */
    void SetParameters(double initialVolumeFraction,
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
                       double vesselGrowthTimstep);

    /**
     * Over-ridden model run methods
     */
    void Run();

private:

    /**
     * Over-ridden muscle send
     */
    void Send();

    /**
     * Over-ridden initialize method
     */
    void Initialize();

    /**
     * Update the solution fields
     * @param speciesIndex the index of the species to be updated
     */
    void UpdateFields(unsigned speciesIndex = 0);
};

#endif /*VESSELSIMULATION_HPP_*/
