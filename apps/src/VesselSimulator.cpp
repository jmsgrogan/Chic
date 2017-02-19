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

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <muscle2/cppmuscle.hpp>
#include "VesselSimulation.hpp"
#include "ExecutableSupport.hpp"
#include "Exception.hpp"
#include "CommandLineArguments.hpp"
#include "PetscTools.hpp"
#include "PetscException.hpp"

using namespace muscle;

int main(int argc, char *argv[])
{
    ExecutableSupport::StandardStartup(&argc, &argv);

    int exit_code = ExecutableSupport::EXIT_OK;

    try
    {
        // Parse the command line options, for now use hard-coded default values if the argument is missing
        // Eventually read these values from an xml config
        unsigned max_timesteps = 1;
        if(CommandLineArguments::Instance()->OptionExists("-max_timesteps"))
        {
            max_timesteps = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-max_timesteps");
        }
        double end_time = 20.0;
        if(CommandLineArguments::Instance()->OptionExists("-end_time"))
        {
            end_time = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-end_time");
        }

        double time_increment = 1.0;
        if(CommandLineArguments::Instance()->OptionExists("-time_increment"))
        {
            time_increment = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-time_increment");
        }

        double current_time = 0.0;
        if(CommandLineArguments::Instance()->OptionExists("-current_time"))
        {
            current_time = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-current_time");
        }

        unsigned output_frequency = 1;
        if(CommandLineArguments::Instance()->OptionExists("-output_frequency"))
        {
            output_frequency = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-output_frequency");
        }

        bool run_standalone_vessel = true;
        if(CommandLineArguments::Instance()->OptionExists("-standalone"))
        {
            run_standalone_vessel = CommandLineArguments::Instance()->GetBoolCorrespondingToOption("-standalone");
        }

        std::string input_file_path;
        if(CommandLineArguments::Instance()->OptionExists("-input"))
        {
            input_file_path = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-input");
        }

        std::string output_file_path;
        if(CommandLineArguments::Instance()->OptionExists("-output"))
        {
            output_file_path = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-output");
        }

        double GC_spacing = 1.0;
        if(CommandLineArguments::Instance()->OptionExists("-GC_spacing"))
        {
            GC_spacing = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-GC_spacing");
        }

        double GC_size_x = 10.0;
        if(CommandLineArguments::Instance()->OptionExists("-GC_size_x"))
        {
            GC_size_x = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-GC_size_x");
        }

        double GC_size_y = 10.0;
        if(CommandLineArguments::Instance()->OptionExists("-GC_size_y"))
        {
            GC_size_y = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-GC_size_y");
        }

        double GC_size_z = 10.0;
        if(CommandLineArguments::Instance()->OptionExists("-GC_size_z"))
        {
            GC_size_z = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-GC_size_z");
        }

        double GC_origin_x = 0.0;
        if(CommandLineArguments::Instance()->OptionExists("-GC_origin_x"))
        {
            GC_origin_x = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("G-C_origin_x");
        }

        double GC_origin_y = 0.0;
        if(CommandLineArguments::Instance()->OptionExists("-GC_origin_y"))
        {
            GC_origin_y = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-GC_origin_y");
        }

        double GC_origin_z = 0.0;
        if(CommandLineArguments::Instance()->OptionExists("-GC_origin_z"))
        {
            GC_origin_z = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-GC_origin_z");
        }

        double initial_vessel_fraction = 0.25;
        if(CommandLineArguments::Instance()->OptionExists("-initial_vessel_fraction"))
        {
            initial_vessel_fraction = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-initial_vessel_fraction");
        }

        double nutrient_diffusivity = 7.2;
        if(CommandLineArguments::Instance()->OptionExists("-nutrient_diffusivity"))
        {
            nutrient_diffusivity = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-nutrient_diffusivity");
        }

        double stimulus_diffusivity = 0.36;
        if(CommandLineArguments::Instance()->OptionExists("-stimulus_diffusivity"))
        {
            stimulus_diffusivity = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-stimulus_diffusivity");
        }

        double stimulus_decay_rate = 0.36;
        if(CommandLineArguments::Instance()->OptionExists("-stimulus_decay_rate"))
        {
            stimulus_decay_rate = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-stimulus_decay_rate");
        }

        double stimulus_release_rate = 1.48;
        if(CommandLineArguments::Instance()->OptionExists("-stimulus_release_rate"))
        {
            stimulus_release_rate = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-stimulus_release_rate");
        }

        double nutrient_consumption_rate = 1.e-6;
        if(CommandLineArguments::Instance()->OptionExists("-nutrient_consumption_rate"))
        {
            nutrient_consumption_rate = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-nutrient_consumption_rate");
        }

        double nutrient_delivery_rate = 0.0;
        if(CommandLineArguments::Instance()->OptionExists("-nutrient_delivery_rate"))
        {
            nutrient_delivery_rate = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-nutrient_delivery_rate");
        }

        double vessel_nutrient_concentration = 40.0;
        if(CommandLineArguments::Instance()->OptionExists("-vessel_nutrient_concentration"))
        {
            vessel_nutrient_concentration = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-vessel_nutrient_concentration");
        }

        double stimulus_concentration_healthy = 0.0;
        if(CommandLineArguments::Instance()->OptionExists("-stimulus_concentration_healthy"))
        {
            stimulus_concentration_healthy = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-stimulus_concentration_healthy");
        }

        double nutrient_concentration_healthy = 40.0;
        if(CommandLineArguments::Instance()->OptionExists("-nutrient_concentration_healthy"))
        {
            nutrient_concentration_healthy = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-nutrient_concentration_healthy");
        }

        double max_vessel_fraction = 0.5;
        if(CommandLineArguments::Instance()->OptionExists("-max_vessel_fraction"))
        {
            max_vessel_fraction = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-max_vessel_fraction");
        }

        double eq_vessel_fraction = 0.25;
        if(CommandLineArguments::Instance()->OptionExists("-eq_vessel_fraction"))
        {
            eq_vessel_fraction = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-eq_vessel_fraction");
        }

        double rate_of_vessel_growth = 0.1;
        if(CommandLineArguments::Instance()->OptionExists("-rate_of_vessel_growth"))
        {
            rate_of_vessel_growth = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-rate_of_vessel_growth");
        }

        double rate_of_vessel_regression = 0.01;
        if(CommandLineArguments::Instance()->OptionExists("-rate_of_vessel_regression"))
        {
            rate_of_vessel_regression = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-rate_of_vessel_regression");
        }

        double vessel_growth_timestep = 1.0;
        if(CommandLineArguments::Instance()->OptionExists("-vessel_growth_timestep"))
        {
            vessel_growth_timestep = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-vessel_growth_timestep");
        }

        // if using muscle set parameters from muscle environment
        if(!run_standalone_vessel)
        {
            // Initialise MUSCLE2 environment
            env::init(&argc, &argv);

            max_timesteps = atoi(cxa::get_property("max_timesteps").c_str());
            end_time = atof(cxa::get_property("end_time").c_str());
            time_increment = atof(cxa::get_property("oncosimulator_vasculature_time_interval").c_str());
            GC_spacing = atof(cxa::get_property("GC_spacing").c_str());
            GC_size_x = atof(cxa::get_property("GC_size_x").c_str());
            GC_size_y = atof(cxa::get_property("GC_size_y").c_str());
            GC_size_z = atof(cxa::get_property("GC_size_z").c_str());
            GC_origin_x = atof(cxa::get_property("GC_origin_x").c_str());
            GC_origin_y = atof(cxa::get_property("GC_origin_y").c_str());
            GC_origin_z = atof(cxa::get_property("GC_origin_z").c_str());

            // Print identity of running instance
            std::cout << "Using Muscle. Kernel Name: " << muscle::cxa::kernel_name() << std::endl;

        }
        // Create a simulation instance
        VesselSimulation simulation;

        // time step data
        simulation.SetMaxIncrements(max_timesteps);
        simulation.SetEndTime(end_time);
        simulation.SetCurrentTime(current_time);
        simulation.SetTargetTimeIncrement(time_increment);
        simulation.SetOutputFrequency(output_frequency);
        if(run_standalone_vessel)
        {
            simulation.SetIsStandalone(true);
            simulation.SetInputFile(input_file_path);
            simulation.SetOutputFile(output_file_path);
        }
        else
        {
            simulation.SetGridSpacing(GC_spacing);
            simulation.SetGridSize(GC_size_x, GC_size_y, GC_size_z);
            simulation.SetGridOrigin(GC_origin_x, GC_origin_y, GC_origin_z);
            simulation.SetIsStandalone(false);
        }
        simulation.SetParameters(initial_vessel_fraction,
                                 nutrient_diffusivity,
                                 stimulus_diffusivity,
                                 stimulus_decay_rate,
                                 stimulus_release_rate,
                                 nutrient_consumption_rate,
                                 nutrient_delivery_rate,
                                 vessel_nutrient_concentration,
                                 stimulus_concentration_healthy,
                                 nutrient_concentration_healthy,
                                 max_vessel_fraction,
                                 eq_vessel_fraction,
                                 rate_of_vessel_growth,
                                 rate_of_vessel_regression,
                                 vessel_growth_timestep);

        // Run the simulation
        simulation.Run();

        // Finalise Muscle environment and cleanup
        if(!run_standalone_vessel)
        {
            env::finalize();
        }
    }
    catch (const Exception& e)
    {
        ExecutableSupport::PrintError(e.GetMessage());
        exit_code = ExecutableSupport::EXIT_ERROR;
    }

    ExecutableSupport::FinalizePetsc();
    return exit_code;
}
