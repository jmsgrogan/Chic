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
#include <climits>
#include <cstring>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <muscle2/cppmuscle.hpp>
#include "CommandLineArguments.hpp"
#include "SmartPointers.hpp"
#include "CellSimulation.hpp"
#include "ExecutableSupport.hpp"
#include "Exception.hpp"
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

        bool run_standalone = true;
        if(CommandLineArguments::Instance()->OptionExists("-standalone"))
        {
            run_standalone = CommandLineArguments::Instance()->GetBoolCorrespondingToOption("-standalone");
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

        if(output_file_path.empty())
        {
            EXCEPTION("Output file name required. Use -output");
        }

        double GC_spacing = 5.0;
        if(CommandLineArguments::Instance()->OptionExists("-GC_spacing"))
        {
            GC_spacing = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-GC_spacing");
        }

        double GC_size_x = 20.0;
        if(CommandLineArguments::Instance()->OptionExists("-GC_size_x"))
        {
            GC_size_x = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-GC_size_x");
        }

        double GC_size_y = 20.0;
        if(CommandLineArguments::Instance()->OptionExists("-GC_size_y"))
        {
            GC_size_y = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-GC_size_y");
        }

        double GC_size_z = 20.0;
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

        double proliferation_rate = 150.0;
        if(CommandLineArguments::Instance()->OptionExists("-proliferation_rate"))
        {
            proliferation_rate = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-proliferation_rate");
        }

        double initial_volume = 144000.0;
        if(CommandLineArguments::Instance()->OptionExists("-initial_volume"))
        {
            initial_volume = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-initial_volume");
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

        double centre_x = 50.0;
        if(CommandLineArguments::Instance()->OptionExists("-centre_x"))
        {
            centre_x = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-centre_x");
        }

        double centre_y = 50.0;
        if(CommandLineArguments::Instance()->OptionExists("-centre_y"))
        {
            centre_y = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-centre_y");
        }

        double centre_z = 50.0;
        if(CommandLineArguments::Instance()->OptionExists("-centre_z"))
        {
            centre_z = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-centre_z");
        }
        c_vector<double, 3> centre;
        centre[0] = centre_x;
        centre[1] = centre_y;
        centre[2] = centre_z;

        double vasc_com_interval = 1.0;
        if(!run_standalone)
        {
            // Initialise MUSCLE2 environment
            env::init(&argc, &argv);

            max_timesteps = atoi(cxa::get_property("max_timesteps").c_str());
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

        // Create a new simulation
        CellSimulation simulation;
        simulation.SetInputFile(input_file_path);
        simulation.SetOutputFile(output_file_path);
        simulation.SetMaxIncrements(max_timesteps);
        simulation.SetTargetTimeIncrement(time_increment);
        simulation.SetGridSpacing(GC_spacing);
        simulation.SetCurrentTime(current_time);
        simulation.SetGridSize(GC_size_x, GC_size_y, GC_size_z);
        simulation.SetGridOrigin(GC_origin_x, GC_origin_y, GC_origin_z);
        simulation.SetIsStandalone(run_standalone);
        simulation.SetParameters(proliferation_rate, initial_volume, centre);
        simulation.SetOutputFrequency(vasc_com_interval);

        // Run the simulation
        simulation.Run();
        if(!run_standalone)
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
