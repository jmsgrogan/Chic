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
#include "MetabolicSimulation.hpp"
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
        bool run_standalone = true;
        if(CommandLineArguments::Instance()->OptionExists("-standalone"))
        {
            run_standalone = CommandLineArguments::Instance()->GetBoolCorrespondingToOption("-standalone");
        }

        double current_time = 0.0;
        if(CommandLineArguments::Instance()->OptionExists("-current_time"))
        {
            current_time = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-current_time");
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

        double max_nutrient = 40.0;
        if(CommandLineArguments::Instance()->OptionExists("-max_nutrient"))
        {
            max_nutrient = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-max_nutrient");
        }

        double min_nutrient = 40.0;
        if(CommandLineArguments::Instance()->OptionExists("-min_nutrient"))
        {
            min_nutrient = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-min_nutrient");
        }

        double vasc_com_interval = 1.0;
        if(!run_standalone)
        {
            // Initialise MUSCLE2 environment
            env::init(&argc, &argv);

            vasc_com_interval = atof(cxa::get_property("oncosimulator_vasculature_time_interval").c_str());

            // Print identity of running instance
            std::cout << "Using Muscle. Kernel Name: " << muscle::cxa::kernel_name() << std::endl;

        }

        // Create a new simulation
        MetabolicSimulation simulation;
        simulation.SetInputFile(input_file_path);
        simulation.SetOutputFile(output_file_path);
        simulation.SetCurrentTime(current_time);
        simulation.SetGridSpacing(GC_spacing);
        simulation.SetGridSize(GC_size_x, GC_size_y, GC_size_z);
        simulation.SetGridOrigin(GC_origin_x, GC_origin_y, GC_origin_z);
        simulation.SetIsStandalone(run_standalone);
        simulation.SetParameters(max_nutrient, min_nutrient);
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
