///*
//
// Copyright (c) 2005-2015, University of Oxford.
// All rights reserved.
//
// University of Oxford means the Chancellor, Masters and Scholars of the
// University of Oxford, having an administrative office at Wellington
// Square, Oxford OX1 2JD, UK.
//
// This file is part of Chaste.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// * Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
// * Neither the name of the University of Oxford nor the names of its
// contributors may be used to endorse or promote products derived from this
// software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
// GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// */
//
//#include <tinyxml2.h>
//#ifdef CHASTE_VTK
//#define _BACKWARD_BACKWARD_WARNING_H 1 //Cut out the strstream deprecated warning for now (gcc4.3)
//#include <vtkXMLImageDataReader.h>
//#include <vtkPointData.h>
//#include <vtkDoubleArray.h>
//#endif // CHASTE_VTK
//
//#include "InputParser.hpp"
//
//InputParser::InputParser() :
//        mInputFile(),
//        mInputData(),
//        mpVtkData(),
//        mpSpatialInputData()
//{
//}
//
//InputParser::~InputParser()
//{
//}
//
//
//std::map<std::string, std::string> InputParser::GetInputData()
//{
//    return mInputData;
//}
//
//void InputParser::ParseCommandLineArgs(std::vector<std::string> arguments)
//{
//
//    for (unsigned idx = 1; idx < arguments.size(); idx++)
//    {
//        if (arguments[idx] == "-h")
//        {
//            std::cout << "Usage: CellSimulator -i <inputfile>";
//            std::cin.get();
//            exit(0);
//        }
//        else if (arguments[idx] == "-i")
//        {
//            mInputFile = arguments[idx + 1];
//        }
//    }
//
//    // Return an error message if an input file is not specified.
//    if (mInputFile.empty())
//    {
//        std::cout << "Input file not specified. \n";
//        std::cout << "Usage: CellSimulator -i <inputfile>";
//        std::cin.get();
//        exit(0);
//    }
//}
//
//void InputParser::ParseInputFile()
//{
//    // Convert the contents of the xml input file to a map
//    tinyxml2::XMLDocument doc;
//    std::map<std::string, std::string> mInputData;
//
//    const char* c = mInputFile.c_str();
//    doc.LoadFile(c);
//
//    for (std::map<std::string, std::string>::iterator it(mInputData.begin()); it != mInputData.end(); ++it)
//    {
//        const char* tag = it->first.c_str();
//        mInputData[it->first] = std::string(doc.FirstChildElement()->FirstChildElement(tag)->GetText());
//    }
//}
//
//#ifdef CHASTE_VTK
//
//vtkSmartPointer<vtkImageData> InputParser::GetInputVtkData()
//{
//    return mpVtkData;
//}
//
//void InputParser::ParseInputVtkFile(std::vector<std::string> keys)
//{
//    vtkSmartPointer<vtkXMLImageDataReader> p_reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
//        p_reader->SetFileName(mInputFile.c_str());
//        p_reader->Update();
//
//        mpVtkData = p_reader->GetOutput();
//}
//
//boost::shared_ptr<std::map<std::string, std::vector<double> > > InputParser::GetInputVtkData()
//{
//    // Set up the data map
//    boost::shared_ptr<std::map<std::string, std::vector<double> > > mpSpatialInputData =
//            boost::shared_ptr<std::map<std::string, std::vector<double> > >(new std::map<std::string, std::vector<double> >());
//    for(unsigned idx=0; idx< keys.size(); idx++)
//    {
//        *mpSpatialInputData[keys(idx)] = std::vector<double>();
//    }
//
//    vtkSmartPointer<vtkPointData> p_point_data = mpVtkData->GetPointData();
//
//    // Iterate over the input template keys, if a corresponding array name is found in the VTK file populate
//    // the corresponding data list in the map.
//    for (std::map<std::string, std::vector<double> >::iterator it(*mSpatialInputData.begin());
//            it != *mSpatialInputData.end(); ++it)
//    {
//        for (int i=0; i < p_point_data->GetNumberOfArrays(); i++)
//        {
//            if (p_point_data->GetArrayName(i) == it->first)
//            {
//                vtkDoubleArray* scalars = vtkDoubleArray::SafeDownCast(p_point_data->GetArray(i));
//                for(vtkIdType j = 0; j < scalars->GetNumberOfTuples(); j++)
//                {
//                    it->second.push_back(scalars->GetValue(j));
//                }
//            }
//        }
//    }
//
//    return mpSpatialInputData;
//}
//#endif // CHASTE_VTK
//
