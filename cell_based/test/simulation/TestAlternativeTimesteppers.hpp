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

#ifndef TEST3DNODEBASEDWITHRKSTEPPER_HPP_
#define TEST3DNODEBASEDWITHRKSTEPPER_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "LogFile.hpp"
#include "SmartPointers.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "PetscSetupAndFinalize.hpp"



class TestAlternativeTimesteppers : public AbstractCellBasedWithTimingsTestSuite
{
public:


    void Test3dNodeBasedWithEulerStepper() throw (Exception)
    {
        RandomNumberGenerator::Instance()->Reseed(0);

        // Create a simple 3D NodeBasedCellPopulation
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false,  0.5, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, -0.5, 0.0, 0.0));
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        std::vector<CellPtr> cells;
        MAKE_PTR(StemCellProliferativeType, pStemType);
        CellsGenerator<StochasticDurationCellCycleModel, 3> cellsGen;
        cellsGen.GenerateBasicRandom(cells, mesh.GetNumNodes(), pStemType);

        NodeBasedCellPopulation<3> cellPopulation(mesh, cells);
        cellPopulation.SetAbsoluteMovementThreshold(0.5);

        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator(cellPopulation, false, true, false, StepperChoice::EULER);
        simulator.SetOutputDirectory("3DNodeBasedWithEulerStepper");
        simulator.SetSamplingTimestepMultiple(250);
        simulator.SetDt(1.0/250.0);
        simulator.SetEndTime(50.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<3>, pLinearForce);
        pLinearForce->SetCutOffLength(1.5);
        simulator.AddForce(pLinearForce);

        std::cout << "Start solve" << std::endl;

        // Run simulation
        simulator.Solve();

        // Avoid memory leak
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }


    void Test3dNodeBasedWithRKStepper() throw (Exception)
    {
        RandomNumberGenerator::Instance()->Reseed(0);
    
        // Create a simple 3D NodeBasedCellPopulation
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false,  0.5, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, -0.5, 0.0, 0.0));
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);
    
        std::vector<CellPtr> cells;
        MAKE_PTR(StemCellProliferativeType, pStemType);
        CellsGenerator<StochasticDurationCellCycleModel, 3> cellsGen;
        cellsGen.GenerateBasicRandom(cells, mesh.GetNumNodes(), pStemType);
    
        NodeBasedCellPopulation<3> cellPopulation(mesh, cells);
        cellPopulation.SetAbsoluteMovementThreshold(0.5);
    
        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator(cellPopulation, false, true, false, StepperChoice::RK);
        simulator.SetOutputDirectory("3DNodeBasedWithRKStepper");
        simulator.SetSamplingTimestepMultiple(250);
        simulator.SetDt(1.0/250.0);
        simulator.SetEndTime(50.0);
    
        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<3>, pLinearForce);
        pLinearForce->SetCutOffLength(1.5);
        simulator.AddForce(pLinearForce);
    
        // Run simulation
        simulator.Solve();
    
        // Avoid memory leak
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }
 
};

#endif /*TEST3DNODEBASEDWITHRKSTEPPER_HPP_*/
