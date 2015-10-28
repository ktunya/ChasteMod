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
#include <ctime>
#include <sstream>

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
#include "DetailedCellTracker.hpp"


class TestAlternativeTimesteppers : public AbstractCellBasedWithTimingsTestSuite
{
public:

    // We want to compare various timestepping methods with forward Euler here, including RK4,
    // to see whether performance improves. RK4 is more work per step than Euler, so we'll only
    // see a benefit if we allow RK to take bigger steps and use a larger absolute movement threshold.
    // Therefore the idea of this test is to run each method several times with different thresholds and
    // step sizes, to identify the point where accuracy is lost. A fairer comparison of the different 
    // methods is then possible

    void Test3dNodeBasedWithEulerStepper() throw (Exception)
    {
        
        double movThresholds[11] = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 1, 2, 3};
        // Predetermined minimum steps per hour for each mov threshold
        int stepsPerHour[11] = {590, 290, 220, 150, 130, 110, 80, 70, 40, 20, 20};


        for(int i=0; i<11; i++){

            double dt = 1.0/stepsPerHour[i];
            double thresh = movThresholds[i];

            std::cout << "dt: " << dt << " thresh:" << thresh << std::endl;

            RandomNumberGenerator::Instance()->Reseed(0);
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);

            // Create a simple 3D NodeBasedCellPopulation
            std::vector<Node<3>*> nodes;
            nodes.push_back(new Node<3>(0, false,  0.5, 0.0, 0.0));
            nodes.push_back(new Node<3>(1, false, -0.5, 0.0, 0.0));
            NodesOnlyMesh<3> mesh;
            mesh.ConstructNodesWithoutMesh(nodes, 2);

            std::vector<CellPtr> cells;
            MAKE_PTR(StemCellProliferativeType, pStemType);
            CellsGenerator<StochasticDurationCellCycleModel, 3> cellsGen;
            cellsGen.GenerateBasicRandom(cells, mesh.GetNumNodes(), pStemType);

            NodeBasedCellPopulation<3> cellPopulation(mesh, cells);
            cellPopulation.SetAbsoluteMovementThreshold(thresh);
            cellPopulation.SetDampingConstantNormal(1.0);

            // Set up cell-based simulation
            OffLatticeSimulation<3> simulator(cellPopulation, false, true, false, StepperChoice::EULER);
            std::stringstream outdir;
            outdir << "EulerStepperThresh" << thresh;
            simulator.SetOutputDirectory(outdir.str());
            simulator.SetSamplingTimestepMultiple(stepsPerHour[i]);
            simulator.SetDt(dt);
            simulator.SetEndTime(100.0);

            // Create a force law and pass it to the simulation
            MAKE_PTR(GeneralisedLinearSpringForce<3>, pLinearForce);
            pLinearForce->SetCutOffLength(1.5);
            simulator.AddForce(pLinearForce);

            // Add a movement tracker 
            MAKE_PTR_ARGS(DetailedCellTracker<3>, pTracking, (stepsPerHour[i],1));
            simulator.AddSimulationModifier(pTracking);

            // Run simulation
            int currentT = time(0);
            std::cout << "MovThreshold = " << thresh << std::endl;
            std::cout << "Sample rand start " << RandomNumberGenerator::Instance()->ranf() << std::endl; 
            simulator.Solve();
            std::cout << "Sample rand end " << RandomNumberGenerator::Instance()->ranf() << std::endl; 
            std::cout << "Time elapsed: " << time(0) - currentT << endl;

            for (unsigned i=0; i<nodes.size(); i++)
            {
                delete nodes[i];
            }
        
        }
    }


    void Test3dNodeBasedWithEulerStepperAdaptive() throw (Exception)
    {
        /*
        double movThresholds[11] = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 1, 2, 3};
        // Predetermined minimum steps per hour for each mov threshold
        int stepsPerHour[11] =     {70, 40, 25, 25, 20, 15, 15, 15, 5, 5, 5}; 


        for(int i=0; i<11; i++){

            double dt = 1.0/stepsPerHour[i];
            double thresh = movThresholds[i];

            std::cout << "dt: " << dt << " thresh:" << thresh << std::endl;

            RandomNumberGenerator::Instance()->Reseed(0);
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);

            // Create a simple 3D NodeBasedCellPopulation
            std::vector<Node<3>*> nodes;
            nodes.push_back(new Node<3>(0, false,  0.5, 0.0, 0.0));
            nodes.push_back(new Node<3>(1, false, -0.5, 0.0, 0.0));
            NodesOnlyMesh<3> mesh;
            mesh.ConstructNodesWithoutMesh(nodes, 10);

            std::vector<CellPtr> cells;
            MAKE_PTR(StemCellProliferativeType, pStemType);
            CellsGenerator<StochasticDurationCellCycleModel, 3> cellsGen;
            cellsGen.GenerateBasicRandom(cells, mesh.GetNumNodes(), pStemType);

            NodeBasedCellPopulation<3> cellPopulation(mesh, cells);
            cellPopulation.SetAbsoluteMovementThreshold(thresh);

            // Set up cell-based simulation
            OffLatticeSimulation<3> simulator(cellPopulation, false, true, true, StepperChoice::EULER);
            std::stringstream outdir;
            outdir << "AdaptiveEulerStepperThresh" << thresh;
            simulator.SetOutputDirectory(outdir.str());
            simulator.SetSamplingTimestepMultiple(stepsPerHour[i]);
            simulator.SetDt(dt);
            simulator.SetEndTime(100.0);

            // Create a force law and pass it to the simulation
            MAKE_PTR(GeneralisedLinearSpringForce<3>, pLinearForce);
            pLinearForce->SetCutOffLength(1.5);
            simulator.AddForce(pLinearForce);

            // Add a movement tracker 
            MAKE_PTR_ARGS(DetailedCellTracker<3>, pTracking, (stepsPerHour[i],1));
            simulator.AddSimulationModifier(pTracking);

            // Run simulation
            int currentT = time(0);
            simulator.Solve();
            std::cout << "Time elapsed: " << time(0) - currentT << " for movThreshold = " << thresh << std::endl;

            for (unsigned i=0; i<nodes.size(); i++)
            {
                delete nodes[i];
            }
        
        }*/
    }


    void Test3dNodeBasedWithRKStepper() throw (Exception)
    {
        
        double movThresholds[11] = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 1, 2, 3};
        // Predetermined minimum steps per hour for each mov threshold
        int stepsPerHour[11] = {590, 290, 220, 150, 130, 110, 80, 70, 40, 20, 20}; 


        for(int i=0; i<11; i++){

            double dt = 1.0/stepsPerHour[i];
            double thresh = movThresholds[i];

            std::cout << "dt: " << dt << " thresh:" << thresh << std::endl;

            RandomNumberGenerator::Instance()->Reseed(0);
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);

            // Create a simple 3D NodeBasedCellPopulation
            std::vector<Node<3>*> nodes;
            nodes.push_back(new Node<3>(0, false,  0.5, 0.0, 0.0));
            nodes.push_back(new Node<3>(1, false, -0.5, 0.0, 0.0));
            NodesOnlyMesh<3> mesh;
            mesh.ConstructNodesWithoutMesh(nodes, 2);

            std::vector<CellPtr> cells;
            MAKE_PTR(StemCellProliferativeType, pStemType);
            CellsGenerator<StochasticDurationCellCycleModel, 3> cellsGen;
            cellsGen.GenerateBasicRandom(cells, mesh.GetNumNodes(), pStemType);

            NodeBasedCellPopulation<3> cellPopulation(mesh, cells);
            cellPopulation.SetAbsoluteMovementThreshold(thresh);
            cellPopulation.SetDampingConstantNormal(1.0);

            // Set up cell-based simulation
            OffLatticeSimulation<3> simulator(cellPopulation, false, true, false, StepperChoice::RK);
            std::stringstream outdir;
            outdir << "RKStepperThresh" << thresh;
            simulator.SetOutputDirectory(outdir.str());
            simulator.SetSamplingTimestepMultiple(stepsPerHour[i]);
            simulator.SetDt(dt);
            simulator.SetEndTime(100.0);

            // Create a force law and pass it to the simulation
            MAKE_PTR(GeneralisedLinearSpringForce<3>, pLinearForce);
            pLinearForce->SetCutOffLength(1.5);
            simulator.AddForce(pLinearForce);

            // Add a movement tracker 
            MAKE_PTR_ARGS(DetailedCellTracker<3>, pTracking, (stepsPerHour[i],1));
            simulator.AddSimulationModifier(pTracking);

            // Run simulation
            int currentT = time(0);
            std::cout << "MovThreshold = " << thresh << std::endl;
            std::cout << "Sample rand start " << RandomNumberGenerator::Instance()->ranf() << std::endl; 
            simulator.Solve();
            std::cout << "Sample rand end " << RandomNumberGenerator::Instance()->ranf() << std::endl; 
            std::cout << "Time elapsed: " << time(0) - currentT << endl;

            for (unsigned i=0; i<nodes.size(); i++)
            {
                delete nodes[i];
            }     
        } 
    }


    void Test3dNodeBasedWithRKStepperAdaptive() throw (Exception)
    {
        /*
        double movThresholds[11] = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 1, 2, 3};
        // Predetermined minimum steps per hour for each mov threshold
        int stepsPerHour[11] = {70, 40, 25, 25, 20, 15, 15, 15, 5, 5, 5}; 


        for(int i=0; i<11; i++){

            double dt = 1.0/stepsPerHour[i];
            double thresh = movThresholds[i];

            std::cout << "dt: " << dt << " thresh:" << thresh << std::endl;

            RandomNumberGenerator::Instance()->Reseed(0);
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);

            // Create a simple 3D NodeBasedCellPopulation
            std::vector<Node<3>*> nodes;
            nodes.push_back(new Node<3>(0, false,  0.5, 0.0, 0.0));
            nodes.push_back(new Node<3>(1, false, -0.5, 0.0, 0.0));
            NodesOnlyMesh<3> mesh;
            mesh.ConstructNodesWithoutMesh(nodes, 2);

            std::vector<CellPtr> cells;
            MAKE_PTR(StemCellProliferativeType, pStemType);
            CellsGenerator<StochasticDurationCellCycleModel, 3> cellsGen;
            cellsGen.GenerateBasicRandom(cells, mesh.GetNumNodes(), pStemType);

            NodeBasedCellPopulation<3> cellPopulation(mesh, cells);
            cellPopulation.SetAbsoluteMovementThreshold(thresh);

            // Set up cell-based simulation
            OffLatticeSimulation<3> simulator(cellPopulation, false, true, true, StepperChoice::RK);
            std::stringstream outdir;
            outdir << "AdaptiveRKStepperThresh" << thresh;
            simulator.SetOutputDirectory(outdir.str());
            simulator.SetSamplingTimestepMultiple(stepsPerHour[i]);
            simulator.SetDt(dt);
            simulator.SetEndTime(100.0);

            // Create a force law and pass it to the simulation
            MAKE_PTR(GeneralisedLinearSpringForce<3>, pLinearForce);
            pLinearForce->SetCutOffLength(1.5);
            simulator.AddForce(pLinearForce);

            // Add a movement tracker 
            MAKE_PTR_ARGS(DetailedCellTracker<3>, pTracking, (stepsPerHour[i],1));
            simulator.AddSimulationModifier(pTracking);

            // Run simulation
            int currentT = time(0);
            simulator.Solve();
            std::cout << "Time elapsed: " << time(0) - currentT << " for movThreshold = " << thresh << std::endl;

            for (unsigned i=0; i<nodes.size(); i++)
            {
                delete nodes[i];
            }     
        }*/
    };
 
};

#endif /*TEST3DNODEBASEDWITHRKSTEPPER_HPP_*/
