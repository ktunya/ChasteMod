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

#ifndef TESTNODEBASEDWITHALTERNATIVETIMESTEPPERS_HPP_
#define TESTNODEBASEDWITHALTERNATIVETIMESTEPPERS_HPP_

#include <cxxtest/TestSuite.h>
#include <ctime>
#include <sstream>

#include "CellBasedSimulationArchiver.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "LogFile.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "RandomNumberGenerator.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "SimpleCellCentrePositionTracker.hpp"
#include "SimplePopulationExtentTracker.hpp"

class TestNodeBasedWithAlternativeTimesteppers : public AbstractCellBasedWithTimingsTestSuite
{

public:
    
    enum RunChoices {EULER, RK4, BACKWARDEULER, ALL};
    static const int toRun = EULER;

    // We want to compare various timestepping methods with forward Euler here, including RK4, and implicit
    // methods, to see whether performance improves. RK4 etc are more work per step than forward Euler, so 
    // we'll only see a benefit if we allow RK to take larger steps and use a larger absolute movement threshold.
    // Therefore the idea of this test is to run the same simulation multiple times for each method, using different
    // thresholds and step sizes, to identify the point where accuracy is lost. A fairer comparison of the different 
    // numerical methods is then possible. See bottom of file for helper methods and test settings.

    void Test3dNodeBasedWithEulerStepper() throw (Exception)
    {
        if(toRun == EULER || toRun == ALL){

            setupTestParameters();

            for(int i = minStepIndex; i <= maxStepIndex; i++){

                double movementThresh = movThresholds[i];
                int stepsPerHour = stepsPerHour_DoesNotExceedMovThresh[i];
                double dt = 1.0/(double)stepsPerHour;

                resetForNewRun();

                std::vector<Node<3>*> nodes;
                nodes.push_back(new Node<3>(0, false,  0.5, 0.0, 0.0));
                nodes.push_back(new Node<3>(1, false, -0.5, 0.0, 0.0));
            
                NodesOnlyMesh<3> mesh;
                mesh.ConstructNodesWithoutMesh(nodes, 3.0);

                std::vector<CellPtr> cells;
                MAKE_PTR(StemCellProliferativeType, pStemType);
                CellsGenerator<StochasticDurationCellCycleModel, 3> cellsGen;
                cellsGen.GenerateBasicRandom(cells, mesh.GetNumNodes(), pStemType);

                NodeBasedCellPopulation<3>* cellPopulation = new NodeBasedCellPopulation<3>(mesh, cells);
                cellPopulation->SetAbsoluteMovementThreshold(movementThresh);
                cellPopulation->SetDampingConstantNormal(1.1);

                OffLatticeSimulation<3> simulation((*cellPopulation), false, true, false, StepperChoice::EULER);

                std::ostringstream outDir;
                outDir << "NB_Euler_Thresh" << movementThresh;

                setupSimulation(&simulation, outDir.str(), stepsPerHour);

                // Add a cell movement tracker 
                MAKE_PTR_ARGS(SimpleCellCentrePositionTracker<3>, pTracking, (stepsPerHour,1));
                simulation.AddSimulationModifier(pTracking);

                int startT = time(0);
                std::cout << "Dt = " << dt << " Absolute Movement Threshold = " << movementThresh << std::endl;
                simulation.Solve();
                std::cout << "Time elapsed: " << time(0) - startT << endl; 

                //if(i==5){ checkAgreementAtFifthStepSize(outDir.str(), pTracking->GetOutputDirectoryFull()); }

                delete nodes[0];
                delete nodes[1];
            }
        }    
    }


    void Test3dNodeBasedWithEulerStepperAdaptive() throw (Exception)
    {
        if(toRun == EULER || toRun == ALL){
            setupTestParameters();

            for(int i = minStepIndex; i <= maxStepIndex; i++){

                double movementThresh = movThresholds[i];
                int stepsPerHour = stepsPerHour_ExceedsMovThresh[i];
                double dt = 1.0/(double)stepsPerHour;

                resetForNewRun();

                std::vector<Node<3>*> nodes;
                nodes.push_back(new Node<3>(0, false,  0.5, 0.0, 0.0));
                nodes.push_back(new Node<3>(1, false, -0.5, 0.0, 0.0));
        
                NodesOnlyMesh<3> mesh;
                mesh.ConstructNodesWithoutMesh(nodes, 3.0);

                std::vector<CellPtr> cells;
                MAKE_PTR(StemCellProliferativeType, pStemType);
                CellsGenerator<StochasticDurationCellCycleModel, 3> cellsGen;
                cellsGen.GenerateBasicRandom(cells, mesh.GetNumNodes(), pStemType);

                NodeBasedCellPopulation<3>* cellPopulation = new NodeBasedCellPopulation<3>(mesh, cells);
                cellPopulation->SetAbsoluteMovementThreshold(movementThresh);
                cellPopulation->SetDampingConstantNormal(1.1);

                OffLatticeSimulation<3> simulation(*cellPopulation, false, true,  true, StepperChoice::EULER);
                std::ostringstream outDir;
                outDir << "NB_AdaptiveEuler_Thresh" << movementThresh;
                setupSimulation(&simulation, outDir.str(), stepsPerHour);

                // Add a cell movement tracker 
                MAKE_PTR_ARGS(SimpleCellCentrePositionTracker<3>, pTracking, (stepsPerHour,1));
                simulation.AddSimulationModifier(pTracking);

                int startT = time(0);
                std::cout << "Dt = " << dt << " Absolute Movement Threshold = " << movementThresh << std::endl;
                simulation.Solve();
                std::cout << "Time elapsed: " << time(0) - startT << endl; 

                //if(i==5){ checkAgreementAtFifthStepSize(outDir.str(), pTracking->GetOutputDirectoryFull()); }

                delete nodes[0];
                delete nodes[1];
            }  
        }  
    }

    void Test3dNodeBasedWithRK4Stepper() throw (Exception)
    {
        if(toRun == RK4 || toRun == ALL){
           
            setupTestParameters();

            for(int i = minStepIndex; i <= maxStepIndex; i++){

                double movementThresh = movThresholds[i];
                int stepsPerHour = stepsPerHour_DoesNotExceedMovThresh[i];
                double dt = 1.0/(double)stepsPerHour;

                resetForNewRun();

                std::vector<Node<3>*> nodes;
                nodes.push_back(new Node<3>(0, false,  0.5, 0.0, 0.0));
                nodes.push_back(new Node<3>(1, false, -0.5, 0.0, 0.0));
            
                NodesOnlyMesh<3> mesh;
                mesh.ConstructNodesWithoutMesh(nodes, 3.0);

                std::vector<CellPtr> cells;
                MAKE_PTR(StemCellProliferativeType, pStemType);
                CellsGenerator<StochasticDurationCellCycleModel, 3> cellsGen;
                cellsGen.GenerateBasicRandom(cells, mesh.GetNumNodes(), pStemType);

                NodeBasedCellPopulation<3>* cellPopulation = new NodeBasedCellPopulation<3>(mesh, cells);
                cellPopulation->SetAbsoluteMovementThreshold(movementThresh);
                cellPopulation->SetDampingConstantNormal(1.1);

                OffLatticeSimulation<3> simulation(*cellPopulation, false, true,  false, StepperChoice::RK4);
                std::ostringstream outDir;
                outDir << "NB_RK4_Thresh" << movementThresh;
                setupSimulation(&simulation, outDir.str(), stepsPerHour);

                // Add a cell movement tracker 
                MAKE_PTR_ARGS(SimpleCellCentrePositionTracker<3>, pTracking, (stepsPerHour,1));
                simulation.AddSimulationModifier(pTracking);

                int startT = time(0);
                std::cout << "Dt = " << dt << " Absolute Movement Threshold = " << movementThresh << std::endl;
                simulation.Solve();
                std::cout << "Time elapsed: " << time(0) - startT << endl; 

                //if(i==5){ checkAgreementAtFifthStepSize(outDir.str(), pTracking->GetOutputDirectoryFull()); }

                delete nodes[0];
                delete nodes[1];
            } 
        }   
    }


    void Test3dNodeBasedWithRK4StepperAdaptive() throw (Exception)
    {
        if(toRun == RK4 || toRun == ALL){
        
            setupTestParameters();

            for(int i = minStepIndex; i <= maxStepIndex; i++){

                double movementThresh = movThresholds[i];
                int stepsPerHour = stepsPerHour_ExceedsMovThresh[i];
                double dt = 1.0/(double)stepsPerHour;

                resetForNewRun();

                std::vector<Node<3>*> nodes;
                nodes.push_back(new Node<3>(0, false,  0.5, 0.0, 0.0));
                nodes.push_back(new Node<3>(1, false, -0.5, 0.0, 0.0));
            
                NodesOnlyMesh<3> mesh;
                mesh.ConstructNodesWithoutMesh(nodes, 3.0);

                std::vector<CellPtr> cells;
                MAKE_PTR(StemCellProliferativeType, pStemType);
                CellsGenerator<StochasticDurationCellCycleModel, 3> cellsGen;
                cellsGen.GenerateBasicRandom(cells, mesh.GetNumNodes(), pStemType);

                NodeBasedCellPopulation<3>* cellPopulation = new NodeBasedCellPopulation<3>(mesh, cells);
                cellPopulation->SetAbsoluteMovementThreshold(movementThresh);
                cellPopulation->SetDampingConstantNormal(1.1);

                OffLatticeSimulation<3> simulation(*cellPopulation, false, true,  true, StepperChoice::RK4);
                std::ostringstream outDir;
                outDir << "NB_AdaptiveRK4_Thresh" << movementThresh;
                setupSimulation(&simulation, outDir.str(), stepsPerHour);

                // Add a cell movement tracker 
                MAKE_PTR_ARGS(SimpleCellCentrePositionTracker<3>, pTracking, (stepsPerHour,1));
                simulation.AddSimulationModifier(pTracking);

                int startT = time(0);
                std::cout << "Dt = " << dt << " Absolute Movement Threshold = " << movementThresh << std::endl;
                simulation.Solve();
                std::cout << "Time elapsed: " << time(0) - startT << endl; 

                //if(i==5){ checkAgreementAtFifthStepSize(outDir.str(), pTracking->GetOutputDirectoryFull()); }

                delete nodes[0];
                delete nodes[1];
            }
        }  
    }

    void Test3dNodeBasedWithBackwardEulerStepper() throw (Exception)
    {
        if(toRun == BACKWARDEULER || toRun == ALL){
            
            setupTestParameters();

            for(int i = minStepIndex; i <= maxStepIndex; i++){

                double movementThresh = movThresholds[i];
                int stepsPerHour = stepsPerHour_DoesNotExceedMovThresh[i];
                double dt = 1.0/(double)stepsPerHour;

                resetForNewRun();

                std::vector<Node<3>*> nodes;
                nodes.push_back(new Node<3>(0, false,  0.5, 0.0, 0.0));
                nodes.push_back(new Node<3>(1, false, -0.5, 0.0, 0.0));
            
                NodesOnlyMesh<3> mesh;
                mesh.ConstructNodesWithoutMesh(nodes, 3.0);

                std::vector<CellPtr> cells;
                MAKE_PTR(StemCellProliferativeType, pStemType);
                CellsGenerator<StochasticDurationCellCycleModel, 3> cellsGen;
                cellsGen.GenerateBasicRandom(cells, mesh.GetNumNodes(), pStemType);

                NodeBasedCellPopulation<3>* cellPopulation = new NodeBasedCellPopulation<3>(mesh, cells);
                cellPopulation->SetAbsoluteMovementThreshold(movementThresh);
                cellPopulation->SetDampingConstantNormal(1.1);

                OffLatticeSimulation<3> simulation(*cellPopulation, false, true,  false, StepperChoice::BACKWARDEULER);
                std::ostringstream outDir;
                outDir << "NB_BackwardEuler_Thresh" << movementThresh;
                setupSimulation(&simulation, outDir.str(), stepsPerHour);

                // Add a cell movement tracker 
                MAKE_PTR_ARGS(SimpleCellCentrePositionTracker<3>, pTracking, (stepsPerHour,1));
                simulation.AddSimulationModifier(pTracking);

                int startT = time(0);
                std::cout << "Dt = " << dt << " Absolute Movement Threshold = " << movementThresh << std::endl;
                simulation.Solve();
                std::cout << "Time elapsed: " << time(0) - startT << endl; 

                //if(i==5){ checkAgreementAtFifthStepSize(outDir.str(), pTracking->GetOutputDirectoryFull()); }

                delete nodes[0];
                delete nodes[1];
            }
        }    
    }


    void Test3dNodeBasedWithBackwardEulerStepperAdaptive() throw (Exception)
    {
        if(toRun == BACKWARDEULER || toRun == ALL){

            setupTestParameters();

            for(int i = minStepIndex; i <= maxStepIndex; i++){

                double movementThresh = movThresholds[i];
                int stepsPerHour = stepsPerHour_ExceedsMovThresh[i];
                double dt = 1.0/(double)stepsPerHour;

                resetForNewRun();

                std::vector<Node<3>*> nodes;
                nodes.push_back(new Node<3>(0, false,  0.5, 0.0, 0.0));
                nodes.push_back(new Node<3>(1, false, -0.5, 0.0, 0.0));
            
                NodesOnlyMesh<3> mesh;
                mesh.ConstructNodesWithoutMesh(nodes, 3.0);

                std::vector<CellPtr> cells;
                MAKE_PTR(StemCellProliferativeType, pStemType);
                CellsGenerator<StochasticDurationCellCycleModel, 3> cellsGen;
                cellsGen.GenerateBasicRandom(cells, mesh.GetNumNodes(), pStemType);

                NodeBasedCellPopulation<3>* cellPopulation = new NodeBasedCellPopulation<3>(mesh, cells);
                cellPopulation->SetAbsoluteMovementThreshold(movementThresh);
                cellPopulation->SetDampingConstantNormal(1.1);

                OffLatticeSimulation<3> simulation(*cellPopulation, false, true,  true, StepperChoice::BACKWARDEULER);
                std::ostringstream outDir;
                outDir << "NB_AdaptiveBackwardEuler_Thresh" << movementThresh;
                setupSimulation(&simulation, outDir.str(), stepsPerHour);

                // Add a cell movement tracker 
                MAKE_PTR_ARGS(SimpleCellCentrePositionTracker<3>, pTracking, (stepsPerHour,1));
                simulation.AddSimulationModifier(pTracking);

                int startT = time(0);
                std::cout << "Dt = " << dt << " Absolute Movement Threshold = " << movementThresh << std::endl;
                simulation.Solve();
                std::cout << "Time elapsed: " << time(0) - startT << endl; 

                //if(i==5){ checkAgreementAtFifthStepSize(outDir.str(), pTracking->GetOutputDirectoryFull()); }

                delete nodes[0];
                delete nodes[1];
            } 
        }   
    }



    //----------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------

    // Test settings. Altering ranges here will drastically alter the run time of the test, since
    // a low MinStepIndex will run the test simulation for very small dt. For my work I need to
    // run at a wide range of dt, for production set: minStepIndex = maxStepIndex = 5, endTime = 100.
    int minStepIndex;
    int maxStepIndex;
    double endTime;
 
    std::vector<double> movThresholds;                   
    std::vector<int> stepsPerHour_DoesNotExceedMovThresh;
    std::vector<int> stepsPerHour_ExceedsMovThresh;            


    void setupTestParameters(){ 
        
        minStepIndex = 0;
        maxStepIndex = 11;
        endTime = 100;
 
        double threshes[12] = {0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
        movThresholds = std::vector<double>(&threshes[0], &threshes[0]+12);
        int steps1[12] = {10000, 5000, 2500, 1800, 1300, 590,  295, 147, 100, 75,  60,  59 };
        stepsPerHour_DoesNotExceedMovThresh = std::vector<int>(&steps1[0], &steps1[0]+12);
        int steps2[12] = {5000,  2500, 1800, 1300, 590,  295,  147, 73,  50,  38,  30,  25 };
        stepsPerHour_ExceedsMovThresh = std::vector<int>(&steps2[0], &steps2[0]+12);
    }


    void resetForNewRun(){
        
        RandomNumberGenerator::Instance()->Reseed(0);
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
    }


    void setupSimulation(OffLatticeSimulation<3>* sim, std::string outDir, int stepsPerHour){

        double dt = 1.0/(double)stepsPerHour; 

        sim->SetOutputDirectory(outDir.c_str());
        sim->SetSamplingTimestepMultiple(stepsPerHour);
        sim->SetDt(dt);
        sim->SetEndTime(endTime);

        // Create a force law
        MAKE_PTR(GeneralisedLinearSpringForce<3>, pLinearForce);
        pLinearForce->SetCutOffLength(3.0);
        sim->AddForce(pLinearForce);

        MAKE_PTR_ARGS(SimplePopulationExtentTracker<3>, pExtentTracker, (stepsPerHour));
        sim->AddSimulationModifier(pExtentTracker);
    }


    void checkAgreementAtFifthStepSize(std::string outDir, std::string fullDir){

        std::string referenceData = "cell_based/test/data/TestNodeBasedWithAlternativeTimesteppers/";
        referenceData += outDir.c_str();
        referenceData += "/PositionData.txt";

        FileComparison(fullDir + "PositionData.txt", referenceData).CompareFiles();
    }
};

#endif /*TESTNODEBASEDWITHALTERNATIVETIMESTEPPERS_HPP_*/
