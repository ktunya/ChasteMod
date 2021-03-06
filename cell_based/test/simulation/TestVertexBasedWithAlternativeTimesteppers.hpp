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

#ifndef TESTVERTEXBASEDWITHALTERNATIVETIMESTEPPERS_HPP_
#define TESTVERTEXBASEDWITHALTERNATIVETIMESTEPPERS_HPP_

#include <cxxtest/TestSuite.h>
#include <ctime>
#include <string>

#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "LogFile.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "RandomNumberGenerator.hpp"
#include "CommandLineArguments.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NagaiHondaForce.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "VertexMeshWriter.hpp"
#include "WildTypeCellMutationState.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "SimpleCellCentrePositionTracker.hpp"
#include "SimplePopulationExtentTracker.hpp"



class TestVertexBasedWithAlternativeTimesteppers : public AbstractCellBasedWithTimingsTestSuite
{
public:

    enum RunChoices {EULER, RK4, BACKWARDEULER, ALL};
    static const int toRun = EULER;


    void Test2dVertexBasedWithEulerStepper() throw (Exception)
    {
        if(toRun == EULER || toRun == ALL){
        
            setupTestParameters();

            for(int i = minStepIndex; i <= maxStepIndex; i++){

                double movementThresh = movThresholds[i];
                int stepsPerHour = stepsPerHour_DoesNotExceedMovThresh[i];
                double dt = 1.0/(double)stepsPerHour;

                resetForNewRun();

                // Create a cell population
                HoneycombVertexMeshGenerator generator(meshSide, meshSide);
                MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

                std::vector<CellPtr> cells;
                MAKE_PTR(StemCellProliferativeType, p_stem_type);
                CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
                cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_stem_type);

                VertexBasedCellPopulation<2> cellPopulation(*p_mesh, cells);

                // Create a simulation
                OffLatticeSimulation<2> simulation(cellPopulation, false, true, false, StepperChoice::EULER);
                std::ostringstream outDir;
                outDir << "VB_Euler_Thresh" << movementThresh;
                setupSimulation(&simulation, outDir.str(), stepsPerHour);

                // Add a cell movement tracker 
                MAKE_PTR_ARGS(SimpleCellCentrePositionTracker<2>, pTracking, (stepsPerHour,1));
                simulation.AddSimulationModifier(pTracking);

                int startT = time(0);
                std::cout << "Dt = " << dt << " Absolute Movement Threshold = " << movementThresh << std::endl;
                simulation.Solve();
                std::cout << "Time elapsed: " << time(0) - startT << endl; 
            }    
        }
    }

    void Test2dVertexBasedWithEulerStepperAdaptive() throw (Exception)
    {
        if(toRun == EULER || toRun == ALL){

            setupTestParameters();

            for(int i = minStepIndex; i <= maxStepIndex; i++){

                double movementThresh = movThresholds[i];
                int stepsPerHour = stepsPerHour_ExceedsMovThresh[i];
                double dt = 1.0/(double)stepsPerHour;

                resetForNewRun();

                // Create a cell population
                HoneycombVertexMeshGenerator generator(meshSide, meshSide);
                MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

                std::vector<CellPtr> cells;
                MAKE_PTR(StemCellProliferativeType, p_stem_type);
                CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
                cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_stem_type);

                VertexBasedCellPopulation<2> cellPopulation(*p_mesh, cells);

                // Create a simulation
                OffLatticeSimulation<2> simulation(cellPopulation, false, true, true, StepperChoice::EULER);
                std::ostringstream outDir;
                outDir << "VB_AdaptiveEuler_Thresh" << movementThresh;
                setupSimulation(&simulation, outDir.str(), stepsPerHour);

                // Add a cell movement tracker 
                MAKE_PTR_ARGS(SimpleCellCentrePositionTracker<2>, pTracking, (stepsPerHour,1));
                simulation.AddSimulationModifier(pTracking);

                int startT = time(0);
                std::cout << "Dt = " << dt << " Absolute Movement Threshold = " << movementThresh << std::endl;
                simulation.Solve();
                std::cout << "Time elapsed: " << time(0) - startT << endl; 
            }
        } 
    }

    void Test2dVertexBasedWithRK4Stepper() throw (Exception)
    {
        if(toRun == RK4 || toRun == ALL){

            setupTestParameters();

            for(int i = minStepIndex; i <= maxStepIndex; i++){

                double movementThresh = movThresholds[i];
                int stepsPerHour = stepsPerHour_DoesNotExceedMovThresh[i];
                double dt = 1.0/(double)stepsPerHour;

                resetForNewRun();

                // Create a cell population
                HoneycombVertexMeshGenerator generator(meshSide, meshSide);
                MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

                std::vector<CellPtr> cells;
                MAKE_PTR(StemCellProliferativeType, p_stem_type);
                CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
                cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_stem_type);

                VertexBasedCellPopulation<2> cellPopulation(*p_mesh, cells);

                // Create a simulation
                OffLatticeSimulation<2> simulation(cellPopulation, false, true, false, StepperChoice::RK4);
                std::ostringstream outDir;
                outDir << "VB_RK4_Thresh" << movementThresh;
                setupSimulation(&simulation, outDir.str(), stepsPerHour);

                // Add a cell movement tracker 
                MAKE_PTR_ARGS(SimpleCellCentrePositionTracker<2>, pTracking, (stepsPerHour,1));
                simulation.AddSimulationModifier(pTracking);

                int startT = time(0);
                std::cout << "Dt = " << dt << " Absolute Movement Threshold = " << movementThresh << std::endl;
                simulation.Solve();
                std::cout << "Time elapsed: " << time(0) - startT << endl; 
            } 
        }
    }

    void Test2dVertexBasedWithRK4StepperAdaptive() throw (Exception)
    {
        if(toRun == RK4 || toRun == ALL){
        
            setupTestParameters();
        
            for(int i = minStepIndex; i <= maxStepIndex; i++){
        
                double movementThresh = movThresholds[i];
                int stepsPerHour = stepsPerHour_ExceedsMovThresh[i];
                double dt = 1.0/(double)stepsPerHour;
        
                resetForNewRun();
        
                // Create a cell population
                HoneycombVertexMeshGenerator generator(meshSide, meshSide);
                MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        
                std::vector<CellPtr> cells;
                MAKE_PTR(StemCellProliferativeType, p_stem_type);
                CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
                cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_stem_type);
        
                VertexBasedCellPopulation<2> cellPopulation(*p_mesh, cells);
        
                // Create a simulation
                OffLatticeSimulation<2> simulation(cellPopulation, false, true, true, StepperChoice::RK4);
                std::ostringstream outDir;
                outDir << "VB_AdaptiveRK4_Thresh" << movementThresh;
                setupSimulation(&simulation, outDir.str(), stepsPerHour);
        
                // Add a cell movement tracker 
                MAKE_PTR_ARGS(SimpleCellCentrePositionTracker<2>, pTracking, (stepsPerHour,1));
                simulation.AddSimulationModifier(pTracking);
        
                int startT = time(0);
                std::cout << "Dt = " << dt << " Absolute Movement Threshold = " << movementThresh << std::endl;
                simulation.Solve();
                std::cout << "Time elapsed: " << time(0) - startT << endl; 
            } 
        }
    }

    void Test2dVertexBasedWithBackwardEulerStepper() throw (Exception)
    {
        if(toRun == BACKWARDEULER || toRun == ALL){

            setupTestParameters();

            for(int i = minStepIndex; i <= maxStepIndex; i++){

                double movementThresh = movThresholds[i];
                int stepsPerHour = stepsPerHour_DoesNotExceedMovThresh[i];
                double dt = 1.0/(double)stepsPerHour;

                resetForNewRun();

                // Create a cell population
                HoneycombVertexMeshGenerator generator(meshSide, meshSide);
                MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

                std::vector<CellPtr> cells;
                MAKE_PTR(StemCellProliferativeType, p_stem_type);
                CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
                cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_stem_type);

                VertexBasedCellPopulation<2> cellPopulation(*p_mesh, cells);

                // Create a simulation
                OffLatticeSimulation<2> simulation(cellPopulation, false, true, false, StepperChoice::BACKWARDEULER);
                std::ostringstream outDir;
                outDir << "VB_BackwardEuler_Thresh" << movementThresh;
                setupSimulation(&simulation, outDir.str(), stepsPerHour);

                // Add a cell movement tracker 
                MAKE_PTR_ARGS(SimpleCellCentrePositionTracker<2>, pTracking, (stepsPerHour,1));
                simulation.AddSimulationModifier(pTracking);

                int startT = time(0);
                std::cout << "Dt = " << dt << " Absolute Movement Threshold = " << movementThresh << std::endl;
                simulation.Solve();
                std::cout << "Time elapsed: " << time(0) - startT << endl; 
            } 
        }
    }

    void Test2dVertexBasedWithBackwardEulerStepperAdaptive() throw (Exception)
    {
        if(toRun == BACKWARDEULER || toRun == ALL){
        
            setupTestParameters();

            for(int i = minStepIndex; i <= maxStepIndex; i++){

                double movementThresh = movThresholds[i];
                int stepsPerHour = stepsPerHour_ExceedsMovThresh[i];
                double dt = 1.0/(double)stepsPerHour;

                resetForNewRun();

                // Create a cell population
                HoneycombVertexMeshGenerator generator(meshSide, meshSide);
                MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

                std::vector<CellPtr> cells;
                MAKE_PTR(StemCellProliferativeType, p_stem_type);
                CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
                cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_stem_type);

                VertexBasedCellPopulation<2> cellPopulation(*p_mesh, cells);

                // Create a simulation
                OffLatticeSimulation<2> simulation(cellPopulation, false, true, true, StepperChoice::BACKWARDEULER);
                std::ostringstream outDir;
                outDir << "VB_AdaptiveBackwardEuler_Thresh" << movementThresh;
                setupSimulation(&simulation, outDir.str(), stepsPerHour);

                // Add a cell movement tracker 
                MAKE_PTR_ARGS(SimpleCellCentrePositionTracker<2>, pTracking, (stepsPerHour,1));
                simulation.AddSimulationModifier(pTracking);

                int startT = time(0);
                std::cout << "Dt = " << dt << " Absolute Movement Threshold = " << movementThresh << std::endl;
                simulation.Solve();
                std::cout << "Time elapsed: " << time(0) - startT << endl; 
            } 
        }
    }

   
    //-------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------------------------------------------

    // Test settings. Altering ranges here will drastically alter the run time of the test, since
    // a low MinStepIndex will run the test simulation for very small dt. For my work I need to
    // run at a wide range of dt, for production set: minStepIndex = maxStepIndex = 5, endTime = 100.
    int minStepIndex;
    int maxStepIndex;
    int meshSide;
    double endTime;
 
    std::vector<double> movThresholds;                   
    std::vector<int> stepsPerHour_DoesNotExceedMovThresh;
    std::vector<int> stepsPerHour_ExceedsMovThresh;   


    void setupTestParameters(){ 
        
        minStepIndex = 0;
        maxStepIndex = 1;
        meshSide = 3;
        endTime = 100;

        double threshes[12] = {0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
        movThresholds = std::vector<double>(&threshes[0], &threshes[0]+12);

        int steps1[12] = {24000,  12000,     1800, 1300, 590,  295,  147, 73,  50,  38,  30,  25 };
        stepsPerHour_DoesNotExceedMovThresh = std::vector<int>(&steps1[0], &steps1[0]+12);

        int steps2[12] = {10000,  6000,      1800, 1300, 590,  295,  147, 73,  50,  38,  30,  25 };
        stepsPerHour_ExceedsMovThresh = std::vector<int>(&steps2[0], &steps2[0]+12);
    }


    void resetForNewRun(){
        
        RandomNumberGenerator::Instance()->Reseed(0);
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
    }


    void setupSimulation(OffLatticeSimulation<2>* sim, std::string outDir, int stepsPerHour){

        double dt = 1.0/(double)stepsPerHour; 

        sim->SetOutputDirectory(outDir.c_str());
        sim->SetSamplingTimestepMultiple(stepsPerHour);
        sim->SetDt(dt);
        sim->SetEndTime(endTime);

        // Create a force law and pass it to the simulation. 
        // A NagaiHondaForce has to be used together with an AbstractTargetAreaModifier 
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        sim->AddForce(p_nagai_honda_force);   
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        sim->AddSimulationModifier(p_growth_modifier);

        MAKE_PTR_ARGS(SimplePopulationExtentTracker<2>, pExtentTracker, (stepsPerHour));
        sim->AddSimulationModifier(pExtentTracker);
    }
};

#endif /*TESTVERTEXBASEDWITHALTERNATIVETIMESTEPPERS_HPP_*/