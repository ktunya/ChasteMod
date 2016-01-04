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

#ifndef TESTMESHBASEDWITHALTERNATIVETIMESTEPPERS_HPP_
#define TESTMESHBASEDWITHALTERNATIVETIMESTEPPERS_HPP_

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
#include "HoneycombMeshGenerator.hpp"

#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "VoronoiDataWriter.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "FixedDivisionTimingsCellCycleModel.hpp"
#include "SimpleCellCentrePositionTracker.hpp"
#include "SimplePopulationExtentTracker.hpp"
#include "StemCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"

class TestMeshBasedWithAlternativeTimesteppers : public AbstractCellBasedWithTimingsTestSuite
{
public:

    enum RunChoices {EULER, RK4, BACKWARDEULER, ADAMSMOULTON, ALL};
    static const int toRun = ADAMSMOULTON;

    enum CycleModel {STOCHASTIC, FIXEDTIMINGS};
    static const int CCmodel = STOCHASTIC;


    void TestMeshBasedWithEulerStepper() throw (Exception)
    {
    	if(toRun == EULER || toRun == ALL){

            setupTestParameters();

            for(int i = minStepIndex; i <= maxStepIndex; i++){

                double movementThresh = movThresholds[i];
                int stepsPerHour = stepsPerHour_DoesNotExceedMovThresh[i];
                double dt = 1.0/(double)stepsPerHour;
                resetForNewRun();

                // Create a 2d mesh-based cell population with 2 rows of ghost nodes
				HoneycombMeshGenerator generator(meshSide, meshSide, 2);
        		MutableMesh<2,2>* p_mesh = generator.GetMesh();
        		std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        		std::vector<CellPtr> cells;
        		GenerateStemCells(cells, location_indices.size(), location_indices);
		
        		MeshBasedCellPopulationWithGhostNodes<2> cellPopulation(*p_mesh, cells, location_indices);
        		cellPopulation.CreateVoronoiTessellation();
        		cellPopulation.AddPopulationWriter<VoronoiDataWriter>();
        		cellPopulation.SetAbsoluteMovementThreshold(movementThresh);
                cellPopulation.SetDampingConstantNormal(1.1);
                
                // Create a simulation-----------------------------------------------------------------------
                OffLatticeSimulation<2> simulation(cellPopulation, false, true, false, StepperChoice::EULER);
                std::ostringstream outDir;
                outDir << "MBG_Euler_Fixed_Rep" << i;
                setupSimulation(&simulation, outDir.str(), stepsPerHour);
                // ------------------------------------------------------------------------------------------

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


    void TestMeshBasedWithRK4Stepper() throw (Exception)
    {
    	if(toRun == RK4 || toRun == ALL){

            setupTestParameters();

            for(int i = minStepIndex; i <= maxStepIndex; i++){

                double movementThresh = movThresholds[i];
                int stepsPerHour = stepsPerHour_DoesNotExceedMovThresh[i];
                double dt = 1.0/(double)stepsPerHour;
                resetForNewRun();

                // Create a 2d mesh-based cell population with 2 rows of ghost nodes
                HoneycombMeshGenerator generator(meshSide, meshSide, 2);
                MutableMesh<2,2>* p_mesh = generator.GetMesh();
                std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

                std::vector<CellPtr> cells;
                GenerateStemCells(cells, location_indices.size(), location_indices);
        
                MeshBasedCellPopulationWithGhostNodes<2> cellPopulation(*p_mesh, cells, location_indices);
                cellPopulation.CreateVoronoiTessellation();
                cellPopulation.AddPopulationWriter<VoronoiDataWriter>();
                cellPopulation.SetAbsoluteMovementThreshold(movementThresh);
                cellPopulation.SetDampingConstantNormal(1.1);
                
                // Create a simulation-----------------------------------------------------------------------
                OffLatticeSimulation<2> simulation(cellPopulation, false, true, false, StepperChoice::RK4);
                std::ostringstream outDir;
                outDir << "MBG_RK4_Fixed_Rep" << i;
                setupSimulation(&simulation, outDir.str(), stepsPerHour);
                // ------------------------------------------------------------------------------------------

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

    void TestMeshBasedWithBackwardEulerStepper() throw (Exception)
    {
    	if(toRun == BACKWARDEULER || toRun == ALL){

            setupTestParameters();

            for(int i = minStepIndex; i <= maxStepIndex; i++){

                double movementThresh = movThresholds[i];
                int stepsPerHour = stepsPerHour_DoesNotExceedMovThresh[i];
                double dt = 1.0/(double)stepsPerHour;
                resetForNewRun();

                // Create a 2d mesh-based cell population with 2 rows of ghost nodes
				HoneycombMeshGenerator generator(meshSide, meshSide, 2);
        		MutableMesh<2,2>* p_mesh = generator.GetMesh();
        		std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        		std::vector<CellPtr> cells;
        		GenerateStemCells(cells, location_indices.size(), location_indices);
		
        		MeshBasedCellPopulationWithGhostNodes<2> cellPopulation(*p_mesh, cells, location_indices);
        		cellPopulation.CreateVoronoiTessellation();
        		cellPopulation.AddPopulationWriter<VoronoiDataWriter>();
        		cellPopulation.SetAbsoluteMovementThreshold(movementThresh);
                cellPopulation.SetDampingConstantNormal(1.1);
                
                // Create a simulation-----------------------------------------------------------------------
                OffLatticeSimulation<2> simulation(cellPopulation, false, true, false, StepperChoice::BACKWARDEULER);
                std::ostringstream outDir;
                outDir << "MBG_BackwardEuler_Rep" << i;
                setupSimulation(&simulation, outDir.str(), stepsPerHour);
                // ------------------------------------------------------------------------------------------

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


    void TestMeshBasedWithAdamsMoultonStepper() throw (Exception)
    {
        if(toRun == ADAMSMOULTON || toRun == ALL){

            setupTestParameters();

            for(int i = minStepIndex; i <= maxStepIndex; i++){

                double movementThresh = movThresholds[i];
                int stepsPerHour = stepsPerHour_DoesNotExceedMovThresh[i];
                double dt = 1.0/(double)stepsPerHour;
                resetForNewRun();

                // Create a 2d mesh-based cell population with 2 rows of ghost nodes
                HoneycombMeshGenerator generator(meshSide, meshSide, 2);
                MutableMesh<2,2>* p_mesh = generator.GetMesh();
                std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

                std::vector<CellPtr> cells;
                GenerateStemCells(cells, location_indices.size(), location_indices);
        
                MeshBasedCellPopulationWithGhostNodes<2> cellPopulation(*p_mesh, cells, location_indices);
                cellPopulation.CreateVoronoiTessellation();
                cellPopulation.AddPopulationWriter<VoronoiDataWriter>();
                cellPopulation.SetAbsoluteMovementThreshold(movementThresh);
                cellPopulation.SetDampingConstantNormal(1.1);
                
                // Create a simulation-----------------------------------------------------------------------
                OffLatticeSimulation<2> simulation(cellPopulation, false, true, false, StepperChoice::ADAMSMOULTON);
                std::ostringstream outDir;
                outDir << "MBG_AdamsMoulton_Rep" << i;
                setupSimulation(&simulation, outDir.str(), stepsPerHour);
                // ------------------------------------------------------------------------------------------

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
        
        minStepIndex = 3;
        maxStepIndex = 8;
        meshSide = 3;
        endTime = 100;
 
        double threshes[9] = {10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0};
        movThresholds = std::vector<double>(&threshes[0], &threshes[0]+9);

        stepsPerHour_DoesNotExceedMovThresh = std::vector<int>();
        stepsPerHour_ExceedsMovThresh       = std::vector<int>();
        for(int i=14; i>5; i--){
            stepsPerHour_DoesNotExceedMovThresh.push_back((int)pow(2,i));
            stepsPerHour_ExceedsMovThresh.push_back( (int)pow(2,(i-1)) );
        } 
    }


    void resetForNewRun(){
        
        RandomNumberGenerator::Instance()->Reseed(0);
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        CellId::ResetMaxCellId();
    }


    void setupSimulation(OffLatticeSimulation<2>* sim, std::string outDir, int stepsPerHour){

        double dt = 1.0/(double)stepsPerHour; 

        sim->SetOutputDirectory(outDir.c_str());
        sim->SetSamplingTimestepMultiple(stepsPerHour);
        sim->SetDt(dt);
        sim->SetEndTime(endTime);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, pLinearForce);
        pLinearForce->SetCutOffLength(3);
        sim->AddForce(pLinearForce);

        MAKE_PTR_ARGS(SimplePopulationExtentTracker<2>, pExtentTracker, (stepsPerHour));
        sim->AddSimulationModifier(pExtentTracker);
    }


    void GenerateStemCells(std::vector<CellPtr>& rCells, unsigned numCells, const std::vector<unsigned> locationIndices)
    {
   
        rCells.clear();
        // If location indices is given, then it needs to match the number of output cells
        if (numCells != locationIndices.size())
        {
            EXCEPTION("The size of the locationIndices vector must match the required number of output cells");
        }
        rCells.reserve(numCells);

        // Create cells
        for (unsigned i=0; i<numCells; i++)
        {   
            AbstractCellCycleModel* p_cell_cycle_model;
            if(CCmodel == STOCHASTIC){
                p_cell_cycle_model = new StochasticDurationCellCycleModel();
            }else if(CCmodel == FIXEDTIMINGS){
                p_cell_cycle_model = new FixedDivisionTimingsCellCycleModel(24);
            }else{
                EXCEPTION("Choose a cell cycle model type");
            }
            p_cell_cycle_model->SetDimension(2);

            boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));
            p_cell->SetCellProliferativeType(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

            double birth_time = 0.0 - locationIndices[i];
            p_cell->SetBirthTime(birth_time);
            rCells.push_back(p_cell);
        }
    }

};

#endif /*TESTMESHBASEDWITHALTERNATIVETIMESTEPPERS_HPP_*/ 