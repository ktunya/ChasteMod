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

#ifndef TESTNODEBASEDWITHDOP853TIMESTEPPER_HPP_
#define TESTNODEBASEDWITHDOP853TIMESTEPPER_HPP_

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
#include "FixedDivisionTimingsCellCycleModel.hpp"
#include "SimpleCellCentrePositionTracker.hpp"
#include "SimplePopulationExtentTracker.hpp"

class TestNodeBasedWithDOP853Timestepper : public AbstractCellBasedWithTimingsTestSuite
{

public:

    enum CycleModel {STOCHASTIC, FIXEDTIMINGS};
    static const int CCmodel = FIXEDTIMINGS;

    // A 3D node based test simulation, run using an adaptive, embedded RK4 scheme
    // (specifically Dormand-Prince 853) 

    void Test3dNodeBasedWithDOP853() throw (Exception)
    {

        double movementThresh = 10;
        int stepsPerHour = 50;
        double dt = 1.0/(double)stepsPerHour;

        resetForNewRun();

        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, false,  0.5, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, false, -0.5, 0.0, 0.0));
        
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 3.0);

        std::vector<CellPtr> cells;
        GenerateStemCells(cells, mesh.GetNumNodes());

        NodeBasedCellPopulation<3>* cellPopulation = new NodeBasedCellPopulation<3>(mesh, cells);
        cellPopulation->SetAbsoluteMovementThreshold(movementThresh);
        cellPopulation->SetDampingConstantNormal(1.1);

        OffLatticeSimulation<3> simulation((*cellPopulation), false, true, true, StepperChoice::DOP853);

        std::ostringstream outDir;
        outDir << "NB_DOP853_Adaptive" << stepsPerHour;

        setupSimulation(&simulation, outDir.str(), stepsPerHour);

        // Add a cell movement tracker 
        MAKE_PTR_ARGS(SimpleCellCentrePositionTracker<3>, pTracking, (stepsPerHour,1));
        simulation.AddSimulationModifier(pTracking);

        int startT = time(0);
        std::cout << "Dt = " << dt << " Absolute Movement Threshold = " << movementThresh << std::endl;
        simulation.Solve();
        std::cout << "Time elapsed: " << time(0) - startT << endl; 

        delete nodes[0];
        delete nodes[1];  
    }


    //----------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------

    void resetForNewRun(){
        
        RandomNumberGenerator::Instance()->Reseed(0);
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        CellId::ResetMaxCellId();
    }


    void setupSimulation(OffLatticeSimulation<3>* sim, std::string outDir, int stepsPerHour){

        double dt = 1.0/(double)stepsPerHour; 

        sim->SetOutputDirectory(outDir.c_str());
        sim->SetSamplingTimestepMultiple(stepsPerHour);
        sim->SetDt(dt);
        sim->SetEndTime(500);

        // Create a force law
        MAKE_PTR(GeneralisedLinearSpringForce<3>, pLinearForce);
        pLinearForce->SetCutOffLength(3.0);
        sim->AddForce(pLinearForce);

        MAKE_PTR_ARGS(SimplePopulationExtentTracker<3>, pExtentTracker, (stepsPerHour));
        sim->AddSimulationModifier(pExtentTracker);
    }


    void GenerateStemCells(std::vector<CellPtr>& rCells, unsigned numCells)
    {
   
        rCells.clear();
        rCells.reserve(numCells);

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
            p_cell_cycle_model->SetDimension(3);

            boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));
            p_cell->SetCellProliferativeType(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

            double birth_time = -p_cell_cycle_model->GetAverageStemCellCycleTime()*RandomNumberGenerator::Instance()->ranf();
            p_cell->SetBirthTime(birth_time);
            rCells.push_back(p_cell);
        }
    }
};

#endif /*TESTNODEBASEDWITHDOP853TIMESTEPPER_HPP_*/
