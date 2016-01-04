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


#ifndef TESTFULLBUSKEMODEL_HPP
#define TESTFULLBUSKEMODEL_HPP

#include <cxxtest/cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
// Must be included before other cell_based headers

#include "CellBasedSimulationArchiver.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "NodeBasedCellPopulationWithBuskeUpdate.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "BuskeAdhesiveForce.hpp"
#include "BuskeElasticForce.hpp"
#include "BuskeCompressionForce.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscTools.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "CellTrackingOutput.hpp"

class TestFullBuskeModel : public AbstractCellBasedTestSuite {

private:

    NodeBasedCellPopulationWithBuskeUpdate<3>* Create3DBuskeCellPair(c_vector<double, 3> pos1, c_vector<double, 3> pos2,
                                                                     double relaxedRadius1, double relaxedRadius2,
                                                                     double startingRadius1, double startingRadius2,
                                                                     boost::shared_ptr<AbstractCellProliferativeType> prolifType){
        //CREATE MESH
        std::vector< Node<3>* > nodes;
        Node<3>* node1 = new Node<3>(0, false, pos1[0], pos1[1], pos1[2]);
        node1->SetRadius(startingRadius1);
        nodes.push_back(node1);
        Node<3>* node2 = new Node<3>(0, false, pos2[0], pos2[1], pos2[2]);
        node2->SetRadius(startingRadius2);
        nodes.push_back(node2);

        NodesOnlyMesh<3>* mesh = new NodesOnlyMesh<3>;
        mesh->ConstructNodesWithoutMesh(nodes, 50);
        delete nodes[0];
        delete nodes[1];

        //CREATE CELLS
        std::vector<CellPtr> cells;
        cells.clear();
        cells.reserve(2);
        for (int i=0; i < 2; i++)
        {
            StochasticDurationCellCycleModel* p_cell_cycle_model = new StochasticDurationCellCycleModel();
            p_cell_cycle_model->SetDimension(3);

            boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
            CellPtr new_cell(new Cell(p_state, p_cell_cycle_model));

            new_cell->SetCellProliferativeType(prolifType);

            double birth_time = -p_cell_cycle_model->GetAverageStemCellCycleTime()*RandomNumberGenerator::Instance()->ranf();
            new_cell->SetBirthTime(birth_time);

            if(i==0){
                new_cell->GetCellData()->SetItem("RelaxedRadius", relaxedRadius1);
                new_cell->GetCellData()->SetItem("Radius", startingRadius1);
                new_cell->GetCellData()->SetItem("IsBuskeKnot", 0);
            }else{
                new_cell->GetCellData()->SetItem("RelaxedRadius", relaxedRadius2);
                new_cell->GetCellData()->SetItem("Radius", startingRadius2);
                new_cell->GetCellData()->SetItem("IsBuskeKnot", 0);
            }
            cells.push_back(new_cell);
        }

        NodeBasedCellPopulationWithBuskeUpdate<3>* population = new NodeBasedCellPopulationWithBuskeUpdate<3>(*mesh, cells);
        return( population );
    }


public:

    void TestCellPairWithAdhesionOnlyNoRadialChanges() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        double Eps = 2000*pow(10,-6);

        c_vector<double, 3> pos1;
        pos1[0]=0; pos1[1]=0; pos1[2]=0;
        c_vector<double, 3> pos2;
        pos2[0]=0; pos2[1]=9*pow(10,-6); pos2[2]=0;

        double relaxedRadius1 = 5*pow(10,-6);
        double relaxedRadius2 = 10*pow(10,-6);
        double startingRadius1 = 5*pow(10,-6);
        double startingRadius2 = 10*pow(10,-6);

        MAKE_PTR(DifferentiatedCellProliferativeType, differentiated);

        NodeBasedCellPopulationWithBuskeUpdate<3>* cellPopulation =
                Create3DBuskeCellPair(pos1, pos2, relaxedRadius1, relaxedRadius2, startingRadius1, startingRadius2, differentiated);

        cellPopulation->SetUseVaryingRadii(false);
        cellPopulation->SetUseVariableRadii(true);
        cellPopulation->SetDampingConstantIntercell(5*pow(10,10));
        cellPopulation->SetDampingConstantMedium(3.2);

        //SETUP SIMULATION
        OffLatticeSimulation<3> simulator(*cellPopulation);
        simulator.SetOutputDirectory("TestBuskeCellPairWithAdhesionNoRadial");
        simulator.SetEndTime(600.0);
        simulator.SetDt(1.0/250.0);
        simulator.SetSamplingTimestepMultiple(250);

        //ADD FORCES
        cellPopulation->EnableAdhesion(Eps);
        MAKE_PTR(BuskeAdhesiveForce<3>, p_buske_adhesive_force);
        p_buske_adhesive_force->SetAdhesionEnergyParameter(Eps);
        simulator.AddForce(p_buske_adhesive_force);

        //RECORD CELL TRACKS
        MAKE_PTR_ARGS(CellTrackingOutput<3>, tracking, (250, 1));
        simulator.AddSimulationModifier(tracking);

        simulator.Solve();
    }


    void TestCellPairWithAdhesionOnlyRadialChangesOnly() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        double Eps = 4000*pow(10,-6);

        c_vector<double, 3> pos1;
        pos1[0]=0; pos1[1]=0; pos1[2]=0;
        c_vector<double, 3> pos2;
        pos2[0]=0; pos2[1]=10*pow(10,-6); pos2[2]=0;

        double relaxedRadius1 = 5*pow(10,-6);
        double relaxedRadius2 = 10*pow(10,-6);
        double startingRadius1 = 5*pow(10,-6);
        double startingRadius2 = 10*pow(10,-6);

        MAKE_PTR(DifferentiatedCellProliferativeType, differentiated);

        NodeBasedCellPopulationWithBuskeUpdate<3>* cellPopulation =
                Create3DBuskeCellPair(pos1, pos2, relaxedRadius1, relaxedRadius2, startingRadius1, startingRadius2, differentiated);

        cellPopulation->SetUseVaryingRadii(true);
        cellPopulation->SetUseVariableRadii(true);
        cellPopulation->SetDisableMovement(true);
        cellPopulation->SetDampingConstantIntercell(5*pow(10,10));
        cellPopulation->SetDampingConstantVolume(4);
        cellPopulation->EnableAdhesion(Eps);

        //SETUP SIMULATION
        OffLatticeSimulation<3> simulator(*cellPopulation);
        simulator.SetOutputDirectory("TestBuskeCellPairWithAdhesionRadialChangesOnly");
        simulator.SetEndTime(4000.0);
        simulator.SetDt(1.0/250.0);
        simulator.SetSamplingTimestepMultiple(250);

        //RECORD CELL TRACKS
        MAKE_PTR_ARGS(CellTrackingOutput<3>, tracking, (250, 1));
        simulator.AddSimulationModifier(tracking);

        simulator.Solve();
    }


    void TestCellPairWithElasticityOnlyNoRadialChanges() throw (Exception) {
        EXIT_IF_PARALLEL;

        double D = 0.1*(4.0 / 3.0) * pow(10, -3);

        c_vector<double, 3> pos1;
        pos1[0] = 0;
        pos1[1] = 0;
        pos1[2] = 0;
        c_vector<double, 3> pos2;
        pos2[0] = 0;
        pos2[1] = 10 * pow(10, -6);
        pos2[2] = 0;

        double relaxedRadius1 = 5 * pow(10, -6);
        double relaxedRadius2 = 10 * pow(10, -6);
        double startingRadius1 = 5 * pow(10, -6);
        double startingRadius2 = 10 * pow(10, -6);

        MAKE_PTR(DifferentiatedCellProliferativeType, differentiated);

        NodeBasedCellPopulationWithBuskeUpdate<3> *cellPopulation =
                Create3DBuskeCellPair(pos1, pos2, relaxedRadius1, relaxedRadius2, startingRadius1, startingRadius2,
                                      differentiated);

        cellPopulation->SetUseVaryingRadii(false);
        cellPopulation->SetUseVariableRadii(true);
        cellPopulation->SetDampingConstantIntercell(5 * pow(10, 10));
        cellPopulation->SetDampingConstantMedium(3.2);

        //SETUP SIMULATION
        OffLatticeSimulation<3> simulator(*cellPopulation);
        simulator.SetOutputDirectory("TestBuskeCellPairWithElasticityNoRadial");
        simulator.SetEndTime(1200.0);
        simulator.SetDt(1.0 / 250.0);
        simulator.SetSamplingTimestepMultiple(250);

        //ADD FORCES
        cellPopulation->EnableElasticity(D);
        MAKE_PTR(BuskeElasticForce<3>, p_buske_elast_force);
        p_buske_elast_force->SetDeformationEnergyParameter(D);
        simulator.AddForce(p_buske_elast_force);

        //RECORD CELL TRACKS
        MAKE_PTR_ARGS(CellTrackingOutput<3>, tracking, (250, 1));
        simulator.AddSimulationModifier(tracking);

        simulator.Solve();
    }



    void TestCellPairWithElasticityOnlyRadialChangesOnly() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        double D = 0.1*(4.0 / 3.0) * pow(10, -3);

        c_vector<double, 3> pos1;
        pos1[0]=0; pos1[1]=0; pos1[2]=0;
        c_vector<double, 3> pos2;
        pos2[0]=0; pos2[1]=10*pow(10,-6); pos2[2]=0;

        double relaxedRadius1 = 5*pow(10,-6);
        double relaxedRadius2 = 10*pow(10,-6);
        double startingRadius1 = 5*pow(10,-6);
        double startingRadius2 = 10*pow(10,-6);

        MAKE_PTR(DifferentiatedCellProliferativeType, differentiated);

        NodeBasedCellPopulationWithBuskeUpdate<3>* cellPopulation =
                Create3DBuskeCellPair(pos1, pos2, relaxedRadius1, relaxedRadius2, startingRadius1, startingRadius2, differentiated);

        cellPopulation->SetUseVaryingRadii(true);
        cellPopulation->SetDisableMovement(true);
        cellPopulation->SetUseVariableRadii(true);
        cellPopulation->SetDampingConstantIntercell(5*pow(10,10));
        cellPopulation->SetDampingConstantVolume(4);
        cellPopulation->EnableElasticity(D);

        //SETUP SIMULATION
        OffLatticeSimulation<3> simulator(*cellPopulation);
        simulator.SetOutputDirectory("TestBuskeCellPairWithElasticityRadialChangesOnly");
        simulator.SetEndTime(1200.0);
        simulator.SetDt(1.0/250.0);
        simulator.SetSamplingTimestepMultiple(1);

        //RECORD CELL TRACKS
        MAKE_PTR_ARGS(CellTrackingOutput<3>, tracking, (250, 1));
        simulator.AddSimulationModifier(tracking);

        simulator.Solve();
    }



    void TestCellPairWithCompressionOnlyNoRadialChanges() throw (Exception) {
        EXIT_IF_PARALLEL;

        double K = 100000;

        c_vector<double, 3> pos1;
        pos1[0] = 0;
        pos1[1] = 0;
        pos1[2] = 0;
        c_vector<double, 3> pos2;
        pos2[0] = 0;
        pos2[1] = 10 * pow(10, -6);
        pos2[2] = 0;

        double relaxedRadius1 = 5 * pow(10, -6);
        double relaxedRadius2 = 10 * pow(10, -6);
        double startingRadius1 = 5 * pow(10, -6);
        double startingRadius2 = 10 * pow(10, -6);

        MAKE_PTR(DifferentiatedCellProliferativeType, differentiated);

        NodeBasedCellPopulationWithBuskeUpdate<3> *cellPopulation =
                Create3DBuskeCellPair(pos1, pos2, relaxedRadius1, relaxedRadius2, startingRadius1, startingRadius2,
                                      differentiated);

        cellPopulation->SetUseVaryingRadii(false);
        cellPopulation->SetUseVariableRadii(true);
        cellPopulation->SetDampingConstantIntercell(5 * pow(10, 10));
        cellPopulation->SetDampingConstantMedium(3.2);

        //SETUP SIMULATION
        OffLatticeSimulation<3> simulator(*cellPopulation);
        simulator.SetOutputDirectory("TestBuskeCellPairWithCompressionNoRadial");
        simulator.SetEndTime(4000.0);
        simulator.SetDt(1.0 / 250.0);
        simulator.SetSamplingTimestepMultiple(250);

        //ADD FORCES
        cellPopulation->EnableCompression(K);
        MAKE_PTR(BuskeCompressionForce<3>, p_buske_comp_force);
        p_buske_comp_force->SetCompressionEnergyParameter(K);
        simulator.AddForce(p_buske_comp_force);

        //RECORD CELL TRACKS
        MAKE_PTR_ARGS(CellTrackingOutput<3>, tracking, (250, 1));
        simulator.AddSimulationModifier(tracking);

        simulator.Solve();
    }



    void TestCellPairWithCompressionOnlyRadialChangesOnly() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        double K = 100000;

        c_vector<double, 3> pos1;
        pos1[0]=0; pos1[1]=0; pos1[2]=0;
        c_vector<double, 3> pos2;
        pos2[0]=0; pos2[1]=10*pow(10,-6); pos2[2]=0;

        double relaxedRadius1 = 5*pow(10,-6);
        double relaxedRadius2 = 10*pow(10,-6);
        double startingRadius1 = 5*pow(10,-6);
        double startingRadius2 = 10*pow(10,-6);

        MAKE_PTR(DifferentiatedCellProliferativeType, differentiated);

        NodeBasedCellPopulationWithBuskeUpdate<3>* cellPopulation =
                Create3DBuskeCellPair(pos1, pos2, relaxedRadius1, relaxedRadius2, startingRadius1, startingRadius2, differentiated);

        cellPopulation->SetUseVaryingRadii(true);
        cellPopulation->SetUseVariableRadii(true);
        cellPopulation->SetDisableMovement(true);
        cellPopulation->SetDampingConstantIntercell(5*pow(10,10));
        cellPopulation->SetDampingConstantVolume(400);
        cellPopulation->EnableCompression(K);

        //SETUP SIMULATION
        OffLatticeSimulation<3> simulator(*cellPopulation);
        simulator.SetOutputDirectory("TestBuskeCellPairWithCompressionRadialChangesOnly");
        simulator.SetEndTime(600.0);
        simulator.SetDt(1.0/250.0);
        simulator.SetSamplingTimestepMultiple(250);

        //RECORD CELL TRACKS
        MAKE_PTR_ARGS(CellTrackingOutput<3>, tracking, (250, 1));
        simulator.AddSimulationModifier(tracking);

        simulator.Solve();
    }

    void TestIsolatedCellsWithCompressionOnlyRadialChangesOnly() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        double K = 100000;

        c_vector<double, 3> pos1;
        pos1[0]=0; pos1[1]=0; pos1[2]=0;
        c_vector<double, 3> pos2;
        pos2[0]=0; pos2[1]=20*pow(10,-6); pos2[2]=0;

        double relaxedRadius1 = 5*pow(10,-6);
        double relaxedRadius2 = 10*pow(10,-6);
        double startingRadius1 = 10*pow(10,-6);
        double startingRadius2 = 5*pow(10,-6);

        MAKE_PTR(DifferentiatedCellProliferativeType, differentiated);

        NodeBasedCellPopulationWithBuskeUpdate<3>* cellPopulation =
                Create3DBuskeCellPair(pos1, pos2, relaxedRadius1, relaxedRadius2, startingRadius1, startingRadius2, differentiated);

        cellPopulation->SetUseVaryingRadii(true);
        cellPopulation->SetUseVariableRadii(true);
        cellPopulation->SetDisableMovement(true);
        cellPopulation->SetDampingConstantIntercell(5*pow(10,10));
        cellPopulation->SetDampingConstantVolume(400);
        cellPopulation->EnableCompression(K);

        //SETUP SIMULATION
        OffLatticeSimulation<3> simulator(*cellPopulation);
        simulator.SetOutputDirectory("TestBuskeCellsInIsolationWithCompressionRadialChangesOnly");
        simulator.SetEndTime(600.0);
        simulator.SetDt(1.0/250.0);
        simulator.SetSamplingTimestepMultiple(250);

        //RECORD CELL TRACKS
        MAKE_PTR_ARGS(CellTrackingOutput<3>, tracking, (1, 1));
        simulator.AddSimulationModifier(tracking);

        simulator.Solve();
    }

};

#endif //TESTFULLBUSKEMODEL_HPP
