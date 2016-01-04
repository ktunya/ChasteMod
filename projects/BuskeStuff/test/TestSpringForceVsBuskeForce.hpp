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


#ifndef TESTSPRINGFORCEVSBUSKEFORCE_HPP
#define TESTSPRINGFORCEVSBUSKEFORCE_HPP

#include <cxxtest/cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"
// Must be included before other cell_based headers

#include "CellBasedSimulationArchiver.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "NodeBasedCellPopulationWithBuskeUpdate.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "BuskeAdhesiveForce.hpp"
#include "BuskeElasticForce.hpp"
#include "BuskeCompressionForce.hpp"
#include "RepulsionForce.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscTools.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "CellTrackingOutput.hpp"

class TestSpringForceVsBuskeForce : public AbstractCellBasedTestSuite {

private:

    NodesOnlyMesh<3>* CreateMeshBuske(){

        double startingRadius1 = 3*pow(10,-6);
        double startingRadius2 = 3*pow(10,-6);

        c_vector<double, 3> pos1;
        pos1[0]=0; pos1[1]=0; pos1[2]=0;
        c_vector<double, 3> pos2;
        pos2[0]=0; pos2[1]=1*pow(10,-6); pos2[2]=0;

        NodesOnlyMesh<3>* mesh = CreateMesh(pos1, pos2, startingRadius1, startingRadius2);
        return(mesh);     
    }

    NodesOnlyMesh<3>* CreateMeshSpring(){

        double relaxedRadius1 = 3;
        double relaxedRadius2 = 3;
        double startingRadius1 = 3;
        double startingRadius2 = 3;

        c_vector<double, 3> pos1;
        pos1[0]=0; pos1[1]=0; pos1[2]=0;
        c_vector<double, 3> pos2;
        pos2[0]=0; pos2[1]=1; pos2[2]=0;

        NodesOnlyMesh<3>* mesh = CreateMesh(pos1, pos2, startingRadius1, startingRadius2);
        return(mesh);      
    }

    NodesOnlyMesh<3>* CreateMesh(c_vector<double,3> pos1, c_vector<double,3>  pos2, double startingRadius1, double startingRadius2){
        
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
        return mesh; 
    }

    std::vector<CellPtr> CreateCellsBuske(){

        double relaxedRadius1 = 3*pow(10,-6);
        double relaxedRadius2 = 3*pow(10,-6);
        double startingRadius1 = 3*pow(10,-6);
        double startingRadius2 = 3*pow(10,-6);

        std::vector<CellPtr> Cells = CreateCells(relaxedRadius1, relaxedRadius2, startingRadius1, startingRadius2, false);
        return(Cells);
    }

    std::vector<CellPtr> CreateCellsSpring(){

        double relaxedRadius1 = 3;
        double relaxedRadius2 = 3;
        double startingRadius1 = 3;
        double startingRadius2 = 3;

        std::vector<CellPtr> Cells = CreateCells(relaxedRadius1, relaxedRadius2, startingRadius1, startingRadius2, true);
        return(Cells);
    }

    std::vector<CellPtr> CreateCells(double relaxedRadius1, double relaxedRadius2,
                                     double startingRadius1, double startingRadius2, bool setRadius){

        MAKE_PTR(DifferentiatedCellProliferativeType, differentiated);

        std::vector<CellPtr> cells;
        cells.clear();
        cells.reserve(2);
        for (int i=0; i < 2; i++)
        {
            StochasticDurationCellCycleModel* p_cell_cycle_model = new StochasticDurationCellCycleModel();
            p_cell_cycle_model->SetDimension(3);

            boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
            CellPtr new_cell(new Cell(p_state, p_cell_cycle_model));

            new_cell->SetCellProliferativeType(differentiated);

            double birth_time = -p_cell_cycle_model->GetAverageStemCellCycleTime()*RandomNumberGenerator::Instance()->ranf();
            new_cell->SetBirthTime(birth_time);

            if(i==0){
                new_cell->GetCellData()->SetItem("RelaxedRadius", relaxedRadius1);
                new_cell->GetCellData()->SetItem("Radius", startingRadius1);
            }else{
                new_cell->GetCellData()->SetItem("RelaxedRadius", relaxedRadius2);
                new_cell->GetCellData()->SetItem("Radius", startingRadius2);
            }
            cells.push_back(new_cell);
        }  
        return( cells );
    };

public:

    void TestSpring() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        std::vector<CellPtr> springcells = CreateCellsSpring();
        NodesOnlyMesh<3>* springmesh = CreateMeshSpring();
        NodeBasedCellPopulation<3>* springCellPopulation = 
            new NodeBasedCellPopulation<3>(*springmesh, springcells);
        springCellPopulation->SetDampingConstantNormal(1);
        springCellPopulation->SetUseVariableRadii(true);

        //SETUP SIMULATION 2
        OffLatticeSimulation<3> simulator(*springCellPopulation);
        simulator.SetOutputDirectory("Spring");
        simulator.SetEndTime(600.0);
        simulator.SetDt(1.0/250.0);
        simulator.SetSamplingTimestepMultiple(250);

        //ADD FORCES
        MAKE_PTR(RepulsionForce<3>, p_spring_force);
        p_spring_force->SetMeinekeSpringStiffness(0.00106);
        simulator.AddForce(p_spring_force);
        //RECORD CELL TRACKS
        MAKE_PTR_ARGS(CellTrackingOutput<3>, tracking2, (250, 1));
        simulator.AddSimulationModifier(tracking2);
        //simulator.Solve();
    }


    void TestBuskeVsSpringNoRadial() throw (Exception)
    {
        EXIT_IF_PARALLEL;
        double D = (4/3)*pow(10,-3);

        std::vector<CellPtr> buskecells = CreateCellsBuske();
        NodesOnlyMesh<3>* buskemesh = CreateMeshBuske();
        NodeBasedCellPopulationWithBuskeUpdate<3>* buskeCellPopulation = 
            new NodeBasedCellPopulationWithBuskeUpdate<3>(*buskemesh, buskecells);
                
        buskeCellPopulation->SetUseVaryingRadii(false);
        buskeCellPopulation->SetUseVariableRadii(true);
        buskeCellPopulation->SetDampingConstantIntercell(1*pow(10,-10));
        buskeCellPopulation->SetDampingConstantMedium(1);

        //SETUP SIMULATION 1
        OffLatticeSimulation<3> simulator(*buskeCellPopulation);
        simulator.SetOutputDirectory("TestBuskeVsSpringNoRadial");
        simulator.SetEndTime(600.0);
        simulator.SetDt(1.0/250.0);
        simulator.SetSamplingTimestepMultiple(250);
        //ADD FORCES
        buskeCellPopulation->EnableElasticity(D);
        MAKE_PTR(BuskeElasticForce<3>, p_buske_elastic_force);
        p_buske_elastic_force->SetDeformationEnergyParameter(D);
        simulator.AddForce(p_buske_elastic_force);
        //RECORD CELL TRACKS
        MAKE_PTR_ARGS(CellTrackingOutput<3>, tracking, (250, 1));
        simulator.AddSimulationModifier(tracking);

        //simulator.Solve();        
    }


    void TestBuskeVsSpringNoRadialwCellCellFriction() throw (Exception)
    {
        EXIT_IF_PARALLEL;
        double D = (4/3)*pow(10,-3);

        std::vector<CellPtr> buskecells = CreateCellsBuske();
        NodesOnlyMesh<3>* buskemesh = CreateMeshBuske();
        NodeBasedCellPopulationWithBuskeUpdate<3>* buskeCellPopulation = 
            new NodeBasedCellPopulationWithBuskeUpdate<3>(*buskemesh, buskecells);
                
        buskeCellPopulation->SetUseVaryingRadii(false);
        buskeCellPopulation->SetUseVariableRadii(true);
        buskeCellPopulation->SetDampingConstantIntercell(5*pow(10,10));
        buskeCellPopulation->SetDampingConstantMedium(1);

        //SETUP SIMULATION 1
        OffLatticeSimulation<3> simulator(*buskeCellPopulation);
        simulator.SetOutputDirectory("TestBuskeVsSpringNoRadialwCellCellFriction");
        simulator.SetEndTime(600.0);
        simulator.SetDt(1.0/250.0);
        simulator.SetSamplingTimestepMultiple(250);
        //ADD FORCES
        buskeCellPopulation->EnableElasticity(D);
        MAKE_PTR(BuskeElasticForce<3>, p_buske_elastic_force);
        p_buske_elastic_force->SetDeformationEnergyParameter(D);
        simulator.AddForce(p_buske_elastic_force);
        //RECORD CELL TRACKS
        MAKE_PTR_ARGS(CellTrackingOutput<3>, tracking, (250, 1));
        simulator.AddSimulationModifier(tracking);

        //simulator.Solve();        
    }


    void TestBuskeVsSpringElasticityOnly() throw (Exception)
    {
        EXIT_IF_PARALLEL;
        double D = (4/3)*pow(10,-3);

        std::vector<CellPtr> buskecells = CreateCellsBuske();
        NodesOnlyMesh<3>* buskemesh = CreateMeshBuske();
        NodeBasedCellPopulationWithBuskeUpdate<3>* buskeCellPopulation = 
            new NodeBasedCellPopulationWithBuskeUpdate<3>(*buskemesh, buskecells);
                
        buskeCellPopulation->SetUseVaryingRadii(true);
        buskeCellPopulation->SetUseVariableRadii(true);
        buskeCellPopulation->SetDampingConstantIntercell(5*pow(10,10));
        buskeCellPopulation->SetDampingConstantMedium(1);
        buskeCellPopulation->SetDampingConstantVolume(400);

        //SETUP SIMULATION 1
        OffLatticeSimulation<3> simulator(*buskeCellPopulation);
        simulator.SetOutputDirectory("TestBuskeVsSpringElasticityOnly");
        simulator.SetEndTime(600.0);
        simulator.SetDt(1.0/250.0);
        simulator.SetSamplingTimestepMultiple(250);
        //ADD FORCES
        buskeCellPopulation->EnableElasticity(D);
        MAKE_PTR(BuskeElasticForce<3>, p_buske_elastic_force);
        p_buske_elastic_force->SetDeformationEnergyParameter(D);
        simulator.AddForce(p_buske_elastic_force);
        //RECORD CELL TRACKS
        MAKE_PTR_ARGS(CellTrackingOutput<3>, tracking, (250, 1));
        simulator.AddSimulationModifier(tracking);

        //simulator.Solve();        
    }


    void TestBuskeVsSpringElasticityPlusCompression() throw (Exception)
    {
        EXIT_IF_PARALLEL;
        double D = (4/3)*pow(10,-3);
        double K = 1000;

        std::vector<CellPtr> buskecells = CreateCellsBuske();
        NodesOnlyMesh<3>* buskemesh = CreateMeshBuske();
        NodeBasedCellPopulationWithBuskeUpdate<3>* buskeCellPopulation = 
            new NodeBasedCellPopulationWithBuskeUpdate<3>(*buskemesh, buskecells);
                
        buskeCellPopulation->SetUseVaryingRadii(true);
        buskeCellPopulation->SetUseVariableRadii(true);
        buskeCellPopulation->SetDampingConstantIntercell(5*pow(10,10));
        buskeCellPopulation->SetDampingConstantMedium(1);
        buskeCellPopulation->SetDampingConstantVolume(400);

        //SETUP SIMULATION 1
        OffLatticeSimulation<3> simulator(*buskeCellPopulation);
        simulator.SetOutputDirectory("TestBuskeVsSpringElasticityPlusCompression");
        simulator.SetEndTime(600.0);
        simulator.SetDt(1.0/250.0);
        simulator.SetSamplingTimestepMultiple(250);
        //ADD FORCES
        buskeCellPopulation->EnableElasticity(D);
        MAKE_PTR(BuskeElasticForce<3>, p_buske_elastic_force);
        p_buske_elastic_force->SetDeformationEnergyParameter(D);
        simulator.AddForce(p_buske_elastic_force);
        buskeCellPopulation->EnableCompression(K);
        MAKE_PTR(BuskeCompressionForce<3>, p_buske_comp_force);
        p_buske_comp_force->SetCompressionEnergyParameter(K);
        simulator.AddForce(p_buske_comp_force);
        //RECORD CELL TRACKS
        MAKE_PTR_ARGS(CellTrackingOutput<3>, tracking, (250, 1));
        simulator.AddSimulationModifier(tracking);

        //simulator.Solve();        
    }


    void TestBuskeVsSpringElasticityPlusAdhesionPlusCompression() throw (Exception)
    {
        EXIT_IF_PARALLEL;
        double D = (4/3)*pow(10,-3);
        double K = 1000;
        double eps = 200*pow(10,-6);

        std::vector<CellPtr> buskecells = CreateCellsBuske();
        NodesOnlyMesh<3>* buskemesh = CreateMeshBuske();
        NodeBasedCellPopulationWithBuskeUpdate<3>* buskeCellPopulation = 
            new NodeBasedCellPopulationWithBuskeUpdate<3>(*buskemesh, buskecells);
                
        buskeCellPopulation->SetUseVaryingRadii(true);
        buskeCellPopulation->SetUseVariableRadii(true);
        buskeCellPopulation->SetDampingConstantIntercell(5*pow(10,10));
        buskeCellPopulation->SetDampingConstantMedium(1);
        buskeCellPopulation->SetDampingConstantVolume(400);

        //SETUP SIMULATION 1
        OffLatticeSimulation<3> simulator(*buskeCellPopulation);
        simulator.SetOutputDirectory("TestBuskeVsSpringElasticityPlusAdhesionPlusCompression");
        simulator.SetEndTime(600.0);
        simulator.SetDt(1.0/250.0);
        simulator.SetSamplingTimestepMultiple(250);
        //ADD FORCES
        buskeCellPopulation->EnableElasticity(D);
        MAKE_PTR(BuskeElasticForce<3>, p_buske_elastic_force);
        p_buske_elastic_force->SetDeformationEnergyParameter(D);
        simulator.AddForce(p_buske_elastic_force);
        buskeCellPopulation->EnableCompression(K);
        MAKE_PTR(BuskeCompressionForce<3>, p_buske_comp_force);
        p_buske_comp_force->SetCompressionEnergyParameter(K);
        simulator.AddForce(p_buske_comp_force);
        buskeCellPopulation->EnableCompression(eps);
        MAKE_PTR(BuskeAdhesiveForce<3>, p_buske_adh_force);
        p_buske_adh_force->SetAdhesionEnergyParameter(eps);
        simulator.AddForce(p_buske_adh_force);
        //RECORD CELL TRACKS
        MAKE_PTR_ARGS(CellTrackingOutput<3>, tracking, (250, 1));
        simulator.AddSimulationModifier(tracking);

        simulator.Solve();        
    }


    

};

#endif //TESTSPRINGFORCEVSBUSKEMODEL_HPP
