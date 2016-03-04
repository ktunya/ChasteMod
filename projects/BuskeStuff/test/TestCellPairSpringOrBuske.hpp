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


#ifndef TESTCELLPAIRSPRINGORBUSKE_HPP_
#define TESTCELLPAIRSPRINGORBUSKE_HPP_

#include <cxxtest/cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

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
#include "GeneralisedLinearSpringForce.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PetscTools.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "CellTrackingOutput.hpp"

class TestCellPairSpringOrBuske : public AbstractCellBasedTestSuite {

public:

    void TestCellPair() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // SETTINGS -------------------------------------------------------------------
        //-----------------------------------------------------------------------------

        std::string outDir("SpringNewMu");

        enum Force {BUSKE, SPRING}; 
        Force ForceLaw = SPRING;

        // buske only
        bool adhesionOn = true;
        bool elasticityOn = true;
        bool compressionOn = true;
        bool varyingRadiiOn = true;

        double eps = 2592;                        //200*pow(10,-6);        
        double d = 1.0288*pow(10,-4);             //(4.0/3.0)*pow(10,-3);         
        double k = 12960;                         //1000;                     
        double dampingIntercell = 180;            //5 * pow(10, 10);
        double dampingMedium = 11520;             //3.2;
        double dampingVolume = 1440000;           //400; 


        // spring only
        double meinekeSpringConstant = 14152;
        double adhesionDecay = 5;
        double damping = 11520;

        double initialSepMicrons = 1;
        double R1Microns = 3;
        double R2Microns = 3;

        //----------------------------------------------------------------------------
        //----------------------------------------------------------------------------

        // CREATE CELL POPULATION
        NodeBasedCellPopulation<3>* cellPopulation;
        
        if(ForceLaw == BUSKE){

            std::vector<CellPtr> buskecells = CreateCellsBuske(R1Microns, R2Microns);
            NodesOnlyMesh<3>* buskemesh = CreateMeshBuske(R1Microns, R2Microns, initialSepMicrons);
            NodeBasedCellPopulationWithBuskeUpdate<3>* bCellPopulation = 
                new NodeBasedCellPopulationWithBuskeUpdate<3>(*buskemesh, buskecells);

            bCellPopulation->SetUseVaryingRadii( varyingRadiiOn );
            bCellPopulation->SetDampingConstantIntercell( dampingIntercell );  
            bCellPopulation->SetDampingConstantMedium( dampingMedium );                 
            bCellPopulation->SetDampingConstantVolume( dampingVolume );                 
            if( adhesionOn ){
                bCellPopulation->EnableAdhesion(eps);
            }
            if( compressionOn ){
                bCellPopulation->EnableCompression(k);
            }
            if( elasticityOn ){
                bCellPopulation->EnableElasticity(d);
            }
            cellPopulation = static_cast<NodeBasedCellPopulation<3>*>(bCellPopulation);

        }else if(ForceLaw == SPRING){

            std::vector<CellPtr> springcells = CreateCellsSpring(R1Microns, R2Microns);
            NodesOnlyMesh<3>* springmesh = CreateMeshSpring(R1Microns, R2Microns, initialSepMicrons);
            cellPopulation = new NodeBasedCellPopulation<3>(*springmesh, springcells);
            
            cellPopulation->SetDampingConstantNormal(damping);
            cellPopulation->SetUseVariableRadii(true);

        }else{
        }

        cellPopulation->SetUseVariableRadii( true );
        cellPopulation->SetAbsoluteMovementThreshold( 0.5 );

        //----------------------------------------------------------------------------
        //----------------------------------------------------------------------------

        // SETUP SIMULATION
        OffLatticeSimulation<3> simulator(*cellPopulation);
        simulator.SetOutputDirectory(outDir.c_str());
        simulator.SetEndTime(0.5);
        simulator.SetDt(1.0/3600.0);
        simulator.SetSamplingTimestepMultiple(1);

        // ADD FORCES
        if(ForceLaw == BUSKE){
            
            if(compressionOn){
                MAKE_PTR(BuskeCompressionForce<3>, pCompression);
                pCompression->SetCompressionEnergyParameter(k);
                simulator.AddForce(pCompression);
            }
            if(elasticityOn){
                MAKE_PTR(BuskeElasticForce<3>, pElastic);
                pElastic->SetDeformationEnergyParameter(d);
                simulator.AddForce(pElastic);
            }
            if(adhesionOn){
                MAKE_PTR(BuskeAdhesiveForce<3>, pAdhesive);
                pAdhesive->SetAdhesionEnergyParameter(eps);
                simulator.AddForce(pAdhesive);
            }

        }else if(ForceLaw == SPRING){
            MAKE_PTR(GeneralisedLinearSpringForce<3>, pGLS);
            pGLS->SetAdhesionDecayAlpha(adhesionDecay);
            pGLS->SetMeinekeSpringStiffness(meinekeSpringConstant);
            simulator.AddForce(pGLS);
        }else{
        }

        // RECORD CELL POSITIONS
        MAKE_PTR_ARGS(CellTrackingOutput<3>, pTracking, (1, 1));
        simulator.AddSimulationModifier(pTracking);
        
        // SOLVE
        simulator.Solve();
    }

private:

    NodesOnlyMesh<3>* CreateMeshBuske(double R1, double R2, double initialSep){

        double startingRadius1 = R1;
        double startingRadius2 = R2;

        c_vector<double, 3> pos1;
        pos1[0]=0; pos1[1]=0; pos1[2]=0;
        c_vector<double, 3> pos2;
        pos2[0]=0; pos2[1]=initialSep; pos2[2]=0;

        NodesOnlyMesh<3>* mesh = CreateMesh(pos1, pos2, startingRadius1, startingRadius2);
        return(mesh);     
    };

    NodesOnlyMesh<3>* CreateMeshSpring(double R1, double R2, double initialSep){

        double relaxedRadius1 = R1;
        double relaxedRadius2 = R2;
        double startingRadius1 = R1;
        double startingRadius2 = R2;

        c_vector<double, 3> pos1;
        pos1[0]=0; pos1[1]=0; pos1[2]=0;
        c_vector<double, 3> pos2;
        pos2[0]=0; pos2[1]=initialSep; pos2[2]=0;

        NodesOnlyMesh<3>* mesh = CreateMesh(pos1, pos2, startingRadius1, startingRadius2);
        return(mesh);      
    };

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
    };

    std::vector<CellPtr> CreateCellsBuske(double R1, double R2){

        double relaxedRadius1 = R1;
        double relaxedRadius2 = R2;
        double startingRadius1 = R1;
        double startingRadius2 = R2;

        std::vector<CellPtr> Cells = CreateCells(relaxedRadius1, relaxedRadius2, startingRadius1, startingRadius2, false);
        return(Cells);
    }

    std::vector<CellPtr> CreateCellsSpring(double R1, double R2){

        double relaxedRadius1 = R1;
        double relaxedRadius2 = R2;
        double startingRadius1 = R1;
        double startingRadius2 = R2;

        std::vector<CellPtr> Cells = CreateCells(relaxedRadius1, relaxedRadius2, startingRadius1, startingRadius2, true);
        return(Cells);
    }

    std::vector<CellPtr> CreateCells(double relaxedRadius1, double relaxedRadius2, double startingRadius1, double startingRadius2, bool setRadius){

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
                new_cell->GetCellData()->SetItem("IsBuskeKnot", 0);
            }else{
                new_cell->GetCellData()->SetItem("RelaxedRadius", relaxedRadius2);
                new_cell->GetCellData()->SetItem("Radius", startingRadius2);
                new_cell->GetCellData()->SetItem("IsBuskeKnot", 0);
            }
            new_cell->GetCellData()->SetItem("volume",0.0); //UNUSED
            cells.push_back(new_cell);
        }  
        return( cells );
    };
};

#endif //TESTCELLPAIRSPRINGORBUSKE_HPP_
