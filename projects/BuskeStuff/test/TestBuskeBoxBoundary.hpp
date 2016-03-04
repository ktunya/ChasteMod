
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

#ifndef TESTBUSKEBOXBOUNDARY_HPP_
#define TESTBUSKEBOXBOUNDARY_HPP_

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "NodeBasedCellPopulationWithBuskeUpdate.hpp"
#include "FixedDivisionTimingsCellCycleModel.hpp"
#include "BuskeAdhesiveForce.hpp"
#include "BuskeElasticForce.hpp"
#include "BuskeCompressionForce.hpp"
#include "BuskeKnotForce.hpp"
#include "CellTrackingOutput.hpp"
#include "BuskePlaneKnotBoundaryCondition.hpp"
#include "BuskePlanePotentialBoundaryCondition.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestBuskeBoxBoundary : public AbstractCellBasedWithTimingsTestSuite
{
public:

    /**
     * Create a simulation of a NodeBasedCellPopulation with all Buske forces.
     * Add 6 plane boundaries forming a cube around the population and check for containment.
     */

    void TestBuskeBoundary(){
 
        EXIT_IF_PARALLEL;   

        // Create nodes and mesh
        std::vector< Node<3>* > nodes;
        Node<3>* node1 = new Node<3>(0, false, 0, 1, 0);
        node1->SetRadius(3);
        nodes.push_back(node1);
        Node<3>* node2 = new Node<3>(0, false, 0, 1.9, 0);
        node2->SetRadius(3);
        nodes.push_back(node2);

        NodesOnlyMesh<3>* realCellsMesh = new NodesOnlyMesh<3>;
        realCellsMesh->ConstructNodesWithoutMesh(nodes, 50);


        // Create cells
        std::vector<CellPtr> cells;
        cells.reserve(2);
        MAKE_PTR(StemCellProliferativeType, pStemType);

        for (unsigned i=0; i<2; i++)
        {
            FixedDivisionTimingsCellCycleModel* p_cell_cycle_model = new FixedDivisionTimingsCellCycleModel(8);
            p_cell_cycle_model->SetDimension(3);

            boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
            CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));

            p_cell->SetCellProliferativeType(pStemType);
            double birth_time = -p_cell_cycle_model->GetAverageStemCellCycleTime()*RandomNumberGenerator::Instance()->ranf();
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        NodeBasedCellPopulationWithBuskeUpdate<3> cellPopulation(*realCellsMesh, cells);

        for(typename AbstractCellPopulation<3,3>::Iterator cell_iter = cellPopulation.Begin();
            cell_iter != cellPopulation.End();
            ++cell_iter){
        	cell_iter->GetCellData()->SetItem("Radius",3);
        	cell_iter->GetCellData()->SetItem("RelaxedRadius",3);
        	cell_iter->GetCellData()->SetItem("IsBuskeKnot",0);
        }


        cellPopulation.SetUseVaryingRadii(true);
        cellPopulation.SetUseVariableRadii(true);

        cellPopulation.SetMeinekeDivisionSeparation(0.95);
        cellPopulation.SetDampingConstantIntercell(180);  
        cellPopulation.SetDampingConstantMedium(11520);                
        cellPopulation.SetDampingConstantVolume(1440000);                
        cellPopulation.SetAbsoluteMovementThreshold(2);  

        double eps = 2592;       
        double d = 1.0288*pow(10,-4);         
        double k = 12960;                    
        cellPopulation.EnableAdhesion(eps);
        //cellPopulation.EnableCompression(k);
        cellPopulation.EnableElasticity(d);



        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator(cellPopulation, false, true);
        simulator.SetOutputDirectory("NewBoxBoundaryNoCompression");
        simulator.SetDt(1/3000.0);
        simulator.SetSamplingTimestepMultiple(300);
        simulator.SetEndTime(40);

        std::cout << "Add boundaries" << std::endl;

        // Add the box boundaries 
        double N = 6;
        double b[18] = {N,-N,-N,  -N,-N,-N,   -N,N,-N,  -N,-N,-N,  -N,-N,N,  -N,-N,-N};
        double t[18] = {N,N,N,    -N,N,N,     N,N,N,    N,-N,N,    N,N,N,    N,N,-N};
        double n[18] = {1,0,0,    -1,0,0,     0,1,0,    0,-1,0,    0,0,1,    0,0,-1};

        for(int side=0; side<6; side++){

        	//c_vector<double, 3> bottomLeft;
        	//bottomLeft[0]=b[side*3];  bottomLeft[1]=b[side*3+1];  bottomLeft[2]=b[side*3+2];
        	
        	//c_vector<double, 3> topRight;
			//topRight[0]=t[side*3];  topRight[1]=t[side*3+1];  topRight[2]=t[side*3+2];
			
			c_vector<double, 3> pointOnPlane;
			pointOnPlane[0]=t[side*3];  pointOnPlane[1]=t[side*3+1];  pointOnPlane[2]=t[side*3+2];
			
			c_vector<double, 3> normal;
			normal[0]=n[side*3];  normal[1]=n[side*3+1];  normal[2]=n[side*3+2];
			
			//double knotSpacing = 1;
        	
            //MAKE_PTR_ARGS(PlaneBoundaryCondition<3>, pBoundary, (&cellPopulation, pointOnPlane, normal));
        	//MAKE_PTR_ARGS(BuskePlaneKnotBoundaryCondition<3>, pBoundary, (&cellPopulation, bottomLeft, topRight, normal, knotSpacing));
        	MAKE_PTR_ARGS(BuskePlanePotentialBoundaryCondition<3>, pBoundary, (pointOnPlane, normal, 1.0, 1.0));
        	simulator.AddForce(pBoundary);
    	}

    	std::cout << "Add forces" << std::endl;

        // Create force laws and pass them to the simulation
        MAKE_PTR(BuskeCompressionForce<3>, pCompression);
        MAKE_PTR(BuskeElasticForce<3>, pElastic);
        MAKE_PTR(BuskeAdhesiveForce<3>, pAdhesive);
        //MAKE_PTR(BuskeKnotForce<3>, pKnot);

        pCompression->SetCompressionEnergyParameter(k);
        pElastic->SetDeformationEnergyParameter(d);
        pAdhesive->SetAdhesionEnergyParameter(eps);
        //pKnot->SetMaxInteractionEnergy(pow(10,-11));
        //pKnot->SetThresholdAdhesionRatio(0.95);
        
        //simulator.AddForce(pCompression);
        simulator.AddForce(pElastic);
        simulator.AddForce(pAdhesive);
        //simulator.AddForce(pKnot);

        MAKE_PTR_ARGS(CellTrackingOutput<3>, pTracking, (3000, 1));
        simulator.AddSimulationModifier(pTracking);

        std::cout << "Begin solve" << std::endl;

        simulator.Solve();
    
        // Cleanup
        delete node1;
        delete node2;
    }

    
};

#endif /*TESTBUSKEBOXBOUNDARY_HPP_*/
