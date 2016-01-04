
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

#ifndef TESTBUSKEDISTALBOUNDARY_HPP_
#define TESTBUSKEDISTALBOUNDARY_HPP_

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "NodeBasedCellPopulationWithBuskeUpdate.hpp"
#include "StochasticDurationCellCycleModel.hpp"
#include "BuskeAdhesiveForce.hpp"
#include "BuskeElasticForce.hpp"
#include "BuskeCompressionForce.hpp"
#include "BuskeDistalPotentialBoundaryCondition.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestBuskeDistalBoundary : public AbstractCellBasedWithTimingsTestSuite
{
public:

    /**
     * Create a simulation of a NodeBasedCellPopulation with all Buske forces.
     * Add a boundary forming a crypt/distal arm shaped potential well around the population 
     * and check for reasonable levels of containment / migration. 
     */

    void TestBuskeDistal(){
 
        EXIT_IF_PARALLEL;   

        // Create nodes and mesh
        double length = 248;
        double radius = 11.3;
        double initialCellRowSpacing = 5;
        double cellsPerRow = 10.0;
        NodesOnlyMesh<3>* mesh = CreateDistalArmCellLocations(length, radius, initialCellRowSpacing, cellsPerRow);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(StemCellProliferativeType, pStemType);
        CellsGenerator<StochasticDurationCellCycleModel, 3> cellsGenerator;
        cellsGenerator.GenerateBasicRandom(cells, mesh->GetNumNodes(), pStemType);
        NodeBasedCellPopulationWithBuskeUpdate<3> cellPopulation(*mesh, cells);

        for(typename AbstractCellPopulation<3,3>::RealCellsIterator cell_iter = cellPopulation.Begin();
            cell_iter != cellPopulation.End();
            ++cell_iter){
        	cell_iter->GetCellData()->SetItem("Radius",2.8);
        	cell_iter->GetCellData()->SetItem("RelaxedRadius",2.8);
        	cell_iter->GetCellData()->SetItem("IsBuskeKnot",0);
        }


        cellPopulation.SetUseVaryingRadii(true);
        cellPopulation.SetUseVariableRadii(true);
        cellPopulation.SetMeinekeDivisionSeparation(0.99);
        cellPopulation.SetDampingConstantIntercell(5 * pow(10, 10));  // 5e10
        cellPopulation.SetDampingConstantMedium(3.2);                 // 3.2
        cellPopulation.SetDampingConstantVolume(400);                 // 400
        cellPopulation.SetAbsoluteMovementThreshold(0.5);  
        double eps = 200*pow(10,-6);         // 200*pow(10,-6)
        double d = (4/3)*pow(10,-3);         // (4/3)*pow(10,-3)
        double k = 1000;                     // 1000;
        cellPopulation.EnableAdhesion(eps);
        cellPopulation.EnableCompression(k);
        cellPopulation.EnableElasticity(d);


        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator(cellPopulation, false, true, true);
        simulator.SetOutputDirectory("TestBuskeDistal");
        simulator.SetDt(1/250.0);
        simulator.SetSamplingTimestepMultiple(25);
        simulator.SetEndTime(500);


        std::cout << "Add boundaries" << std::endl;
        // Add the distal boundary 
        MAKE_PTR_ARGS(BuskeDistalPotentialBoundaryCondition<3>, pBoundary, (length, radius, 1000));
        simulator.AddForce(pBoundary);

        // Add a cell killer at the end of the tube
        c_vector<double, 3> point;
        point[0] = length;
        point[1] = 0;
        point[2] = 0;
        c_vector<double, 3> normal;
        normal[0] = 1;
        normal[1] = 0;
        normal[2] = 0;
        MAKE_PTR_ARGS(PlaneBasedCellKiller<3>, pKiller, (&cellPopulation, point, normal));
        simulator.AddCellKiller(pKiller);    


    	std::cout << "Add forces" << std::endl;
        // Create force laws and pass them to the simulation
        MAKE_PTR(BuskeCompressionForce<3>, pCompression);
        MAKE_PTR(BuskeElasticForce<3>, pElastic);
        MAKE_PTR(BuskeAdhesiveForce<3>, pAdhesive);

        pCompression->SetCompressionEnergyParameter(k);
        pElastic->SetDeformationEnergyParameter(d);
        pAdhesive->SetAdhesionEnergyParameter(eps);
        
        simulator.AddForce(pCompression);
        simulator.AddForce(pElastic);
        simulator.AddForce(pAdhesive);


        std::cout << "Begin solve" << std::endl;
        simulator.Solve();
    }


private:

    NodesOnlyMesh<3>* CreateDistalArmCellLocations(double armLength, double armRadius, double rowSpacing, int cellsPerRow){

        double angleSpacing = 2*M_PI/cellsPerRow;
        double cellRadius = 2.8;
        
        std::vector< Node<3>* > nodes;
        unsigned nodeIndex = 0;

        for(double x = -(armRadius-cellRadius); x < armLength/5.0; x += rowSpacing){
            for(double angle = 0; angle < 2*M_PI; angle += angleSpacing){
                Node<3>* n = new Node<3>(nodeIndex, false, x, (armRadius-cellRadius)*(cos(angle)), (armRadius-cellRadius)*(sin(angle)));
                n->SetRadius(cellRadius);
                nodes.push_back(n);
                nodeIndex++;
            }
        }

        NodesOnlyMesh<3>* distalMesh = new NodesOnlyMesh<3>;
        distalMesh->ConstructNodesWithoutMesh(nodes, 10);
        for (unsigned i = 0; i<nodes.size(); i++)
        {
           delete nodes[i];
        }
        return distalMesh;
    }
    
};

#endif /*TESTBUSKEDISTALBOUNDARY_HPP_*/
