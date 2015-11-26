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

#ifndef TESTFIXEDDIVISIONTIMINGS_HPP_
#define TESTFIXEDDIVISIONTIMINGS_HPP_


#include <cxxtest/TestSuite.h>


#include "CellBasedSimulationArchiver.hpp"
#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "FixedDivisionTimingsCellCycleModel.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "LogFile.hpp"
#include "WildTypeCellMutationState.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"
#include "PetscSetupAndFinalize.hpp"


class TestFixedDivisionTimings: public AbstractCellBasedWithTimingsTestSuite {

public:

void TestNodeBasedWithFixedDivisions(){

	EXIT_IF_PARALLEL; 

	double divisionInterval = 3.7;

    // Create a simple mesh
    int num_cells_depth = 5;
    int num_cells_width = 5;
    HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
    TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

    // Convert this to a NodesOnlyMesh
    NodesOnlyMesh<2> mesh;
    mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);
    int numCells = mesh.GetNumNodes();

    std::vector<CellPtr> rCells;
    rCells.clear();
    rCells.reserve(numCells);

    // Create cells
    for (unsigned i=0; i<numCells; i++)
    {
        FixedDivisionTimingsCellCycleModel* p_cell_cycle_model = new FixedDivisionTimingsCellCycleModel(divisionInterval);
        p_cell_cycle_model->SetDimension(2);

        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));
        p_cell->SetCellProliferativeType(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());
        double birth_time = -p_cell_cycle_model->GetAverageStemCellCycleTime()*RandomNumberGenerator::Instance()->ranf();
        p_cell->SetBirthTime(birth_time);
        rCells.push_back(p_cell);
    }

    // Create a node based cell population
    NodeBasedCellPopulation<2> node_based_cell_population(mesh, rCells);

    // Set up cell-based simulation
    OffLatticeSimulation<2> simulator(node_based_cell_population);
    simulator.SetOutputDirectory("TestNodeBasedWithFixedDivisions");
    simulator.SetDt(0.005);
    simulator.SetEndTime(3.703);

    // Create a force law and pass it to the simulation
    MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
    p_linear_force->SetCutOffLength(1.5);
    simulator.AddForce(p_linear_force);
    
    // Solve
    simulator.Solve();

    TS_ASSERT(simulator.rGetCellPopulation().GetNumNodes() == 2*numCells );

    simulator.SetEndTime(7.403);
    simulator.Solve();

    TS_ASSERT(simulator.rGetCellPopulation().GetNumNodes() == 4*numCells );
} 

};

#endif /*TESTFIXEDDIVISIONTIMINGS_HPP_*/