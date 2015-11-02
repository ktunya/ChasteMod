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

#ifndef TEST3DVERTEXBASEDSIMULATION_HPP_
#define TEST3DVERTEXBASEDSIMULATION_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "NagaiHondaForce.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ShortAxisVertexBasedDivisionRule.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"


class Test3dVertexBasedSimulation : public AbstractCellBasedTestSuite
{
public:

    void Test3dVertexBased() throw(Exception)
    {
        // Create mutable vertex mesh
        std::vector<Node<3>*> nodes;
        nodes.push_back(new Node<3>(0, true, 0.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(1, true, 1.0, 0.0, 0.0));
        nodes.push_back(new Node<3>(2, true, 0.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(3, true, 0.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(4, true, 1.0, 1.0, 0.0));
        nodes.push_back(new Node<3>(5, true, 0.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(6, true, 1.0, 0.0, 1.0));
        nodes.push_back(new Node<3>(7, true, 1.0, 1.0, 1.0));
        nodes.push_back(new Node<3>(8, true, 0.5, 0.5, 1.5));

        std::vector<std::vector<Node<3>*> > nodes_faces(10);
        nodes_faces[0].push_back(nodes[0]);
        nodes_faces[0].push_back(nodes[2]);
        nodes_faces[0].push_back(nodes[4]);
        nodes_faces[0].push_back(nodes[1]);
        nodes_faces[1].push_back(nodes[4]);
        nodes_faces[1].push_back(nodes[7]);
        nodes_faces[1].push_back(nodes[5]);
        nodes_faces[1].push_back(nodes[2]);
        nodes_faces[2].push_back(nodes[7]);
        nodes_faces[2].push_back(nodes[6]);
        nodes_faces[2].push_back(nodes[1]);
        nodes_faces[2].push_back(nodes[4]);
        nodes_faces[3].push_back(nodes[0]);
        nodes_faces[3].push_back(nodes[3]);
        nodes_faces[3].push_back(nodes[5]);
        nodes_faces[3].push_back(nodes[2]);
        nodes_faces[4].push_back(nodes[1]);
        nodes_faces[4].push_back(nodes[6]);
        nodes_faces[4].push_back(nodes[3]);
        nodes_faces[4].push_back(nodes[0]);
        nodes_faces[5].push_back(nodes[7]);
        nodes_faces[5].push_back(nodes[6]);
        nodes_faces[5].push_back(nodes[3]);
        nodes_faces[5].push_back(nodes[5]);
        nodes_faces[6].push_back(nodes[6]);
        nodes_faces[6].push_back(nodes[7]);
        nodes_faces[6].push_back(nodes[8]);
        nodes_faces[7].push_back(nodes[6]);
        nodes_faces[7].push_back(nodes[8]);
        nodes_faces[7].push_back(nodes[3]);
        nodes_faces[8].push_back(nodes[3]);
        nodes_faces[8].push_back(nodes[8]);
        nodes_faces[8].push_back(nodes[5]);
        nodes_faces[9].push_back(nodes[5]);
        nodes_faces[9].push_back(nodes[8]);
        nodes_faces[9].push_back(nodes[7]);

        std::vector<VertexElement<2,3>*> faces;
        for (unsigned i=0; i<10; i++)
        {
            faces.push_back(new VertexElement<2,3>(i, nodes_faces[i]));
        }

        std::vector<VertexElement<2,3>*> faces_element_0, faces_element_1;
        std::vector<bool> orientations_0, orientations_1;
        for (unsigned i=0; i<6; i++)
        {
            faces_element_0.push_back(faces[i]);
            orientations_0.push_back(true);
        }
        for (unsigned i=6; i<10; i++)
        {
            faces_element_1.push_back(faces[i]);
            orientations_1.push_back(true);
        }
        faces_element_1.push_back(faces[5]);
        orientations_1.push_back(false);

        std::vector<VertexElement<3,3>*> elements;
        elements.push_back(new VertexElement<3,3>(0, faces_element_0, orientations_0));
        elements.push_back(new VertexElement<3,3>(1, faces_element_1, orientations_1));

        MutableVertexMesh<3,3> mesh(nodes, elements);

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, mesh.GetNumElements());

        // Create cell population
        VertexBasedCellPopulation<3>* p_cell_population = new VertexBasedCellPopulation<3>(mesh, cells);

        // Cells have been given birth times of 0 and -1.
        // Loop over them to run to time 0.0;
        for (AbstractCellPopulation<3>::Iterator cell_iter = p_cell_population->Begin();
             cell_iter != p_cell_population->End();
             ++cell_iter)
        {
            cell_iter->ReadyToDivide();
        }

        // Set up cell-based simulation
        OffLatticeSimulation<3> simulator(*p_cell_population);
        simulator.SetOutputDirectory("Test3dVertexBased");
        simulator.SetEndTime(1.0);

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<3>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // A NagaiHondaForce has to be used together with an AbstractTargetAreaModifier #2488
        MAKE_PTR(SimpleTargetAreaModifier<3>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Run simulation
        simulator.Solve();

        delete p_cell_population;     
    }
};

#endif /*TEST3DVERTEXBASEDSIMULATION_HPP_*/