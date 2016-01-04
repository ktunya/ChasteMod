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
#include "NodeBasedCellPopulationWithBuskeUpdate.hpp"

#include "ReplicatableVector.hpp"
#include "OdeLinearSystemSolver.hpp"

template<unsigned DIM>
NodeBasedCellPopulationWithBuskeUpdate<DIM>::NodeBasedCellPopulationWithBuskeUpdate(NodesOnlyMesh<DIM>& rMesh,
                                      std::vector<CellPtr>& rCells,
                                      const std::vector<unsigned> locationIndices,
                                      bool deleteMesh)
    : NodeBasedCellPopulation<DIM>(rMesh, rCells, locationIndices, deleteMesh),
    numKnots(0)
{
    Eps = DOUBLE_UNSET;
    K = DOUBLE_UNSET;
    D = DOUBLE_UNSET;
    dampingConstantIntercell = DOUBLE_UNSET;
    dampingConstantVolume = DOUBLE_UNSET;
    dampingConstantMedium = DOUBLE_UNSET;
    allowRadiusVariation = true;
    movementDisabled = false;

    adhesionEnabled = false;
    elasticityEnabled = false;
    compressionEnabled = false;
}


template<unsigned DIM>
NodeBasedCellPopulationWithBuskeUpdate<DIM>::NodeBasedCellPopulationWithBuskeUpdate(NodesOnlyMesh<DIM>& rMesh)
    : NodeBasedCellPopulation<DIM>(rMesh)
{
    // No Validate() because the cells are not associated with the cell population yet in archiving
}


template<unsigned DIM>
void NodeBasedCellPopulationWithBuskeUpdate<DIM>::EnableAdhesion(double eps){
    Eps = eps;
    adhesionEnabled = true;
    std::cout << "Adhesion enabled. Ensure that a BuskeAdhesiveForce is added to the simulator" << std::endl;
};
template<unsigned DIM>
void NodeBasedCellPopulationWithBuskeUpdate<DIM>::EnableElasticity(double d){
    D = d;
    elasticityEnabled = true;
    std::cout << "Elasticity enabled. Ensure that a BuskeElasticForce is added to the simulator" << std::endl;
};
template<unsigned DIM>
void NodeBasedCellPopulationWithBuskeUpdate<DIM>::EnableCompression(double k){
    K = k;
    compressionEnabled = true;
    std::cout << "Compression enabled. Ensure that a BuskeCompressionForce is added to the simulator" << std::endl;
};



template<unsigned DIM>
const double NodeBasedCellPopulationWithBuskeUpdate<DIM>::GetEps() const { return Eps; }
template<unsigned DIM>
const double NodeBasedCellPopulationWithBuskeUpdate<DIM>::GetK() const { return K; }
template<unsigned DIM>
const double NodeBasedCellPopulationWithBuskeUpdate<DIM>::GetD() const { return D; }
template<unsigned DIM>
const double NodeBasedCellPopulationWithBuskeUpdate<DIM>::GetDampingConstantIntercell() const { return dampingConstantIntercell; }
template<unsigned DIM>
const double NodeBasedCellPopulationWithBuskeUpdate<DIM>::GetDampingConstantVolume() const { return dampingConstantVolume; }
template<unsigned DIM>
const double NodeBasedCellPopulationWithBuskeUpdate<DIM>::GetDampingConstantMedium() const { return dampingConstantMedium; }



template<unsigned DIM>
void NodeBasedCellPopulationWithBuskeUpdate<DIM>::SetDampingConstantIntercell(double inter){
    dampingConstantIntercell = inter;
};
template<unsigned DIM>
void NodeBasedCellPopulationWithBuskeUpdate<DIM>::SetDampingConstantVolume(double vol){
    dampingConstantVolume = vol;
};
template<unsigned DIM>
void NodeBasedCellPopulationWithBuskeUpdate<DIM>::SetDampingConstantMedium(double med){
    dampingConstantMedium = med;
};

template<unsigned DIM>
const int NodeBasedCellPopulationWithBuskeUpdate<DIM>::GetNumKnots() const{
    return numKnots;
};
template<unsigned DIM>
void NodeBasedCellPopulationWithBuskeUpdate<DIM>::SetNumKnots(int nKnots){
    numKnots = nKnots;
};


template<unsigned DIM>
const bool NodeBasedCellPopulationWithBuskeUpdate<DIM>::GetUseVaryingRadii() const{
    return allowRadiusVariation;
};
template<unsigned DIM>
void NodeBasedCellPopulationWithBuskeUpdate<DIM>::SetUseVaryingRadii(bool radSetting){
    allowRadiusVariation = radSetting;
};

template<unsigned DIM>
const bool NodeBasedCellPopulationWithBuskeUpdate<DIM>::GetMovementDisabled() const{
    return movementDisabled;
};
template<unsigned DIM>
void NodeBasedCellPopulationWithBuskeUpdate<DIM>::SetDisableMovement(bool movSetting){
    movementDisabled = movSetting;
};


template<unsigned DIM>
void NodeBasedCellPopulationWithBuskeUpdate<DIM>::UpdateNodeLocations(double dt)
{

    int cellCount = this->GetNumNodes() - numKnots;

    // Declare solver and give the size of the system and timestep
    unsigned system_size;
    if(allowRadiusVariation && !movementDisabled){
        system_size = cellCount*(DIM + 1); //3 position eqns and 1 radius eqn per node
    }else if(!allowRadiusVariation && !movementDisabled){
        system_size = cellCount*DIM; //3 position eqns only
    }else if(allowRadiusVariation && movementDisabled){
        system_size = cellCount; // 1 radius eqn per node
    }else{
        EXCEPTION("Neither radial changes nor movement enabled in NodeBasedCellPopulationWithBuskeUpdate");
    }
    OdeLinearSystemSolver solver(system_size, dt);


    // Set up the matrix, initial condition and forcing
    Mat& r_matrix = solver.rGetLhsMatrix();
    Vec initial_condition = PetscTools::CreateAndSetVec(system_size, 0.0);
    Vec& r_vector = solver.rGetForceVector();


    // Iterate over all nodes associated with real cells, and construct the LHS matrix.
    for (typename AbstractCellPopulation<DIM>::RealCellsIterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {

        if(cell_iter->GetCellData()->GetItem("IsBuskeKnot")==1){
            //Knot, ignore
        }else{

            // Get global and local node indices, and the node itself
            unsigned global_node_index = this->GetLocationIndexUsingCell((*cell_iter));
            unsigned node_index = this->rGetMesh().SolveNodeMapping(global_node_index);
            unsigned petscvec_index = node_index;
            if(petscvec_index >= cellCount){
                petscvec_index -= numKnots;
            }
            Node<DIM>* p_node_i = this->GetNode(global_node_index);


            // Get the location and radius
            c_vector<double, DIM> node_i_location = p_node_i->rGetLocation();
            double radius_of_cell_i = this->GetCellUsingLocationIndex(global_node_index)->GetCellData()->GetItem("Radius");


            // Loop over neighbouring cells and calculate surface contact areas, along with contributions to G, VC and dVCdR
            double G = 0;
            double dVCdR = 0;
            double VC = 0;
            std::set<unsigned> neighbouring_node_indices = this->GetNeighbouringNodeIndices(global_node_index);

            for (std::set<unsigned>::iterator iter = neighbouring_node_indices.begin();
                 iter != neighbouring_node_indices.end();
                 ++iter)
            {
                // Get neighbour node indices and node
                unsigned neighbour_node_global_index = *iter;
                unsigned neighbour_node_index = this->rGetMesh().SolveNodeMapping(neighbour_node_global_index);
                unsigned neighbour_petscvec_index = neighbour_node_index;
                if(neighbour_petscvec_index >= cellCount){
                    neighbour_petscvec_index -= numKnots;
                }
                Node<DIM>* p_node_j = this->GetNode(neighbour_node_global_index);

                double isKnot = (this->GetCellUsingLocationIndex(neighbour_node_global_index))->GetCellData()->GetItem("IsBuskeKnot");
                if(isKnot == 1){
                    //Ignore knots
                }else{

                    // Get the neighbour location, and distance between the two cells
                    c_vector<double, DIM> node_j_location = p_node_j->rGetLocation();
                    c_vector<double, DIM> unit_vector = node_j_location - node_i_location;
                    double dij = norm_2(unit_vector);
                    unit_vector /= dij;


                    // Get the radius of the neighbouring node
                    double radius_of_cell_j =
                            (this->GetCellUsingLocationIndex(neighbour_node_global_index))->GetCellData()->GetItem("Radius");


                    // Calculate Aij, the contact area with this neighbour cell
                    double Aij = 0.0;
                    if (dij < radius_of_cell_i + radius_of_cell_j)
                    {

                        // Compute the distance to the contact surface
                        double xij = (radius_of_cell_i*radius_of_cell_i - radius_of_cell_j*radius_of_cell_j + dij*dij)/(2*dij);
                        Aij = M_PI*(radius_of_cell_i*radius_of_cell_i - xij*xij);

                        //std::cout << "x for cell id " << cell_iter->GetCellId() << "\t" << xij << std::endl;


                        // Put the contribution from the sum term in (A7) of Buske's paper into the matrix...
                        if(!movementDisabled) {
                            for (unsigned i = 0; i < DIM; i++) {
                                PetscMatTools::AddToElement(r_matrix, DIM * neighbour_petscvec_index + i,
                                                            DIM * neighbour_petscvec_index + i, -dampingConstantIntercell * Aij);
                                PetscMatTools::AddToElement(r_matrix, DIM * petscvec_index + i, DIM * petscvec_index + i,
                                                            dampingConstantIntercell * Aij);
                            }
                        }
                        // If radial variations are active, add a similar term for the radial equations:
                        if(allowRadiusVariation && !movementDisabled){
                            PetscMatTools::AddToElement(r_matrix, DIM * cellCount + neighbour_petscvec_index,
                                                        DIM * cellCount + neighbour_petscvec_index,
                                                        dampingConstantIntercell * Aij);
                            PetscMatTools::AddToElement(r_matrix, DIM * cellCount + petscvec_index,
                                                        DIM * cellCount + petscvec_index, dampingConstantIntercell * Aij);
                        }else if(allowRadiusVariation && movementDisabled){
                            PetscMatTools::AddToElement(r_matrix, neighbour_petscvec_index, neighbour_petscvec_index,
                                                        dampingConstantIntercell * Aij);
                            PetscMatTools::AddToElement(r_matrix, petscvec_index, petscvec_index, dampingConstantIntercell * Aij);
                        }


                        //Now compute contributions to G as required:
                        if(allowRadiusVariation) {

                            double dx = radius_of_cell_i / dij;

                            if(adhesionEnabled) {
                                double dWAdR = 2 * Eps * M_PI * (radius_of_cell_i - xij * dx);
                                G += dWAdR;

                                //std::cout << "dA for cell id " << cell_iter->GetCellId() << "\t" << dWAdR << std::endl;

                            }

                            if(elasticityEnabled) {
                                double sub1 = pow((radius_of_cell_i + radius_of_cell_j - dij), 0.5);
                                double sub2 = pow((radius_of_cell_i * radius_of_cell_j / (radius_of_cell_i + radius_of_cell_j)),
                                              0.5);
                                double dWDdR = (pow(sub1, 3) / D) * sub2 +
                                               (1.0 / (5.0 * D)) * pow(sub1, 5) * (1.0 / sub2) *
                                               ((radius_of_cell_j * radius_of_cell_j) /
                                                pow((radius_of_cell_i + radius_of_cell_j), 2));
                                G -= dWDdR;

                                //std::cout << "dD for cell id " << cell_iter->GetCellId() << "\t" << dWDdR << std::endl;
                            }

                            if(compressionEnabled) {
                                VC += (M_PI / 3.0) * pow((radius_of_cell_i - xij), 2) * (2 * radius_of_cell_i - xij);
                                dVCdR +=
                                        (M_PI / 3.0) * (2 * (radius_of_cell_i - xij) * (1 - dx) * (2 * radius_of_cell_i - xij) +
                                                        pow(radius_of_cell_i - xij, 2) * (2 - dx));

                            }
                        }
                    }
                }
            }
            if(allowRadiusVariation && compressionEnabled) {
                double targetRadius = this->GetCellUsingLocationIndex(global_node_index)->GetCellData()->GetItem("RelaxedRadius");
                double VT = (4.0 / 3.0) * M_PI * pow(targetRadius, 3);
                double dVAdR = 4 * M_PI * pow(radius_of_cell_i, 2) - dVCdR;
                double VA = (4.0 / 3.0) * M_PI * pow(radius_of_cell_i, 3) - VC;
                double dWKdR = -(K / VT) * (VT-VA) * dVAdR;
                G -= dWKdR;

                //std::cout << "dVC for cell id " << cell_iter->GetCellId() << "\t" << dVCdR << std::endl;
                //std::cout << "VT for cell id " << cell_iter->GetCellId() << "\t" << VT << std::endl;
                //std::cout << "VC for cell id " << cell_iter->GetCellId() << "\t" << VC << std::endl;
                //std::cout << "VA for cell id " << cell_iter->GetCellId() << "\t" << VA << std::endl;
                //std::cout << "dVA for cell id " << cell_iter->GetCellId() << "\t" << dVAdR << std::endl;
            }

            //std::cout << "G for cell id " << cell_iter->GetCellId() << "\t" << G << std::endl;


            // This is the contribution NOT from the sum in (A7)...
            if(!movementDisabled) {

                for (unsigned i = 0; i < DIM; i++) {
                    PetscMatTools::AddToElement(r_matrix, DIM * petscvec_index + i, DIM * petscvec_index + i, dampingConstantMedium);
                }

            }
            if(allowRadiusVariation && !movementDisabled) {
                PetscMatTools::AddToElement(r_matrix,
                                            DIM * cellCount + petscvec_index,
                                            DIM * cellCount + petscvec_index, dampingConstantVolume);
            }else if(allowRadiusVariation && movementDisabled){
                PetscMatTools::AddToElement(r_matrix, petscvec_index, petscvec_index, dampingConstantVolume);
            }


            // Construct initial condition and RHS vectors
            // Note that we define these vectors before setting them as otherwise the profiling build will break (see #2367)
            c_vector<double, DIM> current_location;
            c_vector<double, DIM> forces;
            current_location = this->GetNode(global_node_index)->rGetLocation();
            forces = this->GetNode(global_node_index)->rGetAppliedForce();
            double radius = this->GetCellUsingLocationIndex(global_node_index)->GetCellData()->GetItem("Radius");

            if(!movementDisabled) {
                for (unsigned i = 0; i < DIM; i++) {
                    PetscVecTools::SetElement(initial_condition, DIM * petscvec_index + i, current_location(i));
                    PetscVecTools::SetElement(r_vector, DIM * petscvec_index + i, forces(i));
                }
            }
            if(allowRadiusVariation && !movementDisabled) {
                PetscVecTools::SetElement(initial_condition, DIM * cellCount + petscvec_index, radius);
                PetscVecTools::SetElement(r_vector, DIM * cellCount + petscvec_index, G);
            }else if(allowRadiusVariation && movementDisabled){
                PetscVecTools::SetElement(initial_condition, petscvec_index, radius);
                PetscVecTools::SetElement(r_vector, petscvec_index, G);
            }
        }
    }
    PetscMatTools::Finalise(r_matrix);
    solver.SetInitialConditionVector(initial_condition);


    // Solve to get solution at next timestep
    Vec soln_next_timestep = solver.SolveOneTimeStep();
    ReplicatableVector soln_next_timestep_repl(soln_next_timestep);


    // Iterate over all nodes associated with real cells and update their properties
    for (typename AbstractCellPopulation<DIM>::RealCellsIterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {

        if(cell_iter->GetCellData()->GetItem("IsBuskeKnot")==1){
            //Knot, ignore
        }else{

            // Get indices associated with this node
            unsigned global_node_index = this->GetLocationIndexUsingCell((*cell_iter));
            unsigned node_index = this->rGetMesh().SolveNodeMapping(global_node_index);
            c_vector<double, DIM> old_location = this->GetLocationOfCellCentre(*cell_iter);
            unsigned petscvec_index = node_index;
            if(petscvec_index >= cellCount){
                petscvec_index -= numKnots; 
            }

            if(!movementDisabled) {
                // Setup new node location
                c_vector<double, DIM> new_node_location;
                for (unsigned i = 0; i < DIM; i++) {
                    new_node_location(i) = soln_next_timestep_repl[DIM * petscvec_index + i];
                }

                double displacement = norm_2(new_node_location - old_location);
                if(displacement > this->mAbsoluteMovementThreshold){
                   throw (int)ceil(displacement);
                }

                // Move the node and set new radius as required.
                ChastePoint<DIM> new_point(new_node_location);
                this->SetNode(global_node_index, new_point);
                double currentRadius = cell_iter->GetCellData()->GetItem("Radius");
                this->GetNode(global_node_index)->SetRadius(currentRadius);
            }

            if(allowRadiusVariation && !movementDisabled) {
                double newRadius = soln_next_timestep_repl[DIM * cellCount + petscvec_index];
                this->GetCellUsingLocationIndex(global_node_index)->GetCellData()->SetItem("Radius", newRadius);
                this->GetNode(global_node_index)->SetRadius(newRadius);
            }else if(allowRadiusVariation && movementDisabled){
                double newRadius = soln_next_timestep_repl[petscvec_index];
                this->GetCellUsingLocationIndex(global_node_index)->GetCellData()->SetItem("Radius", newRadius);
                this->GetNode(global_node_index)->SetRadius(newRadius);
                //std::cout << "rad " << newRadius << std::endl;
            }

        }
    }

    // Tidy up
    PetscTools::Destroy(initial_condition);
}



template<unsigned DIM>
void NodeBasedCellPopulationWithBuskeUpdate<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<RadialVariation>" << allowRadiusVariation << "</RadialVariation>\n";
    *rParamsFile << "\t\t<AdhesionEnabled>" << adhesionEnabled << "</AdhesionEnabled>\n";
    *rParamsFile << "\t\t<Eps>" << Eps << "</Eps>\n";
    *rParamsFile << "\t\t<CompressionEnabled>" << compressionEnabled << "</CompressionEnabled>\n";
    *rParamsFile << "\t\t<K>" << K << "</K>\n";
    *rParamsFile << "\t\t<ElasticityEnabled>" << elasticityEnabled << "</ElasticityEnabled>\n";
    *rParamsFile << "\t\t<D>" << D << "</D>\n";
    *rParamsFile << "\t\t<IntercellularDrag>" << dampingConstantIntercell << "</IntercellularDrag>\n";
    *rParamsFile << "\t\t<VolumeChangeDrag>" << dampingConstantVolume << "</VolumeChangeDrag>\n";
    *rParamsFile << "\t\t<MediumDrag>" << dampingConstantMedium << "</MediumDrag>\n";


    // Call method on direct parent class
    NodeBasedCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class NodeBasedCellPopulationWithBuskeUpdate<1>;
template class NodeBasedCellPopulationWithBuskeUpdate<2>;
template class NodeBasedCellPopulationWithBuskeUpdate<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NodeBasedCellPopulationWithBuskeUpdate)
