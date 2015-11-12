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

#include "OffLatticeSimulation.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "T2SwapCellKiller.hpp"
#include "AbstractVertexBasedDivisionRule.hpp"
#include "Cylindrical2dMesh.hpp"
#include "Cylindrical2dVertexMesh.hpp"
#include "AbstractTwoBodyInteractionForce.hpp"
#include "CellBasedEventHandler.hpp"
#include "LogFile.hpp"
#include "Version.hpp"
#include "ExecutableSupport.hpp"
#include "ReplicatableVector.hpp"
#include "GeneralisedLinearSpringForce.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::OffLatticeSimulation(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                                                bool deleteCellPopulationInDestructor,
                                                bool initialiseCells,
                                                bool adaptiveChoice,
                                                int  stepperChoice)
    : AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>(rCellPopulation, deleteCellPopulationInDestructor, initialiseCells),
    adaptive(adaptiveChoice),
    stepper(stepperChoice)
{
    if (!dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation))
    {
        EXCEPTION("OffLatticeSimulations require a subclass of AbstractOffLatticeCellPopulation.");
    }

    // Different time steps are used for cell-centre and vertex-based simulations
    if (bool(dynamic_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation)))
    {
        this->mDt = 1.0/120.0; // 30 seconds
    }
    else if (bool(dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation)))
    {
        this->mDt = 0.002; // smaller time step required for convergence/stability

        // For VertexBasedCellPopulations we automatically add a T2SwapCellKiller. In order to inhibit T2 swaps
        // the user needs to set the threshold for T2 swaps in the mesh to 0.
        VertexBasedCellPopulation<SPACE_DIM>* p_vertex_based_cell_population = dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation);
        MAKE_PTR_ARGS(T2SwapCellKiller<SPACE_DIM>, T2_swap_cell_killer, (p_vertex_based_cell_population));
        this->AddCellKiller(T2_swap_cell_killer);
    }
    else
    {
        // All classes derived from AbstractOffLatticeCellPopulation are covered by the above (except user-derived classes),
        // i.e. if you want to use this class with your own subclass of AbstractOffLatticeCellPopulation, then simply
        // comment out the line below
        NEVER_REACHED;
    }

    // Setup for Adams Moulton if required
    pNewtonSolver = NULL;
    currentAMStep = 0;
    if(stepperChoice == StepperChoice::ADAMSM){
        pNewtonSolver = new SimpleNewtonNonlinearSolver(); 
        pNewtonSolver->SetTolerance(1e-6);
        pNewtonSolver->SetWriteStats(true);
        currentAMStep = this->mDt;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::AddForce(boost::shared_ptr<AbstractForce<ELEMENT_DIM,SPACE_DIM> > pForce)
{
    mForceCollection.push_back(pForce);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::RemoveAllForces()
{
    mForceCollection.clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::AddCellPopulationBoundaryCondition(boost::shared_ptr<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM> > pBoundaryCondition)
{
    mBoundaryConditions.push_back(pBoundaryCondition);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::RemoveAllCellPopulationBoundaryConditions()
{
    mBoundaryConditions.clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::UpdateCellLocationsAndTopology()
{

    double timeAdvanced = 0;
    double targetDt = this->mDt;
    double currentStepSize = targetDt;

    double safetyFactor = 0.95;
    double movementThreshold = dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM, SPACE_DIM>* >(&(this->mrCellPopulation))
                               ->GetAbsoluteMovementThreshold();

    while(timeAdvanced < targetDt){
            
        // Try to update positions
        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::POSITION);

        // Store the initial node positions (these may be needed when applying boundary conditions)    
        std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> > old_node_locations;
        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
            node_iter != this->mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
            ++node_iter)
        {
            old_node_locations[&(*node_iter)] = (node_iter)->rGetLocation();
        }

        try{

            switch( stepper ){
                
                case StepperChoice::EULER :
                {
                    std::vector< c_vector<double, SPACE_DIM> > F = ApplyForces();
                    static_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&(this->mrCellPopulation))
                                                                        ->UpdateNodeLocations(currentStepSize);
                }
                break;
                
                case StepperChoice::RK :
                {
                    std::vector< c_vector<double, SPACE_DIM> > K1 = ApplyForces(); 
                    //K1 is now stored in the population forces
                    static_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&(this->mrCellPopulation))
                                                                        ->UpdateNodeLocations(currentStepSize/2.0);
                
                    std::vector< c_vector<double, SPACE_DIM> > K2 = ApplyForces(); 
                    //K2 is now stored in the population forces                                                    
                    RevertToOldLocations(old_node_locations);
                    static_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&(this->mrCellPopulation))
                                                                        ->UpdateNodeLocations(currentStepSize/2.0);
                    
                    std::vector< c_vector<double, SPACE_DIM> > K3 = ApplyForces(); 
                    //K3 is now stored in the population forces
                    RevertToOldLocations(old_node_locations);
                    static_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&(this->mrCellPopulation))
                                                                        ->UpdateNodeLocations(currentStepSize);
                    
                    std::vector< c_vector<double, SPACE_DIM> > K4 = ApplyForces(); 
                    //K4 is now stored in the population forces
                    //Reapply the other K contributions...
                    RevertToOldLocations(old_node_locations);
                    AddForceVecWithMultiplyingFactor(K1,1);
                    AddForceVecWithMultiplyingFactor(K2,2);
                    AddForceVecWithMultiplyingFactor(K3,2);
                    static_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&(this->mrCellPopulation))
                                                                        ->UpdateNodeLocations(currentStepSize/6.0);                                                        
                }           
                break;

                case StepperChoice::ADAMSM:
                {
                    currentAMStep = currentStepSize;

                    int cellCount = this->mrCellPopulation.GetNumNodes();
                    unsigned systemSize = cellCount * SPACE_DIM;

                    // Setup initial condition
                    Vec initialCondition = PetscTools::CreateAndSetVec(systemSize, 0.0);

                    int tempIteratorIndex = 0;
                    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
                        node_iter != this->mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
                        ++node_iter)
                    {
                        c_vector<double, SPACE_DIM> location = node_iter->rGetLocation();
                        for(int i=0; i<SPACE_DIM; i++){
                            PetscVecTools::SetElement(initialCondition, SPACE_DIM*tempIteratorIndex+i, location[i]);
                        }
                        tempIteratorIndex++;
                    }

                    Vec solnNextTimestep = pNewtonSolver->Solve( &OffLatticeSimulation_AdamsM_ComputeResidual<ELEMENT_DIM, SPACE_DIM>,  /*Compute residual fn ptr*/
                                                                 &OffLatticeSimulation_AdamsM_ComputeJacobian<ELEMENT_DIM, SPACE_DIM>, /*Compute Jacobian fn ptr*/ 
                                                                 initialCondition,   /*Current positions*/
                                                                 UINT_MAX,           /*Triggers automatic preallocation*/
                                                                 this);              /*Optional pContext containing this class*/

                    ReplicatableVector solnNextTimestepRepl(solnNextTimestep);

                    //Copying out results
                    tempIteratorIndex = 0;
                    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
                        node_iter != this->mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
                        ++node_iter)
                    {
                        c_vector<double, SPACE_DIM> oldLocation = node_iter->rGetLocation();
                
                        c_vector<double, SPACE_DIM> newLocation;
                        double displacement = 0;
                        for(int i=0; i<SPACE_DIM; i++){
                            newLocation[i] = solnNextTimestepRepl[SPACE_DIM*tempIteratorIndex+i];
                            displacement += (oldLocation[i] - newLocation[i]) * (oldLocation[i] - newLocation[i]);
                        }
                        
                        displacement = pow(displacement, 0.5);
                        if(displacement > movementThreshold){
                            throw (int)ceil(displacement);
                        }

                        node_iter->rGetModifiableLocation() = newLocation;
                        tempIteratorIndex++;
                    }

                    // Cleanup
                    PetscTools::Destroy(initialCondition);
                    
                }
                break;
            }

            ApplyBoundaries(old_node_locations);

            // Successful timestep. Update totalTimeAdvanced and increase the timestep by 1% if possible 
            timeAdvanced += currentStepSize;
            if(adaptive){
                currentStepSize = fmin(1.01*currentStepSize, targetDt-timeAdvanced);
            }

        }catch(int e){

            if(adaptive){

                // Update failed
                // Revert nodes to their old positions.
                RevertToOldLocations(old_node_locations);
                // Estimate the largest permissible / remaining timestep. 
                currentStepSize = fmin(safetyFactor * currentStepSize * (movementThreshold/(double)e), targetDt-timeAdvanced); 

            }else{
                // Adaptive set to false. Crash with the usual error.
                EXCEPTION("Cells are moving by approximately: " << e <<
                ", which is more than the AbsoluteMovementThreshold. Use a smaller timestep to avoid this exception.");
            }
        }

        CellBasedEventHandler::EndEvent(CellBasedEventHandler::POSITION);
    }

}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::RevertToOldLocations(std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> > old_node_locations){
    
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
        node_iter != this->mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
        ++node_iter)
    {
        (node_iter)->rGetModifiableLocation() = old_node_locations[&(*node_iter)];
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector< c_vector<double, SPACE_DIM> > OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::ApplyForces(){

    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::FORCE);

    // Clear all existing forces
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
         node_iter != this->mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        node_iter->ClearAppliedForce();
    }
    
    // Now add force contributions from each AbstractForce.  
    for (typename std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >::iterator iter = mForceCollection.begin();
         iter != mForceCollection.end();
         ++iter)
    {
        (*iter)->AddForceContribution(this->mrCellPopulation);
    }

    CellBasedEventHandler::EndEvent(CellBasedEventHandler::FORCE);

    // Store a force vector in case the stepper requires it.  
    std::vector< c_vector<double,SPACE_DIM> > forceVector;
    forceVector.reserve(this->mrCellPopulation.GetNumNodes());

    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
         node_iter != this->mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        forceVector.push_back(node_iter->rGetAppliedForce()); 
    }

    return forceVector;
};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::AddForceVecWithMultiplyingFactor(std::vector< c_vector<double, SPACE_DIM> > fVec, int factor){

    int tempIteratorIndex = 0;
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
        node_iter != this->mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
        ++node_iter)
    {
        c_vector<double, SPACE_DIM> force = factor * fVec[tempIteratorIndex];
        (node_iter)->AddAppliedForceContribution( force );
        tempIteratorIndex++;
    }
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::ApplyBoundaries(std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> > old_node_locations)
{

    // Apply any boundary conditions
    for (typename std::vector<boost::shared_ptr<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM> > >::iterator bcs_iter = mBoundaryConditions.begin();
         bcs_iter != mBoundaryConditions.end();
         ++bcs_iter)
    {
        (*bcs_iter)->ImposeBoundaryCondition(old_node_locations);
    }

    // Verify that each boundary condition is now satisfied
    for (typename std::vector<boost::shared_ptr<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM> > >::iterator bcs_iter = mBoundaryConditions.begin();
         bcs_iter != mBoundaryConditions.end();
         ++bcs_iter)
    {
        if (!((*bcs_iter)->VerifyBoundaryCondition()))
        {
            EXCEPTION("The cell population boundary conditions are incompatible.");
        }
    }
}




template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::CalculateCellDivisionVector(CellPtr pParentCell)
{
    /**
     * \todo #2400 Could remove this dynamic_cast by moving the code block below into
     * AbstractCentreBasedCellPopulation::AddCell(), allowing it to be overruled by
     * this method when overridden in subclasses. See also comment on #1093.
     */
    // If it is not vertex based...
    if (bool(dynamic_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&(this->mrCellPopulation))))
    {
        // Location of parent and daughter cells
        c_vector<double, SPACE_DIM> parent_coords = this->mrCellPopulation.GetLocationOfCellCentre(pParentCell);
        c_vector<double, SPACE_DIM> daughter_coords;

        // Get separation parameter
        double separation = static_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&(this->mrCellPopulation))->GetMeinekeDivisionSeparation();

        // Make a random direction vector of the required length
        c_vector<double, SPACE_DIM> random_vector;

        /*
         * Pick a random direction and move the parent cell backwards by 0.5*separation
         * in that direction and return the position of the daughter cell 0.5*separation
         * forwards in that direction.
         */
        switch (SPACE_DIM)
        {
            case 1:
            {
                double random_direction = -1.0 + 2.0*(RandomNumberGenerator::Instance()->ranf() < 0.5);

                random_vector(0) = 0.5*separation*random_direction;
                break;
            }
            case 2:
            {
                double random_angle = 2.0*M_PI*RandomNumberGenerator::Instance()->ranf();

                random_vector(0) = 0.5*separation*cos(random_angle);
                random_vector(1) = 0.5*separation*sin(random_angle);
                break;
            }
            case 3:
            {
                /*
                 * Note that to pick a random point on the surface of a sphere, it is incorrect
                 * to select spherical coordinates from uniform distributions on [0, 2*pi) and
                 * [0, pi) respectively, since points picked in this way will be 'bunched' near
                 * the poles. See #2230.
                 */
                double u = RandomNumberGenerator::Instance()->ranf();
                double v = RandomNumberGenerator::Instance()->ranf();

                double random_azimuth_angle = 2*M_PI*u;
                double random_zenith_angle = std::acos(2*v - 1);

                random_vector(0) = 0.5*separation*cos(random_azimuth_angle)*sin(random_zenith_angle);
                random_vector(1) = 0.5*separation*sin(random_azimuth_angle)*sin(random_zenith_angle);
                random_vector(2) = 0.5*separation*cos(random_zenith_angle);
                break;
            }
            default:
                // This can't happen
                NEVER_REACHED;
        }

        parent_coords = parent_coords - random_vector;
        daughter_coords = parent_coords + random_vector;

        // Set the parent to use this location
        ChastePoint<SPACE_DIM> parent_coords_point(parent_coords);
        unsigned node_index = this->mrCellPopulation.GetLocationIndexUsingCell(pParentCell);
        this->mrCellPopulation.SetNode(node_index, parent_coords_point);

        return daughter_coords;
    }
    else
    {
        // Check this is a Vertex based cell population (in case new types are added later!).
        assert(bool(dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(&(this->mrCellPopulation))));

        VertexBasedCellPopulation<SPACE_DIM>* p_vertex_population = dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(&(this->mrCellPopulation));
        boost::shared_ptr<AbstractVertexBasedDivisionRule<SPACE_DIM> > p_division_rule = p_vertex_population->GetVertexBasedDivisionRule();

        return p_division_rule->CalculateCellDivisionVector(pParentCell, *p_vertex_population);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::WriteVisualizerSetupFile()
{
    if (PetscTools::AmMaster())
    {
        for (unsigned i=0; i<this->mForceCollection.size(); i++)
        {
            // This may cause compilation problems, probably due to AbstractTwoBodyInteractionForce not having two template parameters
            ///\todo Check whether this comment is still valid

            boost::shared_ptr<AbstractForce<ELEMENT_DIM,SPACE_DIM> > p_force = this->mForceCollection[i];
            if (boost::dynamic_pointer_cast<AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM> >(p_force))
            {
                double cutoff = (boost::static_pointer_cast<AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM> >(p_force))->GetCutOffLength();
                *(this->mpVizSetupFile) << "Cutoff\t" << cutoff << "\n";
            }
        }

        // This is a quick and dirty check to see if the mesh is periodic
        if (bool(dynamic_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&this->mrCellPopulation)))
        {
           if (bool(dynamic_cast<Cylindrical2dMesh*>(&(dynamic_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&(this->mrCellPopulation))->rGetMesh()))))
           {
               *this->mpVizSetupFile << "MeshWidth\t" << this->mrCellPopulation.GetWidth(0) << "\n";
           }
        }
        else if (bool(dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(&this->mrCellPopulation)))
        {
           if (bool(dynamic_cast<Cylindrical2dVertexMesh*>(&(dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(&(this->mrCellPopulation))->rGetMesh()))))
           {
               *this->mpVizSetupFile << "MeshWidth\t" << this->mrCellPopulation.GetWidth(0) << "\n";
           }
        }
    }
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::SetupSolve()
{
    // Clear all forces
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
         node_iter != this->mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        node_iter->ClearAppliedForce();
    }
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::OutputAdditionalSimulationSetup(out_stream& rParamsFile)
{
    // Loop over forces
    *rParamsFile << "\n\t<Forces>\n";
    for (typename std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM,SPACE_DIM> > >::iterator iter = mForceCollection.begin();
         iter != mForceCollection.end();
         ++iter)
    {
        // Output force details
        (*iter)->OutputForceInfo(rParamsFile);
    }
    *rParamsFile << "\t</Forces>\n";

    // Loop over cell population boundary conditions
    *rParamsFile << "\n\t<CellPopulationBoundaryConditions>\n";
    for (typename std::vector<boost::shared_ptr<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM> > >::iterator iter = mBoundaryConditions.begin();
         iter != mBoundaryConditions.end();
         ++iter)
    {
        // Output cell Boundary condition details
        (*iter)->OutputCellPopulationBoundaryConditionInfo(rParamsFile);
    }
    *rParamsFile << "\t</CellPopulationBoundaryConditions>\n";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::OutputSimulationParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::OutputSimulationParameters(rParamsFile);
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const bool& OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::GetAdaptive() const{
    return adaptive;
};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const int& OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::GetStepper() const{
    return stepper;
};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void  OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::ComputeResidualAdamsM(const Vec currentGuess, Vec residualVector){

    ReplicatableVector posGuess(currentGuess);
    std::vector< c_vector<double, SPACE_DIM> > currentLocations;
    
    std::vector< c_vector<double, SPACE_DIM> > Fcurr = ApplyForces();

    // Force a position update to the guessed locations. Store the old locations
    int tempIteratorIndex = 0;
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
        node_iter != this->mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
        ++node_iter)
    {  
        currentLocations.push_back(node_iter->rGetLocation());

        c_vector<double, SPACE_DIM> guessLocation; 
        for(int i=0; i<SPACE_DIM; i++){
            guessLocation[i] = posGuess[SPACE_DIM*tempIteratorIndex + i];
        }
        node_iter->rGetModifiableLocation() = guessLocation;

        tempIteratorIndex++;
    }
    std::vector< c_vector<double, SPACE_DIM> > Fguess = ApplyForces();
    

    tempIteratorIndex = 0;
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
         node_iter != this->mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {     
        c_vector<double, SPACE_DIM> posCurr = currentLocations[tempIteratorIndex];
        node_iter->rGetModifiableLocation() = posCurr;

        for(int i=0; i<SPACE_DIM; i++){

            double residual_ith_cpt = posGuess[SPACE_DIM*tempIteratorIndex+i] 
                                      - posCurr[i] 
                                      - 0.5 * currentAMStep * (Fguess[tempIteratorIndex][i] + Fcurr[tempIteratorIndex][i]);

            PetscVecTools::SetElement(residualVector, SPACE_DIM*tempIteratorIndex+i, residual_ith_cpt);
        }
        tempIteratorIndex++;
    }  
};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::ComputeJacobianNumericallyAdamsM(const Vec currentGuess, Mat* pJacobian){

    unsigned num_unknowns = SPACE_DIM * this->mrCellPopulation.GetNumNodes();

    // Set up working vectors
    Vec perturbed_residual_down = PetscTools::CreateVec(num_unknowns);
    Vec perturbed_residual_up   = PetscTools::CreateVec(num_unknowns);
    Vec result                  = PetscTools::CreateVec(num_unknowns);

    // Copy the currentGuess vector; we will perturb the copy
    Vec current_guess_copy_up;
    PETSCEXCEPT( VecDuplicate(currentGuess, &current_guess_copy_up) );
    PETSCEXCEPT( VecCopy(currentGuess, current_guess_copy_up) );
    Vec current_guess_copy_down;
    PETSCEXCEPT( VecDuplicate(currentGuess, &current_guess_copy_down) );
    PETSCEXCEPT( VecCopy(currentGuess, current_guess_copy_down) );

    // Amount to perturb each input element by
    double h = 1e-6;
    double near_hsquared = 1e-12;

    //\todo: is this a sufficient ownership check? Do I need to run this check on current_guess_copy_down?
    PetscInt ilo, ihi;
    VecGetOwnershipRange(current_guess_copy_up, &ilo, &ihi); 
    unsigned lo = ilo;
    unsigned hi = ihi;

    // Iterate over entries in the input vector
    for (unsigned global_index_outer = 0; global_index_outer < num_unknowns; global_index_outer++)
    {
        // Only perturb if we own it. Calculate residuals
        PetscVecTools::AddToElement(current_guess_copy_up, global_index_outer, h);
        ComputeResidualAdamsM(current_guess_copy_up, perturbed_residual_up);
        PetscVecTools::AddToElement(current_guess_copy_down, global_index_outer, -h);
        ComputeResidualAdamsM(current_guess_copy_down, perturbed_residual_down);

        // result = (perturbed_residual_up - perturbed_residual_down) / 2*h
        PetscVecTools::WAXPY(result, -1.0, perturbed_residual_down, perturbed_residual_up);
        PetscVecTools::Scale(result, 1.0/(2*h));

        double* p_result;
        ///\todo This loop is setting the column "global_index_outer" of
        /// the Jacobian matrix to the result vector in a non-sparse way.
        PETSCEXCEPT( VecGetArray(result, &p_result) );

        for (unsigned global_index = lo; global_index < hi; global_index++)
        {
            double result_entry = p_result[global_index - lo];

            if (!CompareDoubles::IsNearZero(result_entry, near_hsquared))
            {
                PetscMatTools::SetElement(*pJacobian, global_index, global_index_outer, result_entry);
            }
        }
        PETSCEXCEPT( VecRestoreArray(result, &p_result) );

        PetscVecTools::AddToElement(current_guess_copy_up,   global_index_outer, -h);
        PetscVecTools::AddToElement(current_guess_copy_down, global_index_outer,  h);

        //std::cout << "Max Jacobian entry: " << maxJacobianEntry << std::endl;
    }

    PetscTools::Destroy(perturbed_residual_up);
    PetscTools::Destroy(perturbed_residual_down);
    PetscTools::Destroy(current_guess_copy_up);
    PetscTools::Destroy(current_guess_copy_down);
    PetscTools::Destroy(result);
    PetscMatTools::Finalise(*pJacobian);
};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::ComputeDefaultJacobianGenLinearSpringForceAdamsM(const Vec currentGuess, Mat* pJacobian){

    unsigned num_unknowns = SPACE_DIM * this->mrCellPopulation.GetNumNodes();
    unsigned num_nodes = this->mrCellPopulation.GetNumNodes();
    ReplicatableVector guessPositions(currentGuess);
    double negligableCpt = 1e-10;

    // Check that the force collection contains only a single GeneralizedLinearSpringForce, and extract its properties.-------------------------
    if(mForceCollection.size()>1){
        EXCEPTION("Adams Moulton method currently only supports a single GeneralizedLinearSpringForce.");
    }
    if(!bool(dynamic_cast< GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>* >(mForceCollection[0].get()) )){
        EXCEPTION("Adams Moulton method currently only has an analytic Jacobian for the GeneralizedLinearSpringForce.");
    }

    GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>* forceLawPtr = 
                                                dynamic_cast< GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>* >(mForceCollection[0].get());

    bool useCutoff = forceLawPtr->GetUseCutOffLength();
    double cutoff = forceLawPtr->GetCutOffLength();
    double springConst = forceLawPtr->GetMeinekeSpringStiffness();
    //------------------------------------------------------------------------------------------------------------------------------------------


    // Loop over all derivative indices, calculating the Jacobian contribution from each.-------------------------------------------------------
    // N is the node whose position we are differentiating with respect to,
    // w is the component of that position under consideration 
    for (unsigned deriv_index = 0; deriv_index < num_unknowns; deriv_index++) 
    {
        int derivCellN = (int)(deriv_index / SPACE_DIM);
        int derivCptW = (int)(deriv_index % SPACE_DIM); 

        // Loop over cell pairs
        int AIndex = 0;
        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
             node_iter != this->mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
             ++node_iter) 
        {   

            // Get properties of cell A
            c_vector<double, SPACE_DIM> ALoc;
            for(int i=0; i<SPACE_DIM; i++){
                ALoc[i] = guessPositions[AIndex*SPACE_DIM + i];
            }
            double ARad = node_iter->GetRadius();
            std::cout << "ARad " << ARad << std::endl;
            double dFA_dN_cptW = 0;

            int BIndex = 0;
            for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator neighbour_iter = this->mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
                 neighbour_iter != this->mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
                 ++neighbour_iter)
            {

                // Get properties of cell B
                c_vector<double, SPACE_DIM> BLoc;
                for(int i=0; i<SPACE_DIM; i++){
                    BLoc[i] = guessPositions[BIndex*SPACE_DIM + i];
                }
                double BRad = neighbour_iter->GetRadius();
                std::cout << "BRad " << BRad << std::endl;

                // Work out whether A or B is the same node as N, in which case there is a contribution to the derivative
                bool N_equals_A = false;
                bool N_equals_B = false;
                if(AIndex == derivCellN){
                    N_equals_A = true;
                }
                if(BIndex == derivCellN){
                    N_equals_B = true;
                }

                if(N_equals_A || N_equals_B){
                    
                    // Work out whether this neighbour is close enough to actually make a contribution to the force on A,
                    // given the currentGuess positions
                    c_vector<double, SPACE_DIM> AtoB = BLoc - ALoc;
                    double separation = norm_2(AtoB); 
                    c_vector<double, SPACE_DIM> AtoB_unit = AtoB/separation;
                    double overlap = separation - ARad - BRad;

                    std::cout << "Sep " << separation << " Overlap " << overlap << std::endl;


                    bool makesContribution = true;
                    if(useCutoff){
                        if(separation > cutoff){
                            makesContribution = false;
                        }
                    }

                    // Determine Jacobian contribution 
                    if(makesContribution){

                        double dOverlap;
                        if(N_equals_A){
                            dOverlap = (ALoc[derivCptW]-BLoc[derivCptW])/separation;   
                        }
                        if(N_equals_B){
                            dOverlap = -(ALoc[derivCptW]-BLoc[derivCptW])/separation; 
                        }

                        double dAtoB_unit_dN_cptW; 
                        if(N_equals_A){
                            dAtoB_unit_dN_cptW = ( -overlap -dOverlap*(BLoc[derivCptW]-ALoc[derivCptW]) ) / (overlap*overlap);   
                        }
                        if(N_equals_B){
                            dAtoB_unit_dN_cptW = ( overlap  -dOverlap*(BLoc[derivCptW]-ALoc[derivCptW]) ) / (overlap*overlap);  
                        }

                        std::cout << "dOverlap " << dOverlap << " dAtoB_unit_dN_cptW " << dAtoB_unit_dN_cptW << std::endl;
                        if(overlap < 0){
                            //Log part of the force law applies
                            dFA_dN_cptW += springConst*(ARad+BRad)* ( dAtoB_unit_dN_cptW * log( 1+(overlap/(ARad+BRad)) ) +
                                                                      AtoB_unit[derivCptW] * dOverlap * ( 1/(ARad+BRad+overlap) ) );
                        }else{
                            //Exponential part of the force law applies
                            dFA_dN_cptW += springConst*(ARad+BRad)* ( dAtoB_unit_dN_cptW * overlap * exp(-5*overlap/(ARad+BRad)) + 
                                                                      AtoB_unit[derivCptW] * ( dOverlap * exp(-5*overlap/(ARad+BRad)) + 
                                                                      (-5*dOverlap/(ARad+BRad)) * overlap * exp(-5*overlap/(ARad+BRad)) ));
                        }
                    }
                }

                BIndex++;
            }

            //Add to Jacobian.
            if (!CompareDoubles::IsNearZero(dFA_dN_cptW, negligableCpt))
            {
                PetscMatTools::SetElement(*pJacobian, AIndex, deriv_index, dFA_dN_cptW);
            }

            AIndex++;
        }    
    }

    PetscMatTools::Finalise(*pJacobian); 

    for (int i=0; i<num_unknowns; i++) {
      for (int j=0; j<num_unknowns; j++) {
        std::cout << "Jcpt " << i << " " << j << ": " << PetscMatTools::GetElement(*pJacobian,i,j) << std::endl;
      }
    };
   
};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PetscErrorCode OffLatticeSimulation_AdamsM_ComputeResidual(SNES snes, Vec currentGuess, Vec residualVector, void* pContext)
{
    //Extract the simulation from the context:
    OffLatticeSimulation<ELEMENT_DIM, SPACE_DIM>* pSim = (OffLatticeSimulation<ELEMENT_DIM, SPACE_DIM>*)pContext; 
    pSim->ComputeResidualAdamsM(currentGuess, residualVector);
    return 0;
};

#if ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=5 )

    template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
    PetscErrorCode OffLatticeSimulation_AdamsM_ComputeJacobianWithComparison(SNES snes, Vec currentGuess, Mat globalJacobian, Mat preconditioner, void* pContext)
    {
        //Extract the simulation from the context:
        OffLatticeSimulation<ELEMENT_DIM, SPACE_DIM>* pSim = (OffLatticeSimulation<ELEMENT_DIM, SPACE_DIM>*)pContext;
        pSim->ComputeDefaultJacobianGenLinearSpringForceAdamsM(currentGuess, &globalJacobian);

        //PetscInt m;
        //PetscInt n;
        //MatGetLocalSize(globalJacobian, &m, &n);
        //Mat numericalJacobian;  
        //MatCreate(PETSC_COMM_WORLD, &numericalJacobian);
        //MatSetSizes(numericalJacobian, m, n, m, n);   
        //pSim->ComputeJacobianNumericallyAdamsM(currentGuess, &numericalJacobian);

        //bool checkResult = PetscMatTools::CheckEquality(globalJacobian, numericalJacobian, 1e-10);
        //if(!checkResult){
        //    std::cout << "Analytic Jacobian does not match numerical Jacobian" << std::endl;
        //    EXCEPTION("Analytic Jacobian does not match numerical Jacobian");
        //}

        return 0;
    };

    template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
    PetscErrorCode OffLatticeSimulation_AdamsM_ComputeJacobian(SNES snes, Vec currentGuess, Mat globalJacobian, Mat preconditioner, void* pContext)
    {
        //Extract the simulation from the context:
        OffLatticeSimulation<ELEMENT_DIM, SPACE_DIM>* pSim = (OffLatticeSimulation<ELEMENT_DIM, SPACE_DIM>*)pContext;
        pSim->ComputeJacobianNumericallyAdamsM(currentGuess, &globalJacobian);
        return 0;
    };

#else

    template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
    PetscErrorCode OffLatticeSimulation_AdamsM_ComputeJacobianWithComparison(SNES snes, Vec currentGuess, Mat* pGlobalJacobian, Mat* pPreconditioner, MatStructure* pMatStructure, void* pContext)
    {
        //Extract the simulation from the context:
        OffLatticeSimulation<ELEMENT_DIM, SPACE_DIM>* pSim = (OffLatticeSimulation<ELEMENT_DIM, SPACE_DIM>*)pContext;
        pSim->ComputeDefaultJacobianGenLinearSpringForceAdamsM(currentGuess, pGlobalJacobian);

        //PetscInt m;
        //PetscInt n;
        //MatGetLocalSize(*pGlobalJacobian, &m, &n);
        //Mat numericalJacobian;  
        //MatCreate(PETSC_COMM_WORLD, &numericalJacobian);
        //MatSetSizes(numericalJacobian, m, n, m, n);
        //pSim->ComputeJacobianNumericallyAdamsM(currentGuess, &numericalJacobian);   
        //
        //bool checkResult = PetscMatTools::CheckEquality(*pGlobalJacobian, numericalJacobian, 1e-10);
        //if(!checkResult){
        //    std::cout << "Analytic Jacobian does not match numerical Jacobian" << std::endl;
        //    EXCEPTION("Analytic Jacobian does not match numerical Jacobian");
        //}

        return 0;
    };

    template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
    PetscErrorCode OffLatticeSimulation_AdamsM_ComputeJacobian(SNES snes, Vec currentGuess, Mat* pGlobalJacobian, Mat* pPreconditioner, MatStructure* pMatStructure, void* pContext)
    {
        //Extract the simulation from the context:
        OffLatticeSimulation<ELEMENT_DIM, SPACE_DIM>* pSim = (OffLatticeSimulation<ELEMENT_DIM, SPACE_DIM>*)pContext;
        pSim->ComputeJacobianNumericallyAdamsM(currentGuess, pGlobalJacobian);
        return 0;
    };

#endif


///////// Explicit instantiation
template class OffLatticeSimulation<1,1>;
template class OffLatticeSimulation<1,2>;
template class OffLatticeSimulation<2,2>;
template class OffLatticeSimulation<1,3>;
template class OffLatticeSimulation<2,3>;
template class OffLatticeSimulation<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(OffLatticeSimulation)