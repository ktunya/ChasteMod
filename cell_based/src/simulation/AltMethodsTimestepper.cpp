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

#include "AltMethodsTimestepper.hpp"
#include "AbstractMesh.hpp"
#include "StepperChoice.hpp"
#include "StepSizeException.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NodeBasedCellPopulationWithParticles.hpp"
#include "CellBasedEventHandler.hpp"
#include "ReplicatableVector.hpp"
#include "PetscVecTools.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>	
AltMethodsTimestepper<ELEMENT_DIM,SPACE_DIM>::AltMethodsTimestepper( AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>&  inputCellPopulation, 
                                                 std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >&  inputForceCollection)
:rCellPopulation( inputCellPopulation ),
rForceCollection( inputForceCollection )
{	
	if(dynamic_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation)){
		nonEulerSteppersEnabled = true;
	}else{
		nonEulerSteppersEnabled = false;
	}

    if(dynamic_cast<MeshBasedCellPopulationWithGhostNodes<SPACE_DIM>*>(&rCellPopulation)){
        ghostNodeForcesEnabled = true;
    }else{
        ghostNodeForcesEnabled = false;
    }

    if(dynamic_cast<NodeBasedCellPopulationWithParticles<SPACE_DIM>*>(&rCellPopulation)){
        particleForcesEnabled = true;
    }else{
        particleForcesEnabled = false;
    }

	pNonlinearSolver = new SimplePetscNonlinearSolver();
	implicitStepSize = 0;
};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AltMethodsTimestepper<ELEMENT_DIM,SPACE_DIM>::~AltMethodsTimestepper(){	
};



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<c_vector<double, SPACE_DIM> > AltMethodsTimestepper<ELEMENT_DIM,SPACE_DIM>::ComputeAndSaveForces(){
	
	CellBasedEventHandler::BeginEvent(CellBasedEventHandler::FORCE);

    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
         node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter)
    {
        node_iter->ClearAppliedForce();
    }
    
    for (typename std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >::iterator iter = rForceCollection.begin();
        iter != rForceCollection.end(); ++iter)
    {
        (*iter)->AddForceContribution(rCellPopulation);
    }

    // Here deal with forces on non-cell nodes (ghosts and particles)
    if(ghostNodeForcesEnabled){
        dynamic_cast<MeshBasedCellPopulationWithGhostNodes<SPACE_DIM>*>(&rCellPopulation)->ApplyGhostForces();
    }
    if(particleForcesEnabled){
        dynamic_cast<NodeBasedCellPopulationWithParticles<SPACE_DIM>*>(&rCellPopulation)->ApplyParticleForces();
    }

    CellBasedEventHandler::EndEvent(CellBasedEventHandler::FORCE);

    // Store the forces in a vector
	std::vector<c_vector<double, SPACE_DIM> > forcesAsVector;
	forcesAsVector.reserve(rCellPopulation.GetNumNodes());

    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
         node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter)
    {
        forcesAsVector.push_back(node_iter->rGetAppliedForce()); 
    }

    return forcesAsVector;
};



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<c_vector<double, SPACE_DIM> > AltMethodsTimestepper<ELEMENT_DIM,SPACE_DIM>::SaveCurrentLocations(){

	std::vector<c_vector<double, SPACE_DIM> > currentLocations;
	currentLocations.reserve(rCellPopulation.GetNumNodes());

	int index = 0;
	for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
         node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd(); 
         ++index, ++node_iter)
    {
		currentLocations[index] = node_iter->rGetLocation();	
	}

	return currentLocations;
}



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AltMethodsTimestepper<ELEMENT_DIM,SPACE_DIM>::UpdateAllNodePositions(double dt, int stepper){

	if(nonEulerSteppersEnabled){

		switch( stepper ){

    		case StepperChoice::EULER :
    		{
    			std::vector<c_vector<double, SPACE_DIM> > F = ComputeAndSaveForces();

    			int index = 0;
    	        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
    	     		 node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    	        {
    	        	double damping = rCellPopulation.GetDampingConstant(node_iter->GetIndex());
                    
                    c_vector<double, SPACE_DIM> oldLocation = node_iter->rGetLocation();
                    c_vector<double, SPACE_DIM> newLocation = oldLocation + dt * F[index]/damping;
                    
                    rCellPopulation.CheckForStepSizeException(norm_2(newLocation-oldLocation), dt, node_iter->GetIndex());
    	        	
                    ChastePoint<SPACE_DIM> new_point(newLocation);
                    rCellPopulation.SetNode(node_iter->GetIndex(), new_point);
    	        }
    	        
    		}
    		break;

    		case StepperChoice::RK4 :
    		{
    			std::vector<c_vector<double, SPACE_DIM> > K1 = ComputeAndSaveForces();

    			int index = 0;
    			for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
    	     		 node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    	        {
    	        	double damping = rCellPopulation.GetDampingConstant(node_iter->GetIndex());
    	        	c_vector<double, SPACE_DIM> newLocation = node_iter->rGetLocation() + dt * K1[index]/(damping*2.0);
                    ChastePoint<SPACE_DIM> new_point(newLocation);
                    rCellPopulation.SetNode(node_iter->GetIndex(), new_point);
    	        }
    	        
    	        std::vector< c_vector<double, SPACE_DIM> > K2 = ComputeAndSaveForces(); 

    	        index = 0;
    	        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
    	     		 node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    	        {
    	        	double damping = rCellPopulation.GetDampingConstant(node_iter->GetIndex());
    	        	c_vector<double, SPACE_DIM> newLocation = node_iter->rGetLocation() + dt * (K2[index] - K1[index])/(damping*2.0); //revert, then update
    	            ChastePoint<SPACE_DIM> new_point(newLocation);
                    rCellPopulation.SetNode(node_iter->GetIndex(), new_point);
                }
    	                
    	        std::vector< c_vector<double, SPACE_DIM> > K3 = ComputeAndSaveForces(); 

    	        index = 0;
    	        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
    	     		 node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    	        {
    	        	double damping = rCellPopulation.GetDampingConstant(node_iter->GetIndex());
    	        	c_vector<double, SPACE_DIM> newLocation = node_iter->rGetLocation() + dt * (K3[index] - K2[index]/2.0) /damping; //revert, then update
                    ChastePoint<SPACE_DIM> new_point(newLocation);
                    rCellPopulation.SetNode(node_iter->GetIndex(), new_point);
                }
    	                
    	        std::vector< c_vector<double, SPACE_DIM> > K4 = ComputeAndSaveForces(); 

    	        index = 0;
    	        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
    	     		 node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    	        {
    	        	double damping = rCellPopulation.GetDampingConstant(node_iter->GetIndex());
    	        	c_vector<double, SPACE_DIM> effectiveForce = (K1[index] + 2*K2[index] + 2*K3[index] + K4[index] )/6.0;

                    c_vector<double, SPACE_DIM> oldLocation = node_iter->rGetLocation() - dt * K3[index]/damping;
    	        	c_vector<double, SPACE_DIM> newLocation = oldLocation + dt * (effectiveForce/damping); 
                    
                    rCellPopulation.CheckForStepSizeException(norm_2(newLocation-oldLocation), dt, node_iter->GetIndex());
                    ChastePoint<SPACE_DIM> new_point(newLocation);
                    rCellPopulation.SetNode(node_iter->GetIndex(), new_point);

                    //Ensure the nodes hold accurate forces, incase they're accessed by some other class
    	        	node_iter->ClearAppliedForce();
    	        	node_iter->AddAppliedForceContribution(effectiveForce);
    	        }
    	 	}
    		break;

    		case StepperChoice::BACKWARDEULER :
    		{
    	        std::vector< c_vector<double,SPACE_DIM> > initialLocations = SaveCurrentLocations();
    	        
    	        implicitStepSize = dt;
    	        unsigned systemSize = rCellPopulation.GetNumNodes() * SPACE_DIM;

    	        // Setup an initial condition consisting of the current node locations + one forward Euler step
    	        this->UpdateAllNodePositions(0.005, StepperChoice::EULER);
    	        Vec initialCondition = PetscTools::CreateAndSetVec(systemSize, 0.0);

    	        int index = 0;
    	        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
    	             node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    	        {
    	            double damping = rCellPopulation.GetDampingConstant(node_iter->GetIndex());
    	            
    	            c_vector<double, SPACE_DIM> location = node_iter->rGetLocation();
    	            for(int i=0; i<SPACE_DIM; i++){
    	                PetscVecTools::SetElement(initialCondition, SPACE_DIM*index + i, location[i]);
    	            }         
    	        }

    	        // Call nonlinear solver
    	        Vec solnNextTimestep = pNonlinearSolver->Solve( &BACKWARDEULER_ComputeResidual<ELEMENT_DIM, SPACE_DIM>,  
    	                                                         SNESComputeJacobianDefault,  
    	                                                         initialCondition,   
    	                                                         UINT_MAX,          
    	                                                         this);              
    	        // Unpack solution. 
    	        ReplicatableVector solnNextTimestepRepl(solnNextTimestep);

    	        index = 0;
    	        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
    	            node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    	        {
    	            double damping = rCellPopulation.GetDampingConstant(node_iter->GetIndex());

    	            c_vector<double, SPACE_DIM> oldLocation = initialLocations[index];
    	            c_vector<double, SPACE_DIM> newLocation;
    	            
    	            for(int i=0; i<SPACE_DIM; i++){
    	                newLocation[i] = solnNextTimestepRepl[SPACE_DIM * index + i];
    	            }
    	            
                    rCellPopulation.CheckForStepSizeException(norm_2(newLocation-oldLocation), dt, node_iter->GetIndex());
                    ChastePoint<SPACE_DIM> new_point(newLocation);
                    rCellPopulation.SetNode(node_iter->GetIndex(), new_point);
    	            
                    node_iter->ClearAppliedForce();
    	            c_vector<double, SPACE_DIM> effectiveForce = (damping/dt)*(newLocation-oldLocation);
    	            node_iter->AddAppliedForceContribution(effectiveForce);
    	        }

    	        PetscTools::Destroy(initialCondition);
    		}
    		break;
		}

	}else{

		ComputeAndSaveForces();
		rCellPopulation.UpdateNodeLocations(dt);		
	}
};	



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AltMethodsTimestepper<ELEMENT_DIM,SPACE_DIM>::BACKWARDEULERComputeResidual(const Vec currentGuess, Vec residualVector){

    std::vector< c_vector<double, SPACE_DIM> > currentLocations = SaveCurrentLocations();

    ReplicatableVector guessPositions(currentGuess);

    int index = 0;
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
        node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    {  
        c_vector<double, SPACE_DIM> guessLocation; 
        for(int i=0; i<SPACE_DIM; i++){
            guessLocation[i] = guessPositions[SPACE_DIM * index + i];
        }
        node_iter->rGetModifiableLocation() = guessLocation;
    }

    // Get force at the guess locations
    std::vector< c_vector<double, SPACE_DIM> > Fguess = ComputeAndSaveForces();
    
    index = 0;
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
         node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    {     
        node_iter->rGetModifiableLocation() = currentLocations[index];

        for(int i=0; i<SPACE_DIM; i++){

            double damping = rCellPopulation.GetDampingConstant(node_iter->GetIndex());

            double residual_ith_cpt = guessPositions[SPACE_DIM * index + i] - currentLocations[index][i] - implicitStepSize * Fguess[index][i]/damping;

            PetscVecTools::SetElement(residualVector, SPACE_DIM * index + i, residual_ith_cpt);
        }
    } 
};



///////// Explicit instantiation
template class AltMethodsTimestepper<1,1>;
template class AltMethodsTimestepper<1,2>;
template class AltMethodsTimestepper<2,2>;
template class AltMethodsTimestepper<1,3>;
template class AltMethodsTimestepper<2,3>;
template class AltMethodsTimestepper<3,3>;