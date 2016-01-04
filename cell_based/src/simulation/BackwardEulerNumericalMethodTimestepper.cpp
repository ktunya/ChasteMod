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

#include "BackwardEulerNumericalMethodTimestepper.hpp"
#include "ReplicatableVector.hpp"
#include "PetscVecTools.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>	
BackwardEulerNumericalMethodTimestepper<ELEMENT_DIM,SPACE_DIM> :: BackwardEulerNumericalMethodTimestepper( 
	                                               AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>&                  inputCellPopulation, 
                                                 std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >&  inputForceCollection)
:AbstractNumericalMethodTimestepper<ELEMENT_DIM,SPACE_DIM> ( inputCellPopulation, inputForceCollection)
{	
    pNonlinearSolver = new SimplePetscNonlinearSolver();
    pNonlinearSolver->SetTolerance(1e-5);
    implicitStepSize = 0;
};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>  
BackwardEulerNumericalMethodTimestepper<ELEMENT_DIM,SPACE_DIM>::~BackwardEulerNumericalMethodTimestepper(){
};



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>  
void BackwardEulerNumericalMethodTimestepper<ELEMENT_DIM,SPACE_DIM>::UpdateAllNodePositions(double dt){

    if(this->nonEulerSteppersEnabled){

      implicitStepSize = dt;
      unsigned systemSize = this->rCellPopulation.GetNumNodes() * SPACE_DIM;
      
      std::vector< c_vector<double,SPACE_DIM> > initialLocations = this->SaveCurrentLocations();
      std::vector< c_vector<double,SPACE_DIM> > initialF = this->ComputeAndSaveForces();                

      // Setup an initial condition consisting of the current node locations + one forward Euler step 
      Vec initialCondition = PetscTools::CreateAndSetVec(systemSize, 0.0);
      double FEStepSize = 0.01; 
      int index = 0;
      for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rCellPopulation.rGetMesh().GetNodeIteratorBegin();
           node_iter != this->rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
      {   
          c_vector<double, SPACE_DIM> location = node_iter->rGetLocation();
          for(int i=0; i<SPACE_DIM; i++){
              PetscVecTools::SetElement(initialCondition, SPACE_DIM*index + i,  location[i] + FEStepSize * initialF[index][i]);
          }         
      }

      // Call nonlinear solver
      Vec solnNextTimestep = pNonlinearSolver->Solve( &BACKWARDEULER_ComputeResidual<ELEMENT_DIM, SPACE_DIM>,  
                                                      &SNESComputeJacobianDefault,  
                                                      initialCondition,   
                                                      UINT_MAX,          
                                                      this);              
      // Unpack solution. 
      ReplicatableVector solnNextTimestepRepl(solnNextTimestep);

      index = 0;
      for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rCellPopulation.rGetMesh().GetNodeIteratorBegin();
          node_iter != this->rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
      {
          c_vector<double, SPACE_DIM> oldLocation = initialLocations[index];
          c_vector<double, SPACE_DIM> newLocation;      
          for(int i=0; i<SPACE_DIM; i++){
              newLocation[i] = solnNextTimestepRepl[SPACE_DIM * index + i];
          }
          
          c_vector<double, SPACE_DIM> displacement = this->rCellPopulation.rGetMesh().GetVectorFromAtoB(oldLocation, newLocation);
          this->HandleStepSizeExceptions(&displacement, dt, node_iter->GetIndex());

          ChastePoint<SPACE_DIM> new_point(newLocation);
          this->rCellPopulation.SetNode(node_iter->GetIndex(), new_point);
        
          node_iter->ClearAppliedForce();
          double damping = this->rCellPopulation.GetDampingConstant(node_iter->GetIndex());
          c_vector<double, SPACE_DIM> effectiveForce = (damping/dt)*displacement;
          node_iter->AddAppliedForceContribution(effectiveForce);
      }

      PetscTools::Destroy(initialCondition);

    }else{

        this->ComputeAndSaveForces();
        this->rCellPopulation.UpdateNodeLocations(dt);    

    }

};



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BackwardEulerNumericalMethodTimestepper<ELEMENT_DIM,SPACE_DIM>::BACKWARDEULERComputeResidual(const Vec currentGuess, Vec residualVector){

    std::vector< c_vector<double, SPACE_DIM> > currentLocations = this->SaveCurrentLocations();

    ReplicatableVector guessPositions(currentGuess);

    int index = 0;
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rCellPopulation.rGetMesh().GetNodeIteratorBegin();
        node_iter != this->rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    {  
        c_vector<double, SPACE_DIM> guessLocation; 
        for(int i=0; i<SPACE_DIM; i++){
            guessLocation[i] = guessPositions[SPACE_DIM * index + i];
        }
        node_iter->rGetModifiableLocation() = guessLocation;
    }

    // Get force at the guess locations
    std::vector< c_vector<double, SPACE_DIM> > Fguess = this->ComputeAndSaveForces();
    
    index = 0;
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rCellPopulation.rGetMesh().GetNodeIteratorBegin();
         node_iter != this->rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    {     
        node_iter->rGetModifiableLocation() = currentLocations[index];

        for(int i=0; i<SPACE_DIM; i++){

            double residual_ith_cpt = guessPositions[SPACE_DIM * index + i] - currentLocations[index][i] - implicitStepSize * Fguess[index][i];

            PetscVecTools::SetElement(residualVector, SPACE_DIM * index + i, residual_ith_cpt);
        }
    } 
};

///////// Explicit instantiation
template class BackwardEulerNumericalMethodTimestepper<1,1>;
template class BackwardEulerNumericalMethodTimestepper<1,2>;
template class BackwardEulerNumericalMethodTimestepper<2,2>;
template class BackwardEulerNumericalMethodTimestepper<1,3>;
template class BackwardEulerNumericalMethodTimestepper<2,3>;
template class BackwardEulerNumericalMethodTimestepper<3,3>;