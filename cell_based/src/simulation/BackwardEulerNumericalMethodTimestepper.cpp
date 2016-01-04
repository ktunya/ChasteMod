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
#include "PetscMatTools.hpp"
#include "GeneralisedLinearSpringForce.hpp"


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




template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void BackwardEulerNumericalMethodTimestepper<ELEMENT_DIM,SPACE_DIM>::BACKWARDEULERComputeSpringJacobian(const Vec input, Mat* pJacobian){

  unsigned num_unknowns = SPACE_DIM * this->rCellPopulation.GetNumNodes();
  unsigned num_nodes = this->rCellPopulation.GetNumNodes();
  
  ReplicatableVector inputPositions(input);

  // Check that the force collection contains only a single GeneralizedLinearSpringForce, and extract its properties.-------------------------
  if(this->rForceCollection.size() != 1){
      EXCEPTION("Backward Euler method currently only supports a single GeneralizedLinearSpringForce.");
  }
  if(!bool(dynamic_cast< GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>* >( this->rForceCollection[0].get()) )){
      EXCEPTION("Backward Euler method currently only has an analytic Jacobian for the GeneralizedLinearSpringForce.");
  }

  GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>* forceLawPtr = 
          dynamic_cast< GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>* >(this->rForceCollection[0].get());

  bool   useCutoff   = forceLawPtr->GetUseCutOffLength();
  double cutoff      = forceLawPtr->GetCutOffLength();
  double springConst = forceLawPtr->GetMeinekeSpringStiffness();
  //------------------------------------------------------------------------------------------------------------------------------------------

  // Loop over rows of the jacobian, calculating the value in each column---------------------------------------------------------------------
  // N is the node whose position we are differentiating with respect to,
  // w is the component of that position under consideration 

  for (unsigned deriv_index = 0; deriv_index < num_unknowns; deriv_index++) 
  {
      int derivCellN = (int)(deriv_index / SPACE_DIM);
      int derivCptW = (int)(deriv_index % SPACE_DIM); 

      // Loop over cells
      int NodeAIndex = 0;
      for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rCellPopulation.rGetMesh().GetNodeIteratorBegin();
           node_iter != this->rCellPopulation.rGetMesh().GetNodeIteratorEnd();
           ++node_iter, ++NodeAIndex) 
      {   

          // Get properties of cell A
          c_vector<double, SPACE_DIM> ALoc;
          for(int i=0; i<SPACE_DIM; i++){
              ALoc[i] = inputPositions[ NodeAIndex*SPACE_DIM + i ];
          }
          double ARad = node_iter->GetRadius();
          
          //Loop over potential neighbours and calculate contribution
          c_vector<double,SPACE_DIM> dFA_dNw;
          dFA_dNw[0] = 0;
          dFA_dNw[1] = 0;
          dFA_dNw[2] = 0;

          int NodeBIndex = 0;
          for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator neighbour_iter = this->rCellPopulation.rGetMesh().GetNodeIteratorBegin();
               neighbour_iter != this->rCellPopulation.rGetMesh().GetNodeIteratorEnd();
               ++neighbour_iter, ++NodeBIndex)
          {

              //Don't compare cells with themselves
              if(neighbour_iter->GetIndex() != node_iter->GetIndex()){

                  // Get properties of cell B
                  c_vector<double, SPACE_DIM> BLoc;
                  for(int i=0; i<SPACE_DIM; i++){
                      BLoc[i] = inputPositions[ NodeBIndex*SPACE_DIM + i ];
                  }
                  double BRad = neighbour_iter->GetRadius();
                  

                  // Dirac delta. There is only a contribution if we're differentiating wrt A or B's position.
                  bool N_equals_A = false;
                  bool N_equals_B = false;
                  if( NodeAIndex == derivCellN ){
                      N_equals_A = true;
                  }
                  if( NodeBIndex == derivCellN ){
                      N_equals_B = true;
                  }

    
                  if(N_equals_A || N_equals_B){
                      
                      // Work out whether this neighbour is close enough to actually contribute to the force on A,
                      // given the current input positions
                      c_vector<double, SPACE_DIM> AtoB = this->rCellPopulation.rGetMesh().GetVectorFromAtoB(ALoc, BLoc);
                      double separation = norm_2(AtoB); 
                      c_vector<double, SPACE_DIM> AtoB_unit = AtoB/separation;
                      double overlap = separation - ARad - BRad;
                      //std::cout << "Overlap " << overlap << std::endl;
    
                      bool makesContribution = true;
                      if(useCutoff){
                          if(separation > cutoff){
                              makesContribution = false;
                              std::cout << "No contribution" << std::endl;
                          }
                      }
    
                      // Determine Jacobian contribution 
                      if(makesContribution){
    
                          c_vector<double,SPACE_DIM> dSeparation;
                          dSeparation[0] = 0;
                          dSeparation[1] = 0;
                          dSeparation[2] = 0;
                          if(N_equals_A){
                              dSeparation[derivCptW] = -AtoB[derivCptW]/separation;   
                          }
                          if(N_equals_B){
                              dSeparation[derivCptW] = AtoB[derivCptW]/separation; 
                          }
    
                          c_vector<double, SPACE_DIM> dAtoB_unit;
                          dAtoB_unit[0] = 0;
                          dAtoB_unit[1] = 0;
                          dAtoB_unit[2] = 0;
                          if(N_equals_A){
                            for(int i=0; i<SPACE_DIM; i++){
                              if(i == derivCptW){
                                dAtoB_unit[i] = (-separation -  dSeparation[i] * AtoB[i]) / (separation*separation);
                              }else{
                                dAtoB_unit[i] = (-dSeparation[i] * AtoB[i]) / (separation*separation);
                              }
                            }
                          }
                          if(N_equals_B){
                            for(int i=0; i<SPACE_DIM; i++){
                              if(i == derivCptW){
                                dAtoB_unit[i] = (separation -  dSeparation[i] * AtoB[i]) / (separation*separation);
                              }else{
                                dAtoB_unit[i] = (-dSeparation[i] * AtoB[i]) / (separation*separation);
                              }
                            }
                          }
    
                          if(overlap < 0){
                            //Log part of the force law applies
                            for(int i=0; i<SPACE_DIM; i++){
                              dFA_dNw[i] += springConst*(ARad+BRad)*( log(1+overlap/(ARad+BRad)) * dAtoB_unit[i] +
                                                                       dSeparation[i] * ((ARad+BRad)/separation) * AtoB_unit[i] );
                            }
                          }else{
                            //Exponential part of the force law applies
                            for(int i=0; i<SPACE_DIM; i++){
                              dFA_dNw[i] += springConst*(ARad+BRad)*( dSeparation[i]*( exp(-5*overlap/(ARad+BRad)) * AtoB_unit[i] ) + 
                                                                        separation*(-5*dSeparation[i] * exp(-5*overlap/(ARad+BRad)) *AtoB_unit[i] +  exp(-5*overlap/(ARad+BRad)) *dAtoB_unit[i])  );
                            }
                          }
                      }
                  }
              }
          }

          for(int i=0; i<SPACE_DIM; i++){
            int JIndex = NodeAIndex*SPACE_DIM + i;
            //std::cout << "i " << JIndex << " j " << deriv_index << " " << dFA_dNw[i] << std::endl;
            PetscMatTools::SetElement(*pJacobian, JIndex, deriv_index, dFA_dNw[i]);
          }
          
        }    
    }

    PetscMatTools::Finalise(*pJacobian); 
};


///////// Explicit instantiation
template class BackwardEulerNumericalMethodTimestepper<1,1>;
template class BackwardEulerNumericalMethodTimestepper<1,2>;
template class BackwardEulerNumericalMethodTimestepper<2,2>;
template class BackwardEulerNumericalMethodTimestepper<1,3>;
template class BackwardEulerNumericalMethodTimestepper<2,3>;
template class BackwardEulerNumericalMethodTimestepper<3,3>;