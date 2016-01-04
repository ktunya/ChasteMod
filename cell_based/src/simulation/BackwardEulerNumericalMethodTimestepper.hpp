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

#ifndef BACKWARDEULERNUMERICALMETHODTIMESTEPPER_HPP_
#define BACKWARDEULERNUMERICALMETHODTIMESTEPPER_HPP_

#include "AbstractNumericalMethodTimestepper.hpp"
#include "SimplePetscNonlinearSolver.hpp"
#include "PetscMatTools.hpp"
#include "PetscVecTools.hpp"
#include "PetscTools.hpp"



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class BackwardEulerNumericalMethodTimestepper : public AbstractNumericalMethodTimestepper<ELEMENT_DIM,SPACE_DIM> {

public:	

	SimplePetscNonlinearSolver* pNonlinearSolver;
	double implicitStepSize;

	
	BackwardEulerNumericalMethodTimestepper(AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>&                  inputCellPopulation, 
                                          std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >&  inputForceCollection);

	virtual ~BackwardEulerNumericalMethodTimestepper();

  virtual void UpdateAllNodePositions(double dt);

	void BACKWARDEULERComputeResidual(const Vec currentGuess, Vec residualVector);

  void BACKWARDEULERComputeSpringJacobian(const Vec currentGuess, Mat* pJacobian);

};





template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PetscErrorCode BACKWARDEULER_ComputeResidual(SNES snes, Vec currentGuess, Vec residualVector, void* pContext){

    BackwardEulerNumericalMethodTimestepper<ELEMENT_DIM, SPACE_DIM>* pStepper = (BackwardEulerNumericalMethodTimestepper<ELEMENT_DIM, SPACE_DIM>*)pContext; 
    pStepper->BACKWARDEULERComputeResidual(currentGuess, residualVector);

    return 0;
};




#if ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=5 )

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PetscErrorCode BACKWARDEULER_ComputeSpringJacobian(SNES snes, Vec input, Mat jacobian, Mat preconditioner, void* pContext){

    BackwardEulerNumericalMethodTimestepper<ELEMENT_DIM, SPACE_DIM>* pStepper = (BackwardEulerNumericalMethodTimestepper<ELEMENT_DIM, SPACE_DIM>*)pContext; 
    pStepper->BACKWARDEULERComputeSpringJacobian(input, &jacobian);

    return 0;
};

/*
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PetscErrorCode BACKWARDEULER_ComputeJacobianComparison(SNES snes, Vec input, Mat* pJacobian, Mat* pPreconditioner, MatStructure* pMatStructure, void* pContext){

    BackwardEulerNumericalMethodTimestepper<ELEMENT_DIM, SPACE_DIM>* pStepper = (BackwardEulerNumericalMethodTimestepper<ELEMENT_DIM, SPACE_DIM>*)pContext; 
    
    pStepper->BACKWARDEULERComputeSpringJacobian(input, pJacobian);

    Mat JacobianReference;
    PetscTools::SetupMat(JacobianReference, num_unknowns, num_unknowns, UINT_MAX, PETSC_DECIDE, PETSC_DECIDE, true, false);
    SNESComputeJacobianDefault(snes, input, &JacobianReference, &JacobianReference, pMatStructure, NULL);   

    unsigned num_unknowns = PetscVecTools::GetSize(input);
    if (!PetscMatTools::CheckEquality(*pJacobian, JacobianReference, 1e-6)) {
      std::cout << "INCORRECT ANALYTIC JACOBIAN" << std::endl;
      for (int i=0; i<num_unknowns; i++) {
        for (int j=0; j<num_unknowns; j++) {
          std::cout << "Jcpt " << i << " " << j << " ref: " << PetscMatTools::GetElement(JacobianReference,i,j) << " actual: " << PetscMatTools::GetElement(*pJacobian,i,j) << std::endl;
        }
      }
      EXCEPTION("Jacobian exception.");
    }

    return 0;
};*/

#else

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PetscErrorCode BACKWARDEULER_ComputeSpringJacobian(SNES snes, Vec input, Mat* pJacobian, Mat* pPreconditioner, MatStructure* pMatStructure, void* pContext){

    BackwardEulerNumericalMethodTimestepper<ELEMENT_DIM, SPACE_DIM>* pStepper = (BackwardEulerNumericalMethodTimestepper<ELEMENT_DIM, SPACE_DIM>*)pContext; 
    pStepper->BACKWARDEULERComputeSpringJacobian(input, pJacobian);

    return 0;
};

/*
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PetscErrorCode BACKWARDEULER_ComputeJacobianComparison(SNES snes, Vec input, Mat* pJacobian, Mat* pPreconditioner, MatStructure* pMatStructure, void* pContext){

    unsigned num_unknowns = PetscVecTools::GetSize(input);

    BackwardEulerNumericalMethodTimestepper<ELEMENT_DIM, SPACE_DIM>* pStepper = (BackwardEulerNumericalMethodTimestepper<ELEMENT_DIM, SPACE_DIM>*)pContext; 
    
    pStepper->BACKWARDEULERComputeSpringJacobian(input, pJacobian);

    Mat JacobianReference;
    PetscTools::SetupMat(JacobianReference, num_unknowns, num_unknowns, UINT_MAX, PETSC_DECIDE, PETSC_DECIDE, true, false);
    SNESComputeJacobianDefault(snes, input, &JacobianReference, &JacobianReference, pMatStructure, pContext);   

    if(!PetscMatTools::CheckEquality(*pJacobian, JacobianReference, 1e-6)){
      std::cout << "INCORRECT ANALYTIC JACOBIAN" << std::endl;
      for (int i=0; i<num_unknowns; i++) {
        for (int j=0; j<num_unknowns; j++) {
          std::cout << "Jcpt " << i << " " << j << " ref: " << PetscMatTools::GetElement(JacobianReference,i,j) << " actual: " << PetscMatTools::GetElement(*pJacobian,i,j) << std::endl;
        }
      }
      //EXCEPTION("Jacobian exception.");
    }

    return 0;
};
*/

#endif


#endif /*BACKWARDEULERNUMERICALMETHODTIMESTEPPER_HPP_*/