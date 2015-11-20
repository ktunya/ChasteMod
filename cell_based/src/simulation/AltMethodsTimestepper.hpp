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

#ifndef ALTMETHODSTIMESTEPPER_HPP_
#define ALTMETHODSTIMESTEPPER_HPP_

#include "AbstractOffLatticeCellPopulation.hpp"
#include "AbstractForce.hpp"
#include "SimplePetscNonlinearSolver.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class AltMethodsTimestepper {

private:

	AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>&                  rCellPopulation;
	std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >&  rForceCollection;

	SimplePetscNonlinearSolver* pNonlinearSolver;
	double implicitStepSize;

	std::vector<c_vector<double, SPACE_DIM> > ComputeAndSaveForces();
	std::vector<c_vector<double, SPACE_DIM> > SaveCurrentLocations();

	bool nonEulerSteppersEnabled; 

public:	
	
	AltMethodsTimestepper(AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>&  inputCellPopulation, 
                std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >&  inputForceCollection);

	~AltMethodsTimestepper();

	void UpdateAllNodePositions(double dt, int stepper);

	void BACKWARDEULERComputeResidual(const Vec currentGuess, Vec residualVector);
	
};


/*
* Global residual functions for use with implicit steppers
*/
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
PetscErrorCode BACKWARDEULER_ComputeResidual(SNES snes, Vec currentGuess, Vec residualVector, void* pContext){

    AltMethodsTimestepper<ELEMENT_DIM, SPACE_DIM>* pStepper = (AltMethodsTimestepper<ELEMENT_DIM, SPACE_DIM>*)pContext; 
    pStepper->BACKWARDEULERComputeResidual(currentGuess, residualVector);

    return 0;
};

#endif /*ALTMETHODSTIMESTEPPER_HPP_*/