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

#ifndef POPULATIONTESTINGFORCE_HPP_
#define POPULATIONTESTINGFORCE_HPP_

#include "AbstractForce.hpp"

/*
* A simple test force used in several off the cell population tests
*/


template<unsigned  ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class PopulationTestingForce : public AbstractForce<ELEMENT_DIM, SPACE_DIM> {

public:

	PopulationTestingForce(): 
	         AbstractForce<ELEMENT_DIM, SPACE_DIM>()
	{
	}

	virtual void AddForceContribution(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation){	

		for (unsigned i=0; i<rCellPopulation.GetNumNodes(); i++)
        {
            c_vector<double, SPACE_DIM> force;

            for(int j=0; j<SPACE_DIM; j++){
            	
            	force[j] = (j+1)*i*0.01;
            }

            rCellPopulation.GetNode(i)->ClearAppliedForce();
            rCellPopulation.GetNode(i)->AddAppliedForceContribution(force);
        }

	};

	virtual void OutputForceParameters(out_stream& rParamsFile){		
	};

};

#endif /*POPULATIONTESTINGFORCE_HPP_*/ 