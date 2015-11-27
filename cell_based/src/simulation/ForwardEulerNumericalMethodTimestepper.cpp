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

#include "ForwardEulerNumericalMethodTimestepper.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>	
ForwardEulerNumericalMethodTimestepper<ELEMENT_DIM,SPACE_DIM> :: ForwardEulerNumericalMethodTimestepper( 
	                                               AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>&                  inputCellPopulation, 
                                                 std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >&  inputForceCollection)
:AbstractNumericalMethodTimestepper<ELEMENT_DIM,SPACE_DIM> ( inputCellPopulation, inputForceCollection)
{	
};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>  
ForwardEulerNumericalMethodTimestepper<ELEMENT_DIM,SPACE_DIM>::~ForwardEulerNumericalMethodTimestepper(){
};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>  
void ForwardEulerNumericalMethodTimestepper<ELEMENT_DIM,SPACE_DIM>::UpdateAllNodePositions(double dt){

    if(this->nonEulerSteppersEnabled){

        std::vector<c_vector<double, SPACE_DIM> > F = this->ComputeAndSaveForces();

        int index = 0;
        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rCellPopulation.rGetMesh().GetNodeIteratorBegin();
             node_iter != this->rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
        {
            c_vector<double, SPACE_DIM> oldLocation = node_iter->rGetLocation();
            c_vector<double, SPACE_DIM> displacement = dt * F[index];

            this->HandleStepSizeExceptions(&displacement, dt, node_iter->GetIndex());

            c_vector<double, SPACE_DIM> newLocation = oldLocation + displacement;
                            
            ChastePoint<SPACE_DIM> new_point(newLocation);
            this->rCellPopulation.SetNode(node_iter->GetIndex(), new_point);
        }

    }else{

        this->ComputeAndSaveForces();
        this->rCellPopulation.UpdateNodeLocations(dt);    

    }

};

///////// Explicit instantiation
template class ForwardEulerNumericalMethodTimestepper<1,1>;
template class ForwardEulerNumericalMethodTimestepper<1,2>;
template class ForwardEulerNumericalMethodTimestepper<2,2>;
template class ForwardEulerNumericalMethodTimestepper<1,3>;
template class ForwardEulerNumericalMethodTimestepper<2,3>;
template class ForwardEulerNumericalMethodTimestepper<3,3>;