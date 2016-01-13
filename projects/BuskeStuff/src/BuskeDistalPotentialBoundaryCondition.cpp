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

#include "BuskeDistalPotentialBoundaryCondition.hpp"
#include "NodeBasedCellPopulationWithBuskeUpdate.hpp"
#include "IsNan.hpp"


template<unsigned DIM>
BuskeDistalPotentialBoundaryCondition<DIM>::BuskeDistalPotentialBoundaryCondition(double inLength, double inRadius, double inStrength)

   : AbstractForce<DIM>(),
     length(inLength),
     radius(inRadius),
     strength(inStrength)
{ 
}


template<unsigned DIM>
const double BuskeDistalPotentialBoundaryCondition<DIM>::GetLength() const{
    return length;
}
template<unsigned DIM>
const double BuskeDistalPotentialBoundaryCondition<DIM>::GetStrength() const{
    return strength;
}
template<unsigned DIM>
const double BuskeDistalPotentialBoundaryCondition<DIM>::GetRadius() const{
    return radius;
}



template<unsigned DIM>
void BuskeDistalPotentialBoundaryCondition<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{

    double maxMagnitude = 0;

    //Loop over cells
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {   
    
        //Retrieve data
        unsigned nodeIndex = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        Node<DIM>* node = rCellPopulation.GetNode(nodeIndex);
        c_vector<double, DIM> cellLocation = node->rGetLocation();
        double cellRadius = node->GetRadius();

        // Apply boundary
        if(cellLocation[0] > length){    
            
            // Not subject to boundary condition, beyond end of tube. Skip.
        
        }else{

            // Determine distance from the boundary midline and outward normal
            double midlineDist;
            c_vector<double,DIM> out;

            if(cellLocation[0] < 0){
            
                // Endcap cell
                midlineDist = norm_2(cellLocation);
                out = cellLocation/midlineDist;

            }else{

                midlineDist = sqrt(cellLocation[1]*cellLocation[1]+cellLocation[2]*cellLocation[2]);
                out[0] = 0;
                out[1] = cellLocation[1]/midlineDist;
                out[2] = cellLocation[2]/midlineDist;
            }

            double targetDist = radius - cellRadius;
            double diff = targetDist - midlineDist;

            if(diff == 0){

                // Do nothing, cell perfectly positioned 
            
            }else{

                // Push toward surface of crypt / arm

                double magnitude = strength*pow(2.71828,( (1/cellRadius)-(1/fabs(diff)) )) * diff/fabs(diff);
                c_vector<double, DIM> force = magnitude*out;
                node->AddAppliedForceContribution(force);
    
                if(magnitude > maxMagnitude){
                    maxMagnitude = magnitude;
                }   

            }
        }
    }

    //std::cout << "Max applied boundary force: " << maxMagnitude << std::endl;
}


template<unsigned DIM>
void BuskeDistalPotentialBoundaryCondition<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<BoundaryLength>" << length << "</BoundaryLength>\n";
    *rParamsFile << "\t\t\t<BoundaryRadius>" << radius << "</BoundaryRadius>\n";
    *rParamsFile << "\t\t\t<BoundaryForceStrength>" << strength << "</BoundaryForceStrength>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class BuskeDistalPotentialBoundaryCondition<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BuskeDistalPotentialBoundaryCondition)
