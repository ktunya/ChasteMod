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

#include "BuskeKnotForce.hpp"
#include "NodeBasedCellPopulationWithBuskeUpdate.hpp"
#include "IsNan.hpp"


template<unsigned DIM>
BuskeKnotForce<DIM>::BuskeKnotForce(double interaction, double omegaRatio)

   : AbstractForce<DIM>(),
     maxInteractionEnergy(interaction),
     thresholdAdhesionRatio(omegaRatio)
{ 
}


template<unsigned DIM>
void BuskeKnotForce<DIM>::SetMaxInteractionEnergy( double interaction )
{
    maxInteractionEnergy = interaction;
}
template<unsigned DIM>
void BuskeKnotForce<DIM>::SetThresholdAdhesionRatio( double omegaRatio )
{
    thresholdAdhesionRatio = omegaRatio;   
}


template<unsigned DIM>
const double BuskeKnotForce<DIM>::GetMaxInteractionEnergy() const{
    return maxInteractionEnergy;
}
template<unsigned DIM>
const double BuskeKnotForce<DIM>::GetThresholdAdhesionRatio() const{
    return thresholdAdhesionRatio;
}



template<unsigned DIM>
void BuskeKnotForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{

    int numKnots = dynamic_cast<NodeBasedCellPopulationWithBuskeUpdate<DIM>* >(&rCellPopulation)->GetNumKnots();

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {   
    
        //If this ia a cell...
        if(cell_iter->GetCellData()->GetItem("IsBuskeKnot")==0){

            c_vector<double,DIM> appliedForce;
            for(int i=0; i<DIM; i++){
                appliedForce[i]=0;
            }
            int interactingKnots = 0;

            // Retrieve neighbours
            unsigned nodeIndex = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            Node<DIM>* node = rCellPopulation.GetNode(nodeIndex);
            c_vector<double, DIM> cellLocation = node->rGetLocation();
            
            std::set<unsigned> neighbours = rCellPopulation.GetNeighbouringLocationIndices(*cell_iter);

            // Loop over neighbours and look for knots
            std::set<unsigned>::iterator it;
            for(it = neighbours.begin(); it != neighbours.end(); ++it){

                double isKnot = (rCellPopulation.GetCellUsingLocationIndex(*it))->GetCellData()->GetItem("IsBuskeKnot");

                if(isKnot == 1){

                    Node<DIM>* knotNode = rCellPopulation.GetNode(*it);
                    c_vector<double, DIM> knotLocation = knotNode->rGetLocation();
                    double radius = knotNode->GetRadius();
                
                    //Check whether the knot is close enough to interact with this cell
                    double separation = norm_2(cellLocation - knotLocation);
    
                    if(separation < radius){

                        interactingKnots++;
    
                        c_vector<double, DIM> unitSep = (cellLocation - knotLocation)/separation;
                        double magnitude = maxInteractionEnergy*(thresholdAdhesionRatio*(1/separation)-(1/radius));
    
                        appliedForce = appliedForce + unitSep * magnitude;                    
                    }

                }
            }

            if(interactingKnots != 0){
                appliedForce = appliedForce / interactingKnots; 
            }
            node->AddAppliedForceContribution(appliedForce); 
        }
    }
}


template<unsigned DIM>
void BuskeKnotForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{

    *rParamsFile << "\t\t\t<MaxInteractionEnergy>" << maxInteractionEnergy << "</MaxInteractionEnergy>\n";
    *rParamsFile << "\t\t\t<ThresholdAdhesionRatio>" << thresholdAdhesionRatio << "</ThresholdAdhesionRatio>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

//Only 3D working for now.
template class BuskeKnotForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BuskeKnotForce)
