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

#include "BuskePlanePotentialBoundaryCondition.hpp"
#include "NodeBasedCellPopulationWithBuskeUpdate.hpp"
#include "IsNan.hpp"


template<unsigned DIM>
BuskePlanePotentialBoundaryCondition<DIM>::BuskePlanePotentialBoundaryCondition(c_vector<double,DIM> pointOnPlane, c_vector<double,DIM> normalToPlane, double interactionE, double threshAdhesionRatio)

   : AbstractForce<DIM>(),
     point(pointOnPlane),
     normal(normalToPlane / norm_2(normalToPlane)),
     interactionEnergy(interactionE),
     thresholdAdhesionRatio(threshAdhesionRatio)
{ 
}


template<unsigned DIM>
const double BuskePlanePotentialBoundaryCondition<DIM>::GetInteractionEnergy() const{
    return interactionEnergy;
}
template<unsigned DIM>
const double BuskePlanePotentialBoundaryCondition<DIM>::GetThresholdAdhesionRatio() const{
    return thresholdAdhesionRatio;
}
template<unsigned DIM>
const c_vector<double,DIM> BuskePlanePotentialBoundaryCondition<DIM>::GetNormal() const{
    return normal;
}
template<unsigned DIM>
const c_vector<double,DIM> BuskePlanePotentialBoundaryCondition<DIM>::GetPoint() const{
    return point;
}



template<unsigned DIM>
void BuskePlanePotentialBoundaryCondition<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{

    //Loop over cells
    for (typename AbstractCellPopulation<DIM>::RealCellsIterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {   
    
        //Retrieve data
        unsigned nodeIndex = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        Node<DIM>* node = rCellPopulation.GetNode(nodeIndex);
        c_vector<double, DIM> cellLocation = node->rGetLocation();
        double radius = node->GetRadius();

        //Determine min distance from plane
        c_vector<double, DIM> planeToCell = cellLocation - point; 
        c_vector<double, DIM> normalProjection;
        double separation = 0;

        for(int i=0; i<DIM; i++){
            normalProjection[i] = planeToCell[i] * normal[i]; 
            separation += planeToCell[i] * normal[i];
        } 

        if(separation > -radius){

            std::cout << "Apply force sep: " << separation << std::endl;

            //Apply force
            double magnitude = 100*pow(2.71828,separation);
            //double magnitude = interactionEnergy * ((thresholdAdhesionRatio/fabs(separation)) - (1/radius));
            c_vector<double, DIM> force = -normal*magnitude;
            node->AddAppliedForceContribution(force);

            std::cout << "force: " << normalProjection[0] << " " << normalProjection[1] << " " << normalProjection[2] << std::endl;

            //Cell has escaped the boundary
            //std::cout << "Outside sep: " << separation << std::endl;
            //std::cout << "Cell has crossed the Buske plane boundary. Try increasing epsilon or reducing the time step." << std::endl;
            //EXCEPTION("Cell has crossed the Buske plane boundary");
        }
    }
}


template<unsigned DIM>
void BuskePlanePotentialBoundaryCondition<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    for(int i=0; i<DIM; i++){
        *rParamsFile << "\t\t\t<PointOnPlane" << i << ">" << point[i] << "</PointOnPlane" << i << ">\n";
    }
    for(int i=0; i<DIM; i++){
        *rParamsFile << "\t\t\t<Normal" << i << ">" << normal[i] << "</Normal" << i << ">\n";
    }
    *rParamsFile << "\t\t\t<InteractionEnergy>" << interactionEnergy << "</InteractionEnergy>\n";
    *rParamsFile << "\t\t\t<ThresholdAdhesionRatio>" << thresholdAdhesionRatio << "</ThresholdAdhesionRatio>\n";

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class BuskePlanePotentialBoundaryCondition<1>;
template class BuskePlanePotentialBoundaryCondition<2>;
template class BuskePlanePotentialBoundaryCondition<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BuskePlanePotentialBoundaryCondition)
