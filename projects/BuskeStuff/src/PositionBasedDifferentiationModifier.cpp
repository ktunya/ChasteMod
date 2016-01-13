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

#include "PositionBasedDifferentiationModifier.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

template<unsigned DIM>
PositionBasedDifferentiationModifier<DIM>::PositionBasedDifferentiationModifier(c_vector<double, DIM> pointOnDiffPlane, 
                                                                                c_vector<double, DIM> normToDiffPlane)
    : AbstractCellBasedSimulationModifier<DIM>()
{
    mPointOnPlane = pointOnDiffPlane;
    mNormalToPlane = normToDiffPlane;
}

template<unsigned DIM>
PositionBasedDifferentiationModifier<DIM>::~PositionBasedDifferentiationModifier(){}



template<unsigned DIM>
void PositionBasedDifferentiationModifier<DIM>::SetPointOnPlane(c_vector<double, DIM> point){
    mPointOnPlane = point;
};

template<unsigned DIM>
void PositionBasedDifferentiationModifier<DIM>::SetNormalToPlane(c_vector<double, DIM> normal){
    mNormalToPlane = normal;
};

template<unsigned DIM>
const c_vector<double, DIM> PositionBasedDifferentiationModifier<DIM>::GetPointOnPlane() const {
    return mPointOnPlane;
};

template<unsigned DIM>
const c_vector<double, DIM> PositionBasedDifferentiationModifier<DIM>::GetNormalToPlane() const {
    return mNormalToPlane;
};



template<unsigned DIM>
void PositionBasedDifferentiationModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure cell locations are updated
    rCellPopulation.Update();

    // Assess whether each cell has crossed the plane and if so differentiate it 
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        unsigned nodeId = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        Node<DIM>* nodePtr = rCellPopulation.GetNode(nodeId);

        c_vector<double, DIM> location = nodePtr->rGetLocation();
        c_vector<double, DIM> diff = location - mPointOnPlane;
        double dotProduct = 0;
        for(int i=0; i<3; i++){
            dotProduct += diff[i] * mNormalToPlane[i]; 
        } 
        if(dotProduct > 0){
            MAKE_PTR(DifferentiatedCellProliferativeType, pDiffState);
            cell_iter->SetCellProliferativeType(pDiffState);
            cell_iter->GetCellData()->SetItem("Differentiated",1);
        } 
    }  
}

template<unsigned DIM>
void PositionBasedDifferentiationModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    // Make sure cell locations are updated
    rCellPopulation.Update();

    // Assess whether each cell has crossed the plane and if so differentiate it 
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        unsigned nodeId = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        Node<DIM>* nodePtr = rCellPopulation.GetNode(nodeId);

        c_vector<double, DIM> location = nodePtr->rGetLocation();
        c_vector<double, DIM> diff = location - mPointOnPlane;
        double dotProduct = 0;
        for(int i=0; i<3; i++){
            dotProduct += diff[i] * mNormalToPlane[i]; 
        } 
        if(dotProduct > 0){
            MAKE_PTR(DifferentiatedCellProliferativeType, pDiffState);
            cell_iter->SetCellProliferativeType(pDiffState);
            cell_iter->GetCellData()->SetItem("Differentiated",1.0);
        } 
    }   
}


template<unsigned DIM>
void PositionBasedDifferentiationModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{

    *rParamsFile << "\t\t\t<PointOnDiffPlane>";
    for (unsigned index=0; index != DIM-1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mPointOnPlane[index] << ",";
    }
    *rParamsFile << mPointOnPlane[DIM-1] << "</PointOnDiffPlane>\n";

    *rParamsFile << "\t\t\t<NormalToDiffPlane>";
    for (unsigned index=0; index != DIM-1U; index++) // Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mNormalToPlane[index] << ",";
    }
    *rParamsFile << mNormalToPlane[DIM-1] << "</NormalToDiffPlane>\n";
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class PositionBasedDifferentiationModifier<1>;
template class PositionBasedDifferentiationModifier<2>;
template class PositionBasedDifferentiationModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PositionBasedDifferentiationModifier)
