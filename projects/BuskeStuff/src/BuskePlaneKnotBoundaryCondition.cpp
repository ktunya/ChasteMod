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

#include "BuskePlaneKnotBoundaryCondition.hpp"
#include "NodeBasedCellPopulationWithBuskeUpdate.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

template<unsigned DIM>
BuskePlaneKnotBoundaryCondition<DIM>::BuskePlaneKnotBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation,
                                                    c_vector<double, DIM> BottomLeft,
                                                    c_vector<double, DIM> TopRight,
                                                    c_vector<double, DIM> Normal,
                                                    double KnotSpacing)
        : AbstractCellPopulationBoundaryCondition<DIM>(pCellPopulation),
          bottomLeft(BottomLeft),
          topRight(TopRight),
          knotSpacing(KnotSpacing)
{
    // Assert checks
    assert(norm_2(Normal) > 0.0);
    normal = Normal/norm_2(Normal);
    assert( (dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(this->mpCellPopulation)) );

    //Setup system of Buske knots
    c_vector<double,DIM> cellPos = bottomLeft;
    int nKnots = 0;

    if(DIM==2){

        c_vector<double,DIM> slopeVec = (topRight-bottomLeft)/norm_2(topRight-bottomLeft);

        while(cellPos[0]<=topRight[0] && cellPos[1]<=topRight[1]){
            
            AddKnot(cellPos);
            nKnots++;
            cellPos += slopeVec*knotSpacing;
        }

    }else if(DIM==3){

        c_vector<double,DIM> slopeVec1 = topRight - bottomLeft;
        slopeVec1 = slopeVec1/norm_2(slopeVec1);

        c_vector<double,DIM> slopeVec2;
        slopeVec2[0] = slopeVec1[1]*normal[2] - slopeVec1[2]*normal[1];  
        slopeVec2[1] = slopeVec1[2]*normal[0] - slopeVec1[0]*normal[2];  
        slopeVec2[2] = slopeVec1[0]*normal[1] - slopeVec1[1]*normal[0];  
        slopeVec2 = slopeVec2/norm_2(slopeVec2);

        int counter = 0;

        while(cellPos[0]<=topRight[0] &&   cellPos[1]<=topRight[1] && cellPos[2]<=topRight[2] &&
              cellPos[0]>=bottomLeft[0] && cellPos[1]>=bottomLeft[1] && cellPos[2]>=bottomLeft[2]){

            c_vector<double, DIM> pos = cellPos;

            while(pos[0]<=topRight[0] && pos[1]<=topRight[1] && pos[2]<=topRight[2] &&
                  pos[0]>=bottomLeft[0] && pos[1]>=bottomLeft[1] && pos[2]>=bottomLeft[2]){    

                AddKnot(pos);
                nKnots++;

                pos += slopeVec2*knotSpacing;
            }

            pos = cellPos-slopeVec2*knotSpacing;

            while(pos[0]<=topRight[0] && pos[1]<=topRight[1] && pos[2]<=topRight[2] &&
                  pos[0]>=bottomLeft[0] && pos[1]>=bottomLeft[1] && pos[2]>=bottomLeft[2]){    

                AddKnot(pos);
                nKnots++;

                pos -= slopeVec2*knotSpacing;
            }
            
            cellPos += slopeVec1*knotSpacing;
        }


    }else{
        EXCEPTION("Illegal dimension in BuskePlaneKnotBoundaryCondition");
    }

    int knots = dynamic_cast<NodeBasedCellPopulationWithBuskeUpdate<DIM>* >(this->mpCellPopulation)->GetNumKnots();
    dynamic_cast<NodeBasedCellPopulationWithBuskeUpdate<DIM>* >(this->mpCellPopulation)->SetNumKnots(knots + nKnots);
    std::cout << "Knots so far: " << knots << " new knots: " << nKnots << std::endl;

}

template<unsigned DIM>
void BuskePlaneKnotBoundaryCondition<DIM>::AddKnot(c_vector<double,DIM> cellPos){

    FixedDurationGenerationBasedCellCycleModel* pCCM = new FixedDurationGenerationBasedCellCycleModel;
    pCCM->SetDimension(DIM);

    boost::shared_ptr<AbstractCellProperty> pState(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
            
    CellPtr pCell(new Cell(pState, pCCM));
           
    pCell->SetCellProliferativeType(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());

    CellPtr result = this->mpCellPopulation->AddCell( pCell, cellPos, *(this->mpCellPopulation->Begin()) );
    result->GetCellData()->SetItem("Radius",knotSpacing/2.0);
    result->GetCellData()->SetItem("RelaxedRadius",knotSpacing/2.0);
    result->GetCellData()->SetItem("IsBuskeKnot",1);
};


template<unsigned DIM>
const c_vector<double, DIM>& BuskePlaneKnotBoundaryCondition<DIM>::rGetBottomLeftOfPlane() const
{
    return bottomLeft;
}

template<unsigned DIM>
const c_vector<double, DIM>& BuskePlaneKnotBoundaryCondition<DIM>::rGetTopRightOfPlane() const
{
    return topRight;
}

template<unsigned DIM>
const c_vector<double, DIM>& BuskePlaneKnotBoundaryCondition<DIM>::rGetNormalToPlane() const
{
    return normal;
}

template<unsigned DIM>
const double& BuskePlaneKnotBoundaryCondition<DIM>::rGetKnotSpacing() const{
    return knotSpacing;
}


template<unsigned DIM>
void BuskePlaneKnotBoundaryCondition<DIM>::ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations)
{
    /* 
    * Does nothing. Provided this class is used alongside BuskeKnotForce, each knot will exert a repulsion on the 
    * real cells, which should keep their population on the right ride of the boundary. Therefore, ImposeBoundaryCondition is
    * blank, but this class will still verify that the boundary is holding.
    */
}


template<unsigned DIM>
bool BuskePlaneKnotBoundaryCondition<DIM>::VerifyBoundaryCondition()
{
    bool condition_satisfied = true;

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
        if(cell_iter->GetCellData()->GetItem("IsBuskeKnot")==1){
        }else{
        
            c_vector<double, DIM> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);

            if (inner_prod(cell_location - bottomLeft, normal) > 0.0)
            {
                condition_satisfied = false;
                break;
            }
        }
    }

    return condition_satisfied;
}


template<unsigned DIM>
void BuskePlaneKnotBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<BottomLeftOfPlane>";
    for (unsigned index=0; index != DIM-1U; index++) 
    {
        *rParamsFile << bottomLeft[index] << ",";
    }
    *rParamsFile << bottomLeft[DIM-1] << "</BottomLeftOfPlane>\n";

    *rParamsFile << "\t\t\t<TopRightOfPlane>";
    for (unsigned index=0; index != DIM-1U; index++) 
    {
        *rParamsFile << topRight[index] << ",";
    }
    *rParamsFile << topRight[DIM-1] << "</TopRightOfPlane>\n";

    *rParamsFile << "\t\t\t<NormalToPlane>";
    for (unsigned index=0; index != DIM-1U; index++) 
    {
        *rParamsFile << normal[index] << ",";
    }
    *rParamsFile << normal[DIM-1] << "</NormalToPlane>\n";

    *rParamsFile << "\t\t\t<KnotSpacing>" << knotSpacing << "</KnotSpacing>\n";

    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class BuskePlaneKnotBoundaryCondition<2>;
template class BuskePlaneKnotBoundaryCondition<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BuskePlaneKnotBoundaryCondition)
