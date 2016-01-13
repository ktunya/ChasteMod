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

#include "PlaneBasedCellKillerWithRecording.hpp"

template<unsigned DIM>
PlaneBasedCellKillerWithRecording<DIM>::PlaneBasedCellKillerWithRecording(AbstractCellPopulation<DIM>* pCellPopulation,
                                                  c_vector<double, DIM> point,
                                                  c_vector<double, DIM> normal,
                                                  std::string outputDir)
    : PlaneBasedCellKiller<DIM>(pCellPopulation, point, normal)
{
  outDir = outputDir;

  OutputFileHandler rOutputFileHandler(outputDir, false);
  outputFile = rOutputFileHandler.OpenOutputFile("PlaneDeathsData.txt");
}




template<unsigned DIM>
void PlaneBasedCellKillerWithRecording<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
        c_vector<double, DIM> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);
        c_vector<double, DIM> diff = cell_location - this->rGetPointOnPlane();
        c_vector<double, DIM> norm = this->rGetNormalToPlane();
        double innerProduct = 0;
        for(int i=0; i<DIM; i++){
          innerProduct += norm[i]*diff[i];
        } 

        if (innerProduct > 0.0)
        {
          *outputFile  << SimulationTime::Instance()->GetTime() << "\t" << cell_iter->GetCellId() << "\n";
          cell_iter->Kill();
        }
    }

    outputFile->flush();

    //If simulation is finished, close the output file.
    if(SimulationTime::Instance()->IsFinished()){
      outputFile->close();
    }
}


template<unsigned DIM>
const std::string PlaneBasedCellKillerWithRecording<DIM>::GetOutputDir() const{
  return outDir;
};


// Explicit instantiation
template class PlaneBasedCellKillerWithRecording<1>;
template class PlaneBasedCellKillerWithRecording<2>;
template class PlaneBasedCellKillerWithRecording<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PlaneBasedCellKillerWithRecording)
