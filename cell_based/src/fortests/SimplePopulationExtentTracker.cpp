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

#include "SimplePopulationExtentTracker.hpp"


//Constructor, initialises sampling interval and sets the output file to null
template<unsigned DIM>
SimplePopulationExtentTracker<DIM>::SimplePopulationExtentTracker(int samplingInterval)
    : AbstractCellBasedSimulationModifier<DIM>(),
      outputFile(NULL),
      outputDir(std::string()),
      interval(samplingInterval)
{}

//Empty destructor 
template<unsigned DIM>
SimplePopulationExtentTracker<DIM>::~SimplePopulationExtentTracker(){}



//Getters 
template<unsigned DIM>
int SimplePopulationExtentTracker<DIM>::GetInterval() const
{
  return interval;
};

template<unsigned DIM>
std::string SimplePopulationExtentTracker<DIM>::GetOutputDirectoryFull(){
  return outputDir;
};



//Open an output file ExtentData.txt in the simulation directory
template<unsigned DIM>
void SimplePopulationExtentTracker<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
  OutputFileHandler rOutputFileHandler(outputDirectory, false);
  outputDir = rOutputFileHandler.GetOutputDirectoryFullPath();
  outputFile = rOutputFileHandler.OpenOutputFile("ExtentData.txt");
}



//At each timestep, if the time is a sampling time, loop over all cells and output data to file.
template<unsigned DIM>
void SimplePopulationExtentTracker<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{

  if(SimulationTime::Instance()->GetTimeStepsElapsed() % GetInterval() == 0){

    c_vector<double, DIM> minLocations;
    for(int i=0; i<DIM; i++){
      minLocations[i] = DBL_MAX;
    }
    c_vector<double, DIM> maxLocations;
    for(int i=0; i<DIM; i++){
      maxLocations[i] = -DBL_MAX;
    }

    for (typename AbstractCellPopulation<DIM,DIM>::Iterator cell_iter = rCellPopulation.Begin();
    cell_iter != rCellPopulation.End(); ++cell_iter)
    { 

      c_vector<double, DIM> cell_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
      
      for(int i=0; i<DIM; i++){
        if(cell_location[i] > maxLocations[i]){
          maxLocations[i] = cell_location[i]; 
        }
        if(cell_location[i] < minLocations[i]){
          minLocations[i] = cell_location[i];
        }
      }

    }

    //Write data
    double meanExtent = 0;
    *outputFile << SimulationTime::Instance()->GetTime() << "\t";
    for(int i=0; i<DIM; i++){
      *outputFile << maxLocations[i] - minLocations[i] << "\t";
      meanExtent += maxLocations[i] - minLocations[i];
    }
    *outputFile << meanExtent/DIM << "\n";

    outputFile->flush();
  }

  //If the simulation is finished, close the output file.
  if(SimulationTime::Instance()->IsFinished()){
    outputFile->close();
  }

}



//Output this class's parameters to a log file
template<unsigned DIM>
void SimplePopulationExtentTracker<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
  *rParamsFile << "\t\t\t<SampleExtentDataEveryXTimesteps>" << GetInterval() << "</SampleExtentDataEveryXTimesteps>\n";
  
  // Call method on direct parent class
  AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile); 
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class SimplePopulationExtentTracker<1>;
template class SimplePopulationExtentTracker<2>;
template class SimplePopulationExtentTracker<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SimplePopulationExtentTracker)
