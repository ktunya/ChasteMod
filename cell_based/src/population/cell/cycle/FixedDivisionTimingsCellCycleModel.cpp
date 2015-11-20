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

#include "FixedDivisionTimingsCellCycleModel.hpp"


FixedDivisionTimingsCellCycleModel::FixedDivisionTimingsCellCycleModel(double timeBetweenDivisions):
targetTimeBetweenDivisions(timeBetweenDivisions),
divideNextStep(false)
{
	std::cout << "WARNING: A FixedDivisionTimingsCellCycleModel manually divides cells every T hours, ";
	std::cout << "allowing variability in division timings to be controlled when comparing simulations. ";
	std::cout << "It does not update or use cell cycle phases, and cell cycle phase getters and setters will ";
	std::cout << "return irrelevant or unset values." << std::cout;
}	

FixedDivisionTimingsCellCycleModel::~FixedDivisionTimingsCellCycleModel()
{
}


bool FixedDivisionTimingsCellCycleModel::ReadyToDivide(){

	double currentTime = SimulationTime::Instance()->GetTime();
	double nextTime = currentTime + SimulationTime::Instance()->GetTimeStep();

	int cycleNumberCurrent = (int)(currentTime/timeBetweenDivisions);
	int cycleNumberNext = (int)(nextTime/timeBetweenDivisions);

	if(cycleNumberNext > cycleNumberCurrent){
		
		// A division should occur between this time step and the next.
		// Find out which is closer to the target time.

		double remainderCurrent = currentTime - cycleNumberCurrent*timeBetweenDivisions;
		double remainderNext = nextTime - cycleNumberNext*timeBetweenDivisions;

		if(remainderCurrent <= remainderNext){
			return true;
		}else{
			divideNextStep = true;
			return false;
		}
	}

	if(divideNextStep == true){
		return true;
	}

	return false;
};


void FixedDivisionTimingsCellCycleModel::UpdateCellCyclePhase(){
};


AbstractCellCycleModel* FixedDivisionTimingsCellCycleModel::CreateCellCycleModel(){
	return new FixedDivisionTimingsCellCycleModel();
};