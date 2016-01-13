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

#include "StochasticCellCycle.hpp"
#include "Exception.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"


StochasticCellCycle::StochasticCellCycle(double meanStemG1, double meanTransitG1, double S,
                                         double meanG2, double M, double rangeG1, double rangeG2)
    : AbstractCellCycleModel()
{  
    // Overwrite the unhelpful default values from the AbstractCellCycleModel.
    // Those values are only useful for colonic crypt simulations and not for other
    // types of biological model.
    SetStemCellG1Duration(meanStemG1);
    SetTransitCellG1Duration(meanTransitG1);
    SetG2Duration(meanG2); 
    SetSDuration(S);           
    SetMDuration(M);    
    mRangeG1 = rangeG1;
    mRangeG2 = rangeG2; 
    mStochasticG1 = DOUBLE_UNSET;
    mStochasticG2 = DOUBLE_UNSET;
}

StochasticCellCycle::StochasticCellCycle(): AbstractCellCycleModel(){ 
}

StochasticCellCycle::~StochasticCellCycle(){};




void StochasticCellCycle::SetRangeG1(double rangeG1){
    mRangeG1 = rangeG1;
};

double StochasticCellCycle::GetRangeG1(){
    return mRangeG1;
};

void StochasticCellCycle::SetRangeG2(double rangeG2){
    mRangeG2 = rangeG2;
};

double StochasticCellCycle::GetRangeG2(){
    return mRangeG2;
};




void StochasticCellCycle::SetStochasticG1Duration(){
    assert(mRangeG1 != DOUBLE_UNSET);
    // Generate a stochastic G1 duration for this cell
    if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
    {
        assert(mStemCellG1Duration != DOUBLE_UNSET);
        mStochasticG1 = mStemCellG1Duration + mRangeG1*(RandomNumberGenerator::Instance()->ranf() - 0.5);
    }
    else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
    {
        assert(mTransitCellG1Duration != DOUBLE_UNSET);
        mStochasticG1 = mTransitCellG1Duration + mRangeG1*(RandomNumberGenerator::Instance()->ranf() - 0.5);
    }
    else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mStochasticG1 = 100000;
    }
    else
    {
        NEVER_REACHED;
    }
};

void StochasticCellCycle::SetStochasticG2Duration(){
    assert(mRangeG2 != DOUBLE_UNSET);
    assert(mG2Duration != DOUBLE_UNSET);
    // Generate a stochastic G2 duration for this cell
    mStochasticG2 = mG2Duration + mRangeG2 * (RandomNumberGenerator::Instance()->ranf() - 0.5);
}; 


void StochasticCellCycle::Initialise(){
    SetStochasticG1Duration();
    SetStochasticG2Duration();
}

void StochasticCellCycle::InitialiseDaughterCell(){   
    this->AbstractCellCycleModel::InitialiseDaughterCell();
    SetStochasticG1Duration();
    SetStochasticG2Duration();
};

void StochasticCellCycle::ResetForDivision(){
    this->AbstractCellCycleModel::ResetForDivision();
    mBirthTime = SimulationTime::Instance()->GetTime();
    SetStochasticG1Duration();
    SetStochasticG2Duration();
}




void StochasticCellCycle::UpdateCellCyclePhase(){

    double time_since_birth = GetAge();
    assert(time_since_birth >= 0 );

    assert(mStochasticG1 != DOUBLE_UNSET);
    assert(mStochasticG2 != DOUBLE_UNSET);
    assert(mSDuration != DOUBLE_UNSET);
    assert(mMDuration != DOUBLE_UNSET); 

    if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mCurrentCellCyclePhase = G_ZERO_PHASE;
    }
    else if (time_since_birth < GetMDuration())
    {
        mCurrentCellCyclePhase = M_PHASE;
    }
    else if (time_since_birth < GetMDuration() + mStochasticG1)
    {
        mCurrentCellCyclePhase = G_ONE_PHASE;
    }
    else if (time_since_birth < GetMDuration() + mStochasticG1 + GetSDuration())
    {
        mCurrentCellCyclePhase = S_PHASE;
    }
    else if (time_since_birth < GetMDuration() + mStochasticG1 + GetSDuration() + mStochasticG2)
    {
        mCurrentCellCyclePhase = G_TWO_PHASE;
    }
};


bool StochasticCellCycle::ReadyToDivide(){
    assert(mpCell != NULL);
    if (!mReadyToDivide)
    {
        UpdateCellCyclePhase();
        if ( (mCurrentCellCyclePhase != G_ZERO_PHASE) &&
             (GetAge() >= GetMDuration() + mStochasticG1 + GetSDuration() + mStochasticG2) )
        {
            mReadyToDivide = true;
        }
    }
    return mReadyToDivide;
};





AbstractCellCycleModel* StochasticCellCycle::CreateCellCycleModel(){
    
    StochasticCellCycle* p_model = new StochasticCellCycle(GetStemCellG1Duration(),
                                                           GetTransitCellG1Duration(),
                                                           GetSDuration(),
                                                           GetG2Duration(),
                                                           GetMDuration(),
                                                           GetRangeG1(),
                                                           GetRangeG2());

    return p_model;
}




void StochasticCellCycle::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<RangeG1>" << GetRangeG1() << "</RangeG1>\n";
    *rParamsFile << "\t\t\t<RangeG2>" << GetRangeG2() << "</RangeG2>\n";
    AbstractCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(StochasticCellCycle)
