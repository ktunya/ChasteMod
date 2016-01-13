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

#ifndef STOCHASTICCELLCYCLEMODEL_HPP_
#define STOCHASTICCELLCYCLEMODEL_HPP_

#include "AbstractCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"

/**
 * A cell cycle model with a fixed duration M and S phase and a stochastic duration
 * G1 and G2 phase.
 */

class StochasticCellCycle : public AbstractCellCycleModel
{

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archive the cell-cycle model and random number generator, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);

        SerializableSingleton<RandomNumberGenerator>* p_wrapper = RandomNumberGenerator::Instance()->GetSerializationWrapper();
        archive & p_wrapper;

        archive & mRangeG1;
        archive & mRangeG2;
        archive & mStochasticG1;
        archive & mStochasticG2;
    }


    /**
    * Range for G1 and G2 about the mean 
    */
    double mRangeG1;
    double mRangeG2;

    /**
    * The G1 and G2 durations with a stochastic contribution added in;
    * valid for this cell cycle only 
    */
    double mStochasticG1;
    double mStochasticG2;


public:


    /**
    * Default constructor. Most values are left as DOUBLE UNSET
    * until the setter methods are called.
    */
    StochasticCellCycle();

    StochasticCellCycle(double meanStemG1, double meanTransitG1, double S, double meanG2,
                        double M, double rangeG1, double rangeG2);

    /**
    * Virtual destructor
    */
    ~StochasticCellCycle();


    /**
    * Setters and getters for the new range parameters
    */
    void   SetRangeG1(double rangeG1);
    double GetRangeG1();
    void   SetRangeG2(double rangeG2);
    double GetRangeG2();
    
    /**
    * Setters for the G1 and G2 durations that include a stochastic contribution
    */
    void SetStochasticG1Duration();
    void SetStochasticG2Duration(); 


    /**
    * Set up ccm's at the start of the simulation and after a division
    */
    virtual void Initialise();  
    virtual void InitialiseDaughterCell();   
    /**
    * Reset parent cell after division triggered
    */
    virtual void ResetForDivision();


    /**
     * Overriden update cell cycle phase method
     */
    virtual void UpdateCellCyclePhase();

    /**
     * Overriden ready to divide method
     */
    virtual bool ReadyToDivide();


    /**
     * Overriden method to create a new ccm 
     */
    virtual AbstractCellCycleModel* CreateCellCycleModel();


    /**
     * Outputs cell cycle model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(StochasticCellCycle)

#endif /*STOCHASTICCELLCYCLEMODEL_HPP_*/