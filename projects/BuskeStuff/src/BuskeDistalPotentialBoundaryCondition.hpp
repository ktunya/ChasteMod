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

#ifndef BUSKEDISTALPOTENTIALBOUNDARYCONDITION_HPP_
#define BUSKEDISTALPOTENTIALBOUNDARYCONDITION_HPP_

#include "AbstractForce.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>


/*
 * Adds a distal arm-like / colonic crypt-like boundary to a Buske force simulation by applying an appropriate energy potential.
 * The force on cells scales exponentially with distance from the boundary surface.
 */

template<unsigned DIM>
class BuskeDistalPotentialBoundaryCondition : public AbstractForce<DIM>
{
    friend class TestForcesNotForRelease;


private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & strength;
        archive & length;
        archive & radius;
    }

    /*
     * Strength of the exponential force (multiplying factor)
     */
    double strength;


    /*
     * Length of the boundary tube
     */
    double length;
    

    /**
     * Radius of the boundary tube
     */
    double radius;


public:

    /**
     * Constructor.
     */
    BuskeDistalPotentialBoundaryCondition(double inLength, double inRadius, double inStrength = 100);

    /**
     * Getter methods for parameters
     */
    const double GetLength() const;
    const double GetRadius() const;
    const double GetStrength() const;


    /*
    * Adds a force contribution to cells due to the boundary
    */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BuskeDistalPotentialBoundaryCondition)



namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a BuskeDistalPotentialBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const BuskeDistalPotentialBoundaryCondition<DIM>* t, const BOOST_PFTO unsigned int file_version)
{

    double len = t->GetLength();
    ar << len;

    double rad = t->GetRadius();
    ar << rad;

    double str = t->GetStrength();
    ar << str;
}

/**
 * De-serialize constructor parameters and initialize a BuskeDistalPotentialBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, BuskeDistalPotentialBoundaryCondition<DIM>* t, const unsigned int file_version)
{

    double len;
    ar >> len;

    double rad;
    ar >> rad;

    double str;
    ar >> str;

    // Invoke inplace constructor to initialise instance
    ::new(t)BuskeDistalPotentialBoundaryCondition<DIM>(len, rad, str);
}
}
} // namespace ...

#endif /*BUSKEDISTALPOTENTIALBOUNDARYCONDITION_HPP_*/
