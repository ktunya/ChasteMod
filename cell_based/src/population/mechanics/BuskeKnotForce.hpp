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

#ifndef BuskeKnotForce_HPP_
#define BuskeKnotForce_HPP_

#include "AbstractForce.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * Adds a Buske force exerted by Buske knots (see doi:10.1371/journal.pcbi.1001045).
 * For use in conjunction with a Buske type boundary condition.
 */


template<unsigned DIM>
class BuskeKnotForce : public AbstractForce<DIM>
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
        archive & maxInteractionEnergy;
        archive & thresholdAdhesionRatio;
    }

    /**
     * Max interaction energy of a knot
     */
    double maxInteractionEnergy;
    
    /**
     * Defines the distance at which adhesion is overtaken by repulsion
     */
    double thresholdAdhesionRatio;


public:

    /**
     * Constructor.
     */
    BuskeKnotForce(double interaction = 1e-11, double omegaRatio = 0.95);

    /**
     * Setter/Getter methods for parameters
     */
    void SetMaxInteractionEnergy(double interaction);
    void SetThresholdAdhesionRatio(double omegaRatio);
    const double GetMaxInteractionEnergy() const;
    const double GetThresholdAdhesionRatio() const;

    
    /*
    * Adds force contributions due to the knots in the cell population
    */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);


    /**
     * Overridden OutputForceParameters() method.
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BuskeKnotForce)



namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a BuskeKnotForce.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const BuskeKnotForce<DIM>* t, const BOOST_PFTO unsigned int file_version)
{

    double maxInteract = t->GetMaxInteractionEnergy();
    ar << maxInteract;

    double threshAdhesion = t->GetThresholdAdhesionRatio();
    ar << threshAdhesion;
}

/**
 * De-serialize constructor parameters and initialize a BuskeKnotForce.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, BuskeKnotForce<DIM>* t, const unsigned int file_version)
{

    double maxInteract;
    ar >> maxInteract;

    double threshAdhesion;
    ar >> threshAdhesion;

    // Invoke inplace constructor to initialise instance
    ::new(t)BuskeKnotForce<DIM>(maxInteract, threshAdhesion);
}
}
} // namespace ...

#endif /*BuskeKnotForce_HPP_*/
