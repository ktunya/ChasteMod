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

#ifndef BUSKEPLANEPOTENTIALBOUNDARYCONDITION_HPP_
#define BUSKEPLANEPOTENTIALBOUNDARYCONDITION_HPP_

#include "AbstractForce.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>


/**
 * Adds a plane boundary to a Buske force simulation by applying an appropriate energy potential.
 * The force law chosen matches doi:10.1371/journal.pcbi.1001045 but does not use discrete knots.
 */

template<unsigned DIM>
class BuskePlanePotentialBoundaryCondition : public AbstractForce<DIM>
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
        archive & interactionEnergy;
        archive & thresholdAdhesionRatio;
    }
    
    /*
    * A point on the boundary plane
    */
    c_vector<double,DIM> point;

    /*
    * The outward facing normal to the boundary plane
    */
    c_vector<double,DIM> normal;

    /**
     * Interaction energy epsilon associated with this plane
     */
    double interactionEnergy;
    
    /**
     * Defines the distance from the plane at which adhesion is overtaken by repulsion
     */
    double thresholdAdhesionRatio;


public:

    /**
     * Constructor.
     */
    BuskePlanePotentialBoundaryCondition(c_vector<double,DIM> pointOnPlane, c_vector<double,DIM> normalToPlane, double interactionE = 1e-11, double threshAdhesionRatio = 0.95);

    /**
     * Getter methods for parameters
     */
    const double GetInteractionEnergy() const;
    const double GetThresholdAdhesionRatio() const;
    const c_vector<double,DIM> GetPoint() const;
    const c_vector<double,DIM> GetNormal() const;

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
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BuskePlanePotentialBoundaryCondition)



namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a BuskePlanePotentialBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const BuskePlanePotentialBoundaryCondition<DIM>* t, const BOOST_PFTO unsigned int file_version)
{

    c_vector<double,DIM> point = t->GetPoint();
    for(int i=0; i<DIM; i++){
        ar << point[i];
    }

    c_vector<double,DIM> norm = t->GetNormal();
    for(int i=0; i<DIM; i++){
        ar << norm[i];
    }

    double interact = t->GetInteractionEnergy();
    ar << interact;

    double threshAdhesion = t->GetThresholdAdhesionRatio();
    ar << threshAdhesion;
}

/**
 * De-serialize constructor parameters and initialize a BuskePlanePotentialBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, BuskePlanePotentialBoundaryCondition<DIM>* t, const unsigned int file_version)
{
    
    c_vector<double,DIM> point;
    for(int i=0; i<DIM; i++){
        ar >> point[i];
    }

    c_vector<double,DIM> norm;
    for(int i=0; i<DIM; i++){
        ar >> norm[i];
    }

    double interact;
    ar >> interact;

    double threshAdhesion;
    ar >> threshAdhesion;

    // Invoke inplace constructor to initialise instance
    ::new(t)BuskePlanePotentialBoundaryCondition<DIM>(point, norm, interact, threshAdhesion);
}
}
} // namespace ...

#endif /*BUSKEPLANEPOTENTIALBOUNDARYCONDITION_HPP_*/
