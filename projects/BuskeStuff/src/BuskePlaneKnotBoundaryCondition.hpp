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

#ifndef BUSKEPLANEKNOTBOUNDARYCONDITION_HPP_
#define BUSKEPLANEKNOTBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
 * A plane cell population boundary condition class, which prevents nodes moving through
 * a specified plane in the domain. The plane is formed by a collection of tightly packed "knots"
 * (see ), which are immobile and exert a repulsion force on real cells. Designed for use with 
 * Buske forces, specifically requires BuskeKnotForce to function.
 */
template<unsigned DIM>
class BuskePlaneKnotBoundaryCondition : public AbstractCellPopulationBoundaryCondition<DIM>
{
private:

    /**
     * Bottom left corner of the boundary plane
     */
    c_vector<double, DIM> bottomLeft;

    /**
     * Top right corner of the boundary plane
     */
    c_vector<double, DIM> topRight;

    /**
     * The outward-facing unit normal vector to the boundary plane.
     */
    c_vector<double, DIM> normal;

    /**
     * Separation of the Buske knots that make up the boundary.
     */
    double knotSpacing;


    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<DIM> >(*this);
    }

public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population
     * @param bottomLeft bottom left of the boundary plane
     * @param topRight top right of the boundary plane
     * @param normal the outward-facing unit normal vector to the boundary plane
     * @param knotSpacing the separation between Buske knots
     */
    BuskePlaneKnotBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation,
                           c_vector<double, DIM> BottomLeft,
                           c_vector<double, DIM> TopRight,
                           c_vector<double, DIM> Normal,
                           double KnotSpacing = 1e-7);

    /*
    * Add a knot to the simulation cell population at position cellPos
    */
    void AddKnot(c_vector<double,DIM> cellPos);

    /**
     * @return #bottomLeft.
     */
    const c_vector<double, DIM>& rGetBottomLeftOfPlane() const;

    /**
     * @return #topRight.
     */
    const c_vector<double, DIM>& rGetTopRightOfPlane() const;

    /**
     * @return #normal.
     */
    const c_vector<double, DIM>& rGetNormalToPlane() const;

    /**
     * @return #knotSpacing.
     */
    const double& rGetKnotSpacing() const;

    /**
     * Overridden ImposeBoundaryCondition() method.
     *
     * Apply the cell population boundary conditions.
     *
     * @param rOldLocations the node locations before any boundary conditions are applied
     */
    void ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations);

    /**
     * Overridden VerifyBoundaryCondition() method.
     * Verify the boundary conditions have been applied.
     * This is called after ImposeBoundaryCondition() to ensure the condition is still satisfied.
     *
     * @return whether the boundary conditions are satisfied.
     */
    bool VerifyBoundaryCondition();

    /**
     * Overridden OutputCellPopulationBoundaryConditionParameters() method.
     * Output cell population boundary condition parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BuskePlaneKnotBoundaryCondition)



namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a BuskePlaneKnotBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const BuskePlaneKnotBoundaryCondition<DIM>* t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;

    // Archive c_vectors one component at a time
    c_vector<double, DIM> bl = t->rGetBottomLeftOfPlane();
    for (unsigned i=0; i<DIM; i++)
    {
        ar << bl[i];
    }
    c_vector<double, DIM> tr = t->rGetTopRightOfPlane();
    for (unsigned i=0; i<DIM; i++)
    {
        ar << tr[i];
    }
    c_vector<double, DIM> n = t->rGetNormalToPlane();
    for (unsigned i=0; i<DIM; i++)
    {
        ar << n[i];
    }

    ar << t->rGetKnotSpacing();
}

/**
 * De-serialize constructor parameters and initialize a BuskePlaneKnotBoundaryCondition.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, BuskePlaneKnotBoundaryCondition<DIM>* t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Archive c_vectors one component at a time
    c_vector<double, DIM> bl;
    for (unsigned i=0; i<DIM; i++)
    {
        ar >> bl[i];
    }
    c_vector<double, DIM> tr;
    for (unsigned i=0; i<DIM; i++)
    {
        ar >> tr[i];
    }
    c_vector<double, DIM> n;
    for (unsigned i=0; i<DIM; i++)
    {
        ar >> n[i];
    }

    double knotSpacing;
    ar >> knotSpacing;

    // Invoke inplace constructor to initialise instance
    ::new(t)BuskePlaneKnotBoundaryCondition<DIM>(p_cell_population, bl, tr, n, knotSpacing);
}
}
} // namespace ...

#endif /*BUSKEPLANEKNOTBOUNDARYCONDITION_HPP_*/
