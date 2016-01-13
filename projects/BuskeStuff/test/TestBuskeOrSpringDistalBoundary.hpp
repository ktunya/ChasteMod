
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

#ifndef TESTBUSKEORSPRINGDISTALBOUNDARY_HPP_
#define TESTBUSKEORSPRINGDISTALBOUNDARY_HPP_

#include <cxxtest/TestSuite.h>

#include "AbstractCellBasedTestSuite.hpp"
#include "OffLatticeSimulation.hpp"
#include "NodeBasedCellPopulationWithBuskeUpdate.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StochasticCellCycle.hpp"
#include "StochasticDurationCellCycleModel.hpp"

#include "PositionBasedDifferentiationModifier.hpp"
#include "VolumeTrackingModifier.hpp"
#include "PlaneBasedCellKillerWithRecording.hpp"
#include "CellTrackingOutput.hpp"


#include "BuskeAdhesiveForce.hpp"
#include "BuskeElasticForce.hpp"
#include "BuskeCompressionForce.hpp"
#include "BuskeDistalPotentialBoundaryCondition.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "StandardDistalArmBoundaryCondition.hpp"


class TestBuskeOrSpringDistalBoundary : public AbstractCellBasedTestSuite
{
public:

    /*
     * Creates a simulation with a crypt/distal arm shaped boundary and a proliferative zone at one end.
     * Cell positions are tracked.
     * 
     * The test can be run either with Buske forces or with a Generalised Linear Spring force.
    */


    void TestBuskeOrSpringDistal() throw (Exception){
 
        EXIT_IF_PARALLEL;   

        std::cout << "FIXED SEED!" << std::endl;
        RandomNumberGenerator::Instance()->Reseed(0);

        //-----------------------------------------------------------------------------------
        //-----------------------------------------------------------------------------------

        // Force law settings

        enum Force {BUSKE, SPRING}; 
        Force ForceLaw = BUSKE;

        // buske only
        bool adhesionOn = true;
        bool elasticityOn = true;
        bool compressionOn = true;
        bool varyingRadiiOn = true;
        double eps = 200*pow(10,-6);         
        double d = (4/3)*pow(10,-3);         
        double k = 1000;                     
        double dampingIntercell = 5 * pow(10, 10);
        double dampingMedium = 3.2;
        double dampingVolume = 400; 

        // spring only
        double meinekeSpringConstant = 15;
        double adhesionDecay = 5;

        // both
        double meinekeDivisionSeparation = 0.95; 

        //-----------------------------------------------------------------------------------
        //-----------------------------------------------------------------------------------

        // Boundary condition properties

        enum Boundary {ARTIFICIAL, FORCEBASED}; 
        Boundary BoundaryType = FORCEBASED;

        double armLength = 248;
        double armRadius = 11.3;
        double cellRowSpacing = 5;
        double cellsPerRow = 10.0;
        double cellRadius = 2.8;
        double prolifZoneLength = 70;

        //-----------------------------------------------------------------------------------
        //-----------------------------------------------------------------------------------

        // Cell cycle properties
        
        double rangeG1 = 0.016;
        double rangeG2 = 0.312;
        double transitG1 = 0.16;
        double stemG1 = 0.16;
        double G2 = 3.12;
        double S = 4.56;
        double M = 0.16;

        //-----------------------------------------------------------------------------------
        //-----------------------------------------------------------------------------------

        // Create mesh and cell population

        NodesOnlyMesh<3>* mesh = CreateDistalArmCellLocations(armLength, armRadius, cellRowSpacing, cellsPerRow, cellRadius);

        MAKE_PTR(StemCellProliferativeType, pStemType);
        MAKE_PTR(DifferentiatedCellProliferativeType, pDiffType);

        int nCells = mesh->GetNumNodes();
        std::vector<CellPtr> cells;
        cells.clear();
        cells.reserve(nCells);

        for (unsigned i=0; i<nCells; i++)
        {
            // Add cell cycle
            StochasticCellCycle* ccm = new StochasticCellCycle(stemG1, transitG1, S, G2, M, rangeG1, rangeG2);
            ccm->SetDimension(3);

            // Make WT cell
            boost::shared_ptr<AbstractCellProperty> wtState(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
            CellPtr pCell(new Cell(wtState, ccm));
            
            // Make cells outside the proliferative zone differentiate
            if( (mesh->GetNode(i))->rGetLocation()[0] < prolifZoneLength ){
                pCell->SetCellProliferativeType(pStemType);
                pCell->GetCellData()->SetItem("Differentiated",0);
            }else{
                pCell->SetCellProliferativeType(pDiffType);
                pCell->GetCellData()->SetItem("Differentiated",1);
            }

            // Set random birth name 
            double birth_time = ccm->GetAverageStemCellCycleTime()*RandomNumberGenerator::Instance()->ranf();
            pCell->SetBirthTime(-birth_time);

            // Set cell data  
            pCell->GetCellData()->SetItem("Radius", cellRadius);
            pCell->GetCellData()->SetItem("RelaxedRadius", cellRadius);
            pCell->GetCellData()->SetItem("IsBuskeKnot", 0);
            pCell->GetCellData()->SetItem("volume", 0);

            cells.push_back(pCell);
        }

        NodeBasedCellPopulationWithBuskeUpdate<3> cellPopulation(*mesh, cells);

        //-----------------------------------------------------------------------------------
        //-----------------------------------------------------------------------------------

        // Set some cell population options

        // population options for both forces
        cellPopulation.SetUseVariableRadii( true );
        cellPopulation.SetMeinekeDivisionSeparation( meinekeDivisionSeparation );
        cellPopulation.SetAbsoluteMovementThreshold(2);

        // population options for Buske only
        if(ForceLaw == BUSKE){
            cellPopulation.SetUseVaryingRadii( varyingRadiiOn );
            cellPopulation.SetDampingConstantIntercell( dampingIntercell );  
            cellPopulation.SetDampingConstantMedium( dampingMedium );                 
            cellPopulation.SetDampingConstantVolume( dampingVolume );                 
            if( adhesionOn ){
                cellPopulation.EnableAdhesion(eps);
            }
            if( compressionOn ){
                cellPopulation.EnableCompression(k);
            }
            if( elasticityOn ){
                cellPopulation.EnableElasticity(d);
            }
        }

        //-----------------------------------------------------------------------------------
        //-----------------------------------------------------------------------------------

        // Set up the simulation
        
        OffLatticeSimulation<3> simulator(cellPopulation);
        
        std::stringstream filename;
        if(ForceLaw == BUSKE){
            filename << "TestBuskeDistal_";
            if(adhesionOn){
                filename << "eps" << eps << "_";
            }
            if(elasticityOn){
                filename << "d" << d << "_";
            }
            if(compressionOn){
                filename << "k" << k << "_";
            }
        }
        if(ForceLaw == SPRING){
            filename << "TestSpringDistal_mu" << meinekeSpringConstant << "_alpha" << adhesionDecay << "_";
        }
        simulator.SetOutputDirectory(filename.str());
        float stepsPerHour = 2500;
        simulator.SetDt(1/stepsPerHour); // Don't forget to try a run lowering dt, to check convergence
        simulator.SetSamplingTimestepMultiple((int)(stepsPerHour));
        simulator.SetEndTime(200);

        //-----------------------------------------------------------------------------------
        //-----------------------------------------------------------------------------------

        // Add a cell killer at the end of the tube

        c_vector<double, 3> point;
        point[0] = armLength;
        point[1] = 0;
        point[2] = 0;
        c_vector<double, 3> normal;
        normal[0] = 1;
        normal[1] = 0;
        normal[2] = 0;

        MAKE_PTR_ARGS(PlaneBasedCellKillerWithRecording<3>, pKiller, (&cellPopulation, point, normal, filename.str()));
        simulator.AddCellKiller(pKiller);     

        //-----------------------------------------------------------------------------------
        //-----------------------------------------------------------------------------------
        
        // Add a distal arm boundary condition

        if(BoundaryType == FORCEBASED){
            double forceStrength = 10000;
            MAKE_PTR_ARGS(BuskeDistalPotentialBoundaryCondition<3>, pBoundary, (armLength, armRadius, forceStrength));
            simulator.AddForce(pBoundary);
        }else if(BoundaryType == ARTIFICIAL){
            MAKE_PTR_ARGS(StandardDistalArmBoundaryCondition<3>, pBoundary, (&cellPopulation, armRadius, armLength));
            simulator.AddCellPopulationBoundaryCondition(pBoundary);
        }else{
        }

        //-----------------------------------------------------------------------------------
        //-----------------------------------------------------------------------------------

        // Add forces
        if(ForceLaw == BUSKE){
            
            if(compressionOn){
                MAKE_PTR(BuskeCompressionForce<3>, pCompression);
                pCompression->SetCompressionEnergyParameter(k);
                simulator.AddForce(pCompression);
            }
            if(elasticityOn){
                MAKE_PTR(BuskeElasticForce<3>, pElastic);
                pElastic->SetDeformationEnergyParameter(d);
                simulator.AddForce(pElastic);
            }
            if(adhesionOn){
                MAKE_PTR(BuskeAdhesiveForce<3>, pAdhesive);
                pAdhesive->SetAdhesionEnergyParameter(eps);
                simulator.AddForce(pAdhesive);
            }

        }else if(ForceLaw == SPRING){
            MAKE_PTR(GeneralisedLinearSpringForce<3>, pGLS);
            pGLS->SetAdhesionDecayAlpha(adhesionDecay);
            pGLS->SetMeinekeSpringStiffness(meinekeSpringConstant);
            simulator.AddForce(pGLS);
        }else{
        }

        //-----------------------------------------------------------------------------------
        //-----------------------------------------------------------------------------------

        // Add a differentiation modifier

        c_vector<double, 3> pointDiff;
        pointDiff[0] = prolifZoneLength;
        pointDiff[1] = 0;
        pointDiff[2] = 0;
        c_vector<double, 3> normalDiff;
        normalDiff[0] = 1;
        normalDiff[1] = 0;
        normalDiff[2] = 0;

        MAKE_PTR_ARGS(PositionBasedDifferentiationModifier<3>, pDiff, (pointDiff, normalDiff));
        simulator.AddSimulationModifier(pDiff);

        //-----------------------------------------------------------------------------------
        //-----------------------------------------------------------------------------------

        // Add cell position output
        
        MAKE_PTR_ARGS(CellTrackingOutput<3>, pTracking, ((int)stepsPerHour, 1));
        simulator.AddSimulationModifier(pTracking);
        
        // For compression comparisons?
        if(ForceLaw == SPRING){
            MAKE_PTR(VolumeTrackingModifier<3>, pVolumes);
            simulator.AddSimulationModifier(pVolumes);
        } 

        //-----------------------------------------------------------------------------------
        //-----------------------------------------------------------------------------------

        std::cout << "Solver started" << std::endl;
        try{
            simulator.Solve();
        }catch(int e){
            std::cout << "Error " << e << std::endl;
        }
    }





private:

    NodesOnlyMesh<3>* CreateDistalArmCellLocations(double armLength, 
                                                   double armRadius, 
                                                   double rowSpacing, 
                                                   int cellsPerRow, 
                                                   double cellRadius)
    {

        double angleSpacing = 2*M_PI / cellsPerRow;
        
        std::vector< Node<3>* > nodes;
        unsigned nodeIndex = 0;

        for(double x = -(armRadius-cellRadius); x < armLength; x += rowSpacing){
            for(double angle = 0; angle < 2*M_PI; angle += angleSpacing){
                Node<3>* n = new Node<3>(nodeIndex, false, x, (armRadius-cellRadius)*(cos(angle)), (armRadius-cellRadius)*(sin(angle)));
                n->SetRadius(cellRadius);
                nodes.push_back(n);
                nodeIndex++;
            }
        }

        NodesOnlyMesh<3>* distalMesh = new NodesOnlyMesh<3>;
        distalMesh->ConstructNodesWithoutMesh(nodes, 10);

        for (unsigned i = 0; i<nodes.size(); i++)
        {
           delete nodes[i];
        }

        return distalMesh;
    }
    
};

#endif /*TESTBUSKEORSPRINGDISTALBOUNDARY_HPP_*/
