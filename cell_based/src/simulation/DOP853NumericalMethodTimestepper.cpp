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

#include "DOP853NumericalMethodTimestepper.hpp"
#include "StepSizeException.hpp"
using boost::numeric::ublas::norm_2;


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>	
DOP853NumericalMethodTimestepper<ELEMENT_DIM,SPACE_DIM> :: DOP853NumericalMethodTimestepper( 
	                                                          AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>&  inputCellPopulation, 
                                                            std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM, SPACE_DIM> > >&  inputForceCollection)
:AbstractNumericalMethodTimestepper<ELEMENT_DIM,SPACE_DIM> ( inputCellPopulation, inputForceCollection )
{	
  //Set some default error tolerances
  absoluteErrorTolerance = 1e-5;
  relativeErrorTolerance = 1e-5;
};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>  
DOP853NumericalMethodTimestepper<ELEMENT_DIM,SPACE_DIM>::~DOP853NumericalMethodTimestepper(){
};



template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>  
void DOP853NumericalMethodTimestepper<ELEMENT_DIM,SPACE_DIM>::UpdateAllNodePositions(double dt){

  if(this->nonEulerSteppersEnabled){


    // Check the current step-size. If it's become unreasonably low, kill the simulation. 
    double currentTime = SimulationTime::Instance()->GetTime();
    if (dt <= currentTime * std::numeric_limits<double>::epsilon()) {
      EXCEPTION("Time step underflow");
    }

    int systemDimension = this->rCellPopulation.GetNumNodes() * SPACE_DIM;
    std::vector< c_vector<double,SPACE_DIM> > initialLocations = this->SaveCurrentLocations();             

    // Evaluate Ks----------------------------------------------------------------------------------------------------------------------
    // ---------------------------------------------------------------------------------------------------------------------------------

    std::vector< c_vector<double,SPACE_DIM> > K1  =  this->ComputeAndSaveForces(); 

    // SET TO initialLocations + dt*(a21*K1)--------------------------------------------------------------------------------------------
    int index = 0;
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rCellPopulation.rGetMesh().GetNodeIteratorBegin();
        node_iter != this->rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    {
      c_vector<double, SPACE_DIM> newLocation = initialLocations[index] + dt * a21*K1[index]; 
      ChastePoint<SPACE_DIM> new_point(newLocation);
      this->rCellPopulation.SetNode(node_iter->GetIndex(), new_point);
    }        
    std::vector< c_vector<double,SPACE_DIM> > K2  =  this->ComputeAndSaveForces();
    

    // SET TO initialLocations + dt*(a31*K1+a32*K2)-------------------------------------------------------------------------------------
    index = 0;
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rCellPopulation.rGetMesh().GetNodeIteratorBegin();
        node_iter != this->rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    {
      c_vector<double, SPACE_DIM> newLocation = initialLocations[index] + dt * (a31*K1[index] + a32*K2[index]);
      ChastePoint<SPACE_DIM> new_point(newLocation);
      this->rCellPopulation.SetNode(node_iter->GetIndex(), new_point);
    }  
    std::vector< c_vector<double,SPACE_DIM> > K3  =  this->ComputeAndSaveForces();


    // SET TO initialLocations + dt*(a41*K1+a43*K3)--------------------------------------------------------------------------------------
    index = 0;
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rCellPopulation.rGetMesh().GetNodeIteratorBegin();
        node_iter != this->rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    {
      c_vector<double, SPACE_DIM> newLocation = initialLocations[index] + dt * (a41*K1[index] + a43*K3[index]);
      ChastePoint<SPACE_DIM> new_point(newLocation);
      this->rCellPopulation.SetNode(node_iter->GetIndex(), new_point);
    }  
    std::vector< c_vector<double,SPACE_DIM> > K4  =  this->ComputeAndSaveForces();


    // SET TO initialLocations + dt*(a51*K1+a53*K3+a54*K4)------------------------------------------------------------------------------
    index = 0;
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rCellPopulation.rGetMesh().GetNodeIteratorBegin();
        node_iter != this->rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    {
      c_vector<double, SPACE_DIM> newLocation = initialLocations[index] + dt * (a51*K1[index] + a53*K3[index] + a54*K4[index]);
      ChastePoint<SPACE_DIM> new_point(newLocation);
      this->rCellPopulation.SetNode(node_iter->GetIndex(), new_point);
    }  
    std::vector< c_vector<double,SPACE_DIM> > K5  =  this->ComputeAndSaveForces();


    // SET TO initialLocations + dt*(a61*K1+a64*K4+a65*K5)------------------------------------------------------------------------------
    index = 0;
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rCellPopulation.rGetMesh().GetNodeIteratorBegin();
        node_iter != this->rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    {
      c_vector<double, SPACE_DIM> newLocation = initialLocations[index] + dt * (a61*K1[index] + a64*K4[index] + a65*K5[index]);
      ChastePoint<SPACE_DIM> new_point(newLocation);
      this->rCellPopulation.SetNode(node_iter->GetIndex(), new_point);
    } 
    std::vector< c_vector<double,SPACE_DIM> > K6  =  this->ComputeAndSaveForces();


    // SET TO initialLocations + dt*(a71*K1+a74*K4+a75*K5+a76*K6)----------------------------------------------------------------------
    index = 0;
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rCellPopulation.rGetMesh().GetNodeIteratorBegin();
        node_iter != this->rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    {
      c_vector<double, SPACE_DIM> newLocation = initialLocations[index] + dt * (a71*K1[index] + a74*K4[index] + a75*K5[index] + a76*K6[index]);
      ChastePoint<SPACE_DIM> new_point(newLocation);
      this->rCellPopulation.SetNode(node_iter->GetIndex(), new_point);
    } 
    std::vector< c_vector<double,SPACE_DIM> > K7  =  this->ComputeAndSaveForces();


    // SET TO initialLocations + dt*(a81*K1+a84*K4+a85*K5+a86*K6+a87*K7)---------------------------------------------------------------
    index = 0;
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rCellPopulation.rGetMesh().GetNodeIteratorBegin();
        node_iter != this->rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    {
      c_vector<double, SPACE_DIM> newLocation = initialLocations[index] + dt * (a81*K1[index] + a84*K4[index] + a85*K5[index] + a86*K6[index] + a87*K7[index]);
      ChastePoint<SPACE_DIM> new_point(newLocation);
      this->rCellPopulation.SetNode(node_iter->GetIndex(), new_point);
    } 
    std::vector< c_vector<double,SPACE_DIM> > K8  =  this->ComputeAndSaveForces();


    // SET TO initialLocations + dt*(a91*K1+a94*K4+a95*K5+a96*K6+a97*K7+a98*K8)--------------------------------------------------------
    index = 0;
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rCellPopulation.rGetMesh().GetNodeIteratorBegin();
        node_iter != this->rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    {
      c_vector<double, SPACE_DIM> newLocation = initialLocations[index] + dt * (a91*K1[index] + a94*K4[index] + a95*K5[index] + a96*K6[index] + a97*K7[index] + a98*K8[index]);
      ChastePoint<SPACE_DIM> new_point(newLocation);
      this->rCellPopulation.SetNode(node_iter->GetIndex(), new_point);
    } 
    std::vector< c_vector<double,SPACE_DIM> > K9  =  this->ComputeAndSaveForces();


    // SET TO initialLocations + dt*(a101*K1+a104*K4+a105*K5+a106*K6+a107*K7+a108*K8+a109*K9)-----------------------------------------
    index = 0;
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rCellPopulation.rGetMesh().GetNodeIteratorBegin();
        node_iter != this->rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    {
      c_vector<double, SPACE_DIM> newLocation = initialLocations[index] + dt * (a101*K1[index] + a104*K4[index] + a105*K5[index] + a106*K6[index] + a107*K7[index] + a108*K8[index] + a109*K9[index]);
      ChastePoint<SPACE_DIM> new_point(newLocation);
      this->rCellPopulation.SetNode(node_iter->GetIndex(), new_point);
    } 
    std::vector< c_vector<double,SPACE_DIM> > K10 =  this->ComputeAndSaveForces();
    

    // SET TO initialLocations + dt*(a111*K1+a114*K4+a115*K5+a116*K6+a117*K7+a118*K8+a119*K9+a1110*K10)-------------------------------
    index = 0;
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rCellPopulation.rGetMesh().GetNodeIteratorBegin();
        node_iter != this->rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    {
      c_vector<double, SPACE_DIM> newLocation = initialLocations[index] + dt * (a111*K1[index] + a114*K4[index] + a115*K5[index] + a116*K6[index] + a117*K7[index] + a118*K8[index] + a119*K9[index] + a1110*K10[index]);
      ChastePoint<SPACE_DIM> new_point(newLocation);
      this->rCellPopulation.SetNode(node_iter->GetIndex(), new_point);
    } 
    K2  =  this->ComputeAndSaveForces();


    // SET TO initialLocations + dt*(a121*K1+a124*K4+a125*K5+a126*K6+a127*K7+a128*K8+a129*K9+a1210*K10+a1211*K2)---------------------
    index = 0;
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rCellPopulation.rGetMesh().GetNodeIteratorBegin();
        node_iter != this->rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    {
      c_vector<double, SPACE_DIM> newLocation = initialLocations[index] + dt * (a121*K1[index] + a124*K4[index] + a125*K5[index] + a126*K6[index] + a127*K7[index] + a128*K8[index] + a129*K9[index] + a1210*K10[index] + a1211*K2[index]);
      ChastePoint<SPACE_DIM> new_point(newLocation);
      this->rCellPopulation.SetNode(node_iter->GetIndex(), new_point);
    } 
    K3  =  this->ComputeAndSaveForces();

    // Make an error estimate:--------------------------------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------------------------

    std::vector< c_vector<double,SPACE_DIM> > delta;
    delta.reserve(K1.size());
    std::vector< c_vector<double,SPACE_DIM> > err5thOrder;
    err5thOrder.reserve(K1.size());
    std::vector< c_vector<double,SPACE_DIM> > err3rdOrder;
    err3rdOrder.reserve(K1.size());

    for(int i=0; i<K1.size(); i++){
      delta.push_back( dt*(b1*K1[i] + b6*K6[i] + b7*K7[i] + b8*K8[i] + b9*K9[i] + b10*K10[i] + b11*K2[i] + b12*K3[i]) );
      err5thOrder.push_back( delta[i] - dt*(bhh1*K1[i] + bhh2*K9[i] + bhh3*K3[i]) );
      err3rdOrder.push_back( dt*(er1*K1[i] + er6*K6[i] + er7*K7[i] + er8*K8[i] + er9*K9[i] + er10*K10[i] + er11*K2[i] + er12*K3[i]) );
    }

    double err3 = 0.0;
    double err5 = 0.0;
    double maxDisplacement = 0.0;
    index = 0;
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rCellPopulation.rGetMesh().GetNodeIteratorBegin();
        node_iter != this->rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
    {
      double scaling = absoluteErrorTolerance + relativeErrorTolerance * std::max(norm_2(node_iter->rGetLocation()),  norm_2(node_iter->rGetLocation() + delta[index]));
      err3 = std::max(err3, norm_2(err3rdOrder[index])/scaling);
      err5 = std::max(err5, norm_2(err5thOrder[index])/scaling);
      maxDisplacement = std::max(norm_2(delta[index]), maxDisplacement);
    }

    double denominator = err3 * err3  +  0.01 * err5 * err5;
    if(denominator <= 0.0){
      denominator = 1.0;
    } 

    double MaxNormError = err5 * err5 / std::sqrt(denominator);

    // If the error has become too large, reduce dt; otherwise complete the step:-----------------------------------------------------
    //--------------------------------------------------------------------------------------------------------------------------------

    if (MaxNormError > 1.0) {

      double newStepSize = 0.5 * dt;
      std::string message("Error estimate exceeded the acceptable tolerance"); 
      bool isTerminalIfNonAdaptive = true;
      throw new StepSizeException(maxDisplacement, newStepSize, message, isTerminalIfNonAdaptive);

    } else {

      index = 0;
      for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->rCellPopulation.rGetMesh().GetNodeIteratorBegin();
          node_iter != this->rCellPopulation.rGetMesh().GetNodeIteratorEnd(); ++node_iter, ++index)
      {
        c_vector<double, SPACE_DIM> newLocation = initialLocations[index] + delta[index];
        ChastePoint<SPACE_DIM> new_point(newLocation);
        this->rCellPopulation.SetNode(node_iter->GetIndex(), new_point);
      }
      K1  =  this->ComputeAndSaveForces();

    }

  }else{

    this->ComputeAndSaveForces();
    this->rCellPopulation.UpdateNodeLocations(dt);    
  }

};


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>  
void DOP853NumericalMethodTimestepper<ELEMENT_DIM,SPACE_DIM>::SetAbsoluteErrorTolerance(double tol){
  absoluteErrorTolerance = tol;
};

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>  
void DOP853NumericalMethodTimestepper<ELEMENT_DIM,SPACE_DIM>::SetRelativeErrorTolerance(double tol){
  relativeErrorTolerance = tol;
};


///////// Explicit instantiation
template class DOP853NumericalMethodTimestepper<1,1>;
template class DOP853NumericalMethodTimestepper<1,2>;
template class DOP853NumericalMethodTimestepper<2,2>;
template class DOP853NumericalMethodTimestepper<1,3>;
template class DOP853NumericalMethodTimestepper<2,3>;
template class DOP853NumericalMethodTimestepper<3,3>;