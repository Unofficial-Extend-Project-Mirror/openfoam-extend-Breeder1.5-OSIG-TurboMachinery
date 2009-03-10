/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008 Franz Blaim All rights reserved
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "primitiveMixingPlanePatch.H"
#include "addToRunTimeSelectionTable.H"
#include "SubField.H"
#include "entry.H"
#include "OFstream.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(primitiveMixingPlanePatch, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

primitiveMixingPlanePatch::primitiveMixingPlanePatch
(
   const List<face>& faces, 
   const pointField& points,
   const primitivePatch& OrigPatch,
   const point& CCSCenter,
   const vector& CCSaxis,
   const vector& CCSdirection,
   const word& name,
   bool blockMeshDict
)
: 
primitiveRotationalPatch
    (
        faces,
        points,
        OrigPatch,
        CCSCenter,
        CCSaxis,
        CCSdirection, 
        name
),
masterToPatchPtr_(NULL),
tensorListList_(0)
{   
       
       if(!blockMeshDict)
           calcRadialGeometry();
       
       Info << "In constructor" << endl;
}
   
   
// TODO The simple constructor must be implemented too

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

primitiveMixingPlanePatch::~primitiveMixingPlanePatch()
{
    //clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//bool primitiveMixingPlanePatch::constraintType(const word& pt)
//{
//    return pointPatchField<scalar>::PointPatchConstructorTablePtr_->found(pt);
//}

void primitiveMixingPlanePatch::movePoints(const pointField& p)
{
    primitiveAuxPatch::movePoints(p);
}

const List<face>& primitiveMixingPlanePatch::GetFaces() const
{
    return ( AllFaces_ );
}

void primitiveMixingPlanePatch::writePatches(const word& patchName)
{ 
  OFstream osNew("faces.obj");
           
  forAll(localPoints_, i)
  {
       writeOBJ(osNew, localPoints_[i]);
  }
  forAll(AllFaces_,i)
  {
       osNew << "f " ; 
       forAll(AllFaces_[i],j)
       {
          osNew << AllFaces_[i][j] + 1 << " " ;                                           
       }
       osNew << endl;
  }   
}


// * * * * * * * * * * * * * * *  protected Functions * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * *  template Specializations  * * * * * * * * * * * * * //
Field<vector> primitiveMixingPlanePatch::Average
(
     const List<vector>& sliceAvg
)
{
    Field<vector> returnField( 0 );

    // To adopt the calculated average to the local face direction 
    // this is especially necessary for the radial and mixed flow calculations 
    // a tensor rotatio must be applied for every local face
    // Therefore a List<tensor> pointer which stores the difference angle rotations is 
    // generated after the first averaging.
    vector AvgNormal(0,0,0);
    
    forAll( AllFaces_, i)
    {
        AvgNormal += AllFaces_[i].normal(localPoints_);
    }
    Info << AvgNormal << endl;
    AvgNormal /= AllFaces_.size();
    
    
    forAll(sliceFaces_,i)
    {      
      forAll(sliceFaces_[i], j)
      {   
           if( sliceAvg.size() <= i )
           { 
               FatalError
                      << "List sliceAvg size does not match with Rotationalpatch" << endl 
                      << " please correct geometry " << exit(FatalError);;   
           }
           returnField.setSize(returnField.size()+1);
           returnField[sliceFaces_[i][j]] = tensorListList_[i][j] & sliceAvg[i];//;
      }
   }
    
     
    
   return (returnField);
}

Field<vector> primitiveMixingPlanePatch::Average
(
    const Field<vector>& tpf,
    List<vector>& sliceAvg,
    scalarField massField
)const
{
    bool bdebug = false;
    Field<vector> returnField( tpf.size(), pTraits<vector>::zero );
       
    sliceAvg.clear();
    sliceAvg.setSize(0);
          
    vector avg =  pTraits<vector>::zero; 
    
    
    forAll(sliceFaces_,i)
    {   
              
          scalar Mass = 0;
          avg =  pTraits<vector>::zero;
          forAll(sliceFaces_[i],j)
          {
              // It is important to transform into the correct coordinate system 
              // it is the local cylindrical coordinate system
              //point faceCenter = AllFaces_[sliceFaces_[i][j]].centre(localPoints_) ;
             // Info << "FaceCentre"<< faceCenter << endl;           
              //tensor  tensorCartCCS = tensorFromCartesianToCCS( faceCenter );
                                                  
              avg = avg +    (tensorListList_[i][j].T() &  tpf[sliceFaces_[i][j]])
                      *AllFaces_[sliceFaces_[i][j]].mag(localPoints_)*massField[sliceFaces_[i][j]];  
              Mass += AllFaces_[sliceFaces_[i][j]].mag(localPoints_)*massField[sliceFaces_[i][j]];
                  
          }
           
   //       Info << avg << endl; 
          if (bdebug)
             Info << "Area " << Mass << endl;
              
         if(Mass == 0)
         { 
            FatalError
              << "Mass of RotationalPatch is equal zero" << " please correct geometry " << exit(FatalError);;   
         }
         else
         {
             
             
            avg = avg / ( Mass );
            sliceAvg.setSize(i+1);
            sliceAvg[i] = avg;
         }
              
         forAll(sliceFaces_[i], j)
         {
              returnField[sliceFaces_[i][j]] = avg;
         }
              
              
   }
   //   Info << "Return field " << returnField << endl;
      
   return (returnField);    
    
}
     
Field<vector> primitiveMixingPlanePatch::Average
(
   const Field<vector>& tpf,
   List<vector>& sliceAvg
)const    
{
    bool bdebug = false;
       Field<vector> returnField( tpf.size(), pTraits<vector>::zero );
             
       sliceAvg.clear();
       sliceAvg.setSize(0);
          
       vector avg =  pTraits<vector>::zero; 
       forAll(sliceFaces_,i)
       {   
              
          scalar Area = 0;
          avg =  pTraits<vector>::zero;
          forAll(sliceFaces_[i],j)
          {              
              //point faceCenter = AllFaces_[sliceFaces_[i][j]].centre(localPoints_);
             // Info << "FaceCentre" <<  faceCenter << endl;                          
              //tensor  tensorCartCCS = tensorFromCartesianToCCS( faceCenter );
                                                                
              avg = avg + (tensorListList_[i][j].T() & tpf[sliceFaces_[i][j]] ) /*cylindricalCS_.localVector()*/
                      *AllFaces_[sliceFaces_[i][j]].mag(localPoints_);// This is the area of one face
              
              Area += AllFaces_[sliceFaces_[i][j]].mag(localPoints_); // the whole area of the patch                  
          }
              
         if (bdebug)
             Info << "Area " << Area << endl;
      
         if(Area == 0)
         { 
            FatalError
              << "Area of RotationalPatch is equal zero" << " please correct geometry " << exit(FatalError);;   
         }
         else
         {
             
             
            avg = avg / ( Area );
            sliceAvg.setSize(i+1);
            sliceAvg[i] = avg;
         }
              
         forAll(sliceFaces_[i], j)
         {
              returnField[sliceFaces_[i][j]] = avg;
         }
              
              
      }
       // This function is called
   //  Info << "Return field " << returnField << endl;
      
      return (returnField);        
}



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "primitiveMixingPlaneCalcGeometry.C"

// ************************************************************************* //


