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
#include "PrimitivePatch.H"
#include "face.H"
#include "pointField.H"
#include "List.H"
#include "mathematicalConstants.H"

#include "OFstream.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Public Member Templates * * * * * * * * * * * * * //



template<class Type, class Patch> 
tmp<Field<Type> > primitiveMixingPlanePatch::MakeCircumferentialAverage
(
       const Field<Type>& masterField, 
       const Patch& masterPatch,
       List<Type>& sliceAvg,
       const vectorField& masterFieldMass         
) 
{
      tmp<Field<Type> > tFirst
      (
           new Field<Type>
           (
                AllFaces_.size(),
                pTraits<Type>::zero
           )
       );
      
         
      Field<Type>&  firstField = tFirst();

       // TODO FB Must be called otherwise the interpolation does not work correctly
//       pointField points = localPoints_;
      
       if(!primPatch_)
           primPatch_ = new primitiveAuxPatch(AllFaces_, localPoints_ );
            
       
       // TODO FB This should be removed to one central place because right now 
       // the masterToPatchPtr is only built once but it is required that 
       // a argument with a patch has to be given every time 

       
       if(mType_ != Axial)
           Info << "Not axial " << endl;
       if(!masterToPatchPtr_)
       {
            masterToPatchPtr_ = new ggiMMInterpolation
          (
                  *primPatch_,
                   masterPatch,         
                   tensorField(0),
                   tensorField(0),
                   vectorField(0)
           );

       }
       
       const vectorField& faceNormals = primPatch_->faceNormals(); 

       forAll(faceNormals,i)
       {
           faceNormals[i];        
       }
       
       // To get the flow in normal direction the masterFieldMass vectorField is mutiplied with
       // the faceNormals
       const vectorField& masterNormals = masterPatch.faceNormals();
       
       scalarField multipliedField = masterFieldMass & masterNormals;
                 
                 
      // the R() is the local to global Transformation tensor 
      // its transpond gives the global to local transformation tensor which is needed in 
      // that case
       
       firstField = masterToPatchPtr_->slaveToMaster(masterField);
      
                
       // Transfomation into the cylindrical coordinate system before averaging 
       // therefore the direction is not averaged in a cartesian coordinate frame 
       // which would lead to a wrong direction in the whole averaged field.
       
       scalarField massFieldInterpolated = masterToPatchPtr_->slaveToMaster(multipliedField); 
       
       
       firstField = Average(firstField, sliceAvg, massFieldInterpolated);

       
       tmp<Field<Type> > tResult
       (  
           new Field<Type>
           (
               masterPatch.size(),
               pTraits<Type>::zero
           )
       );
         
       Field<Type>& resultField = tResult();
       
       resultField = masterToPatchPtr_->masterToSlave(firstField);
       

       return (tResult);
}

template<class Type, class Patch> 
tmp<Field<Type> > primitiveMixingPlanePatch::MakeCircumferentialAverage
(
       const Field<Type>& masterField, 
       const Patch& masterPatch,
       List<Type>& sliceAvg  
) 
{
      tmp<Field<Type> > tFirst
      (
           new Field<Type>
           (
                AllFaces_.size(),
                pTraits<Type>::zero
           )
       );
      
         
      Field<Type>&  firstField = tFirst();

       // TODO FB Must be called otherwise the interpolation does not work correctly
//       pointField points = localPoints_;
      
      if(!primPatch_)
      {
             Info << "Make aux pointer "<< endl;
             primPatch_ = new primitiveAuxPatch(AllFaces_, localPoints_ );                    
      }
       
       if(!masterToPatchPtr_)
       {
  
               Info << "masterToPatchPtr_ " << endl;
               masterToPatchPtr_ = new ggiMMInterpolation
               (
                  *primPatch_,
                  masterPatch,         
                   tensorField(0),
                   tensorField(0),
                   vectorField(0)
               );
         }
       
       const vectorField& faceNormals = primPatch_->faceNormals(); 

       forAll(faceNormals,i)
       {
           faceNormals[i];        
       }
//       Info << "FaceNormals " << masterPatch.faceNormals() <<  faceNormals << endl;
//       Info << "Before Interpolation" << masterField << endl;
//       
       firstField = masterToPatchPtr_->slaveToMaster(masterField);
       
                           
       
//       Info <<  "After interPolation" << firstField << endl;
       firstField = Average(firstField,sliceAvg);

       
       tmp<Field<Type> > tResult
       (  
           new Field<Type>
           (
               masterPatch.size(),
               pTraits<Type>::zero
           )
       );
         
       Field<Type>& resultField = tResult();
       
       resultField = masterToPatchPtr_->masterToSlave(firstField);
       
//       Info << sliceAvg << endl;

       return (tResult);
}

template<class Type, class Patch> 
tmp<Field<Type> > primitiveMixingPlanePatch::MakeCircumferentialAverage
(
       const Patch& masterPatch,
       const List<Type>& sliceAvg  
) 
{
      tmp<Field<Type> > tFirst
      (
           new Field<Type>
           (
                AllFaces_.size(),
                pTraits<Type>::zero
           )
       );
      
         
      Field<Type>&  firstField = tFirst();
      
       // TODO FB Must be called otherwise the interpolation does not work correctly

      if(!primPatch_)
      {
          Info << "Make aux pointer "<< endl;
         primPatch_ = new primitiveAuxPatch(AllFaces_, localPoints_ );
                
      }
            
            
      
       if(!masterToPatchPtr_)
       {
  
               masterToPatchPtr_ = new ggiMMInterpolation
               (
                   *primPatch_,
                   masterPatch,  
                   tensorField(0),
                   tensorField(0),
                   vectorField(0)
                   // Contact_Sphere checks for a intersection on a sphere
                   // delivers very good results for a cylindrical surface
                   // but for the axial case it delivers bad results
                   // because it can happen that it points into the same direction as another 
                  
                );  
       }
       
       
       const vectorField& faceNormals = primPatch_->faceNormals(); 
//
       forAll(faceNormals,i)
       {
           faceNormals[i];        
       }
      
       firstField = Average(sliceAvg);

       
       tmp<Field<Type> > tResult
       (  
           new Field<Type>
           (
               masterPatch.size(),
               pTraits<Type>::zero
           )
       );
         
       Field<Type>& resultField = tResult();
       Info << "AllFaces " << AllFaces_.size() << endl;
       Info << firstField.size() << "slave size " <<  (*primPatch_).size() << endl;
       
       resultField = masterToPatchPtr_->masterToSlave(firstField);
   
       return (tResult);
}

// Starting of Average template functions
template<class Type > 
Field<Type> primitiveMixingPlanePatch::Average
(
       const Field<Type>& tpf
)const 
{
    bool bdebug = false;
    Field<Type> returnField( tpf.size(), pTraits<Type>::zero );
    
    Type avg =  pTraits<Type>::zero; 
    forAll(sliceFaces_,i)
    {   
        
        scalar Area = 0;
        avg =  pTraits<Type>::zero;
        forAll(sliceFaces_[i],j)
        {
            
            avg = avg + tpf[sliceFaces_[i][j]]
                   *AllFaces_[sliceFaces_[i][j]].mag(localPoints_);  
            Area += AllFaces_[sliceFaces_[i][j]].mag(localPoints_);
            
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
        }
        
        forAll(sliceFaces_[i], j)
        {
            returnField[sliceFaces_[i][j]] = avg;
        }
        
        
    }
//    Info << "Return field " << returnField << endl;
    
    return (returnField);
}


template<class Type > 
Field<Type> primitiveMixingPlanePatch::Average
(
          const Field<Type>& tpf,
          List<Type>& sliceAvg,
          scalarField massField
)const 
{
    bool bdebug = false;
    Field<Type> returnField( tpf.size(), pTraits<Type>::zero );
    
    sliceAvg.clear();
    sliceAvg.setSize(0);
       
    Type avg =  pTraits<Type>::zero; 
    forAll(sliceFaces_,i)
    {   
           
       scalar Mass = 0;
       avg =  pTraits<Type>::zero;
       forAll(sliceFaces_[i],j)
       {
           avg = avg + tpf[sliceFaces_[i][j]]
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


template<class Type > 
Field<Type> primitiveMixingPlanePatch::Average
(
          const Field<Type>& tpf,
          List<Type>& sliceAvg
)const 
{
    bool bdebug = false;
    Field<Type> returnField( tpf.size(), pTraits<Type>::zero );
    
//    Info << tpf << endl;
    
    sliceAvg.clear();
    sliceAvg.setSize(0);
       
    Type avg =  pTraits<Type>::zero; 
    forAll(sliceFaces_,i)
    {   
           
       scalar Area = 0;
       avg =  pTraits<Type>::zero;
       forAll(sliceFaces_[i],j)
       {
           avg = avg + tpf[sliceFaces_[i][j]]
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

template<class Type > 
Field<Type> primitiveMixingPlanePatch::Average
(
     const List<Type>& sliceAvg
)const
{
    Field<Type> returnField( 0 );

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
           returnField[sliceFaces_[i][j]] = sliceAvg[i];
      }
   }
    
//    Info << "returnField in Average " << returnField << endl; 
    
   return (returnField);
}

template<class Type, class Source, class Target> 
tmp<Field<Type> > primitiveMixingPlanePatch::InterpolateBetweenPatches
(
       const Field<Type>& tpf,
       const Source& sourcePatch,
       const Target& targetPatch
) const
{

    // TODO FB Must be called otherwise the interpolation does not work correctly
    typedef GGIInterpolation
    <
        Source,
        Target
    >    ggiMMInterpolation;

    
    ggiMMInterpolation*  patchToPatchPtr(NULL);
    
   // TODO Warum keine Abfrage ob der Pointer besteht ?
    patchToPatchPtr = new ggiMMInterpolation
    (
              sourcePatch,
              targetPatch,         
            tensorField(targetPatch.size(),
		   tensor(1,0,0,0,1,0,0,0,1) ),
                   tensorField( sourcePatch.size(),
		   tensor(1,0,0,0,1,0,0,0,1) )
    );
   
    
    return ( patchToPatchPtr->slaveToMaster(tpf) );
}





// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
