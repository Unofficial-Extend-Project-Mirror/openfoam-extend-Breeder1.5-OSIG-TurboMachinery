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


#include "primitiveDomainScalingPatch.H"
#include "PrimitivePatch.H"
#include "face.H"
#include "pointField.H"
#include "List.H"
#include "ggiInterpolation.H"
#include "GGIInterpolation.H"
#include "OFstream.H"
#include "transform.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Public Member Templates * * * * * * * * * * * * * //

template<class Type> 
tmp< Field<Type> > primitiveDomainScalingPatch::cloneField
(
        const Field<Type>& sourceField
) const
{
      tmp<Field<Type> > tFirst
      (
           new Field<Type>
           (
                AllFaces_.size(),
                pTraits<Type>::zero
           )
       );
      
      label numOfSubPatches = label(360)/label(diffPhi_);
      
      // Check if the geometry is correct for the field
      Field<Type>&  firstField = tFirst();

      
      // A rotation tensor  must be made to adopt the field for the rotation
      // every nPP > 0 is Rotated around the diffPhi_ 
      tensor rotationTensor = makeRotationTensor(diffPhi_);
      tensor incrementTensor = rotationTensor;
      
  
      for (label nPP = 0; nPP < numOfSubPatches; nPP++)
      {
//          // tranform must be used
          if(nPP > 0)
          {
              forAll(sourceField,i)
              {
                  firstField[i + nPP * sourceField.size()] = sourceField[i];
//                  transform
//                  (
//                      rotationTensor , sourceField[i]
//                  );
                   
              }
              // Increment the rotation for the 
              rotationTensor = incrementTensor & rotationTensor;
          }
          else
          {
            
          
             forAll(sourceField,i)
              {                  
                   firstField[i + nPP * sourceField.size()] = sourceField[i];
              }              
          }
      }
//      
//      if(rotationTensor_)
//      {
//          forAll(firstField,i)
//          {
//          firstField[i] = transform
//             (
//                 (*rotationTensor_) , firstField[i]
//             );
//          }
//      }
        
     // Info << sourceField  << endl<< firstField << endl;
      return (tFirst);
}

template<class Type> 
tmp<Field<Type> > primitiveDomainScalingPatch::cloneField
(
     const tmp<Field<Type> >& sourceField
) const
{
    tmp<Field<Type> > tFirst
    (
       new Field<Type>
       (
           AllFaces_.size(),
           pTraits<Type>::zero
       )
    );
    
    const Field<Type>& sourceFieldRef = sourceField();
    
  
    label numOfSubPatches = label(360)/label(diffPhi_);
         
    // Check if the geometry is correct for the field
         
    Field<Type>&  firstField = tFirst();

    tensor rotationTensor = makeRotationTensor(diffPhi_);
    tensor incrementTensor = rotationTensor;
        
         
    for (label nPP = 0; nPP < numOfSubPatches; nPP++)
    {
//        // tranform must be used
//        if(nPP > 0)
//        {
//           forAll(sourceField,i)
//           {
//               firstField[i + nPP * sourceFieldRef.size()]  =
//                transform
//                (
//                   rotationTensor , sourceFieldRef[i]
//                );        
//           }
//           // Increment the rotation for the 
//           rotationTensor = incrementTensor & rotationTensor;
//        }
//        else
//        {
            forAll(sourceField,i)
            {
                firstField[i + nPP * sourceFieldRef.size()] = sourceFieldRef[i];
            }              
//        }
    }
   
//    if(rotationTensor_)
//    {
//        forAll(firstField,i)
//        {
//            firstField[i] = transform
//            (
//                (*rotationTensor_) , firstField[i]
//            );
//        }
//    }
  //  Info << sourceField  << endl<< firstField << endl;
    
   return (tFirst);
}


template<class Patch>
Patch& buildSlave360Patch(const Patch& targetPatch,scalar& lastFieldAdress/*, scalar& diffPhi*/)
{ 
    //label numOfSubPatches = label(360)/label(diffPhi);
/*    
    List<face> sliceFaces;
    sliceFaces.setSize(0);
    List<face> AllFaces;
    AllFaces.setSize(0);
    
    // Patch faceNormals gives back a normalized normal.
    vectorField faceNormals = targetPatch.faceNormals();
    
    pointField origPoints = targetPatch.localPoints();
    
    List<face> localFaces = targetPatch.localFaces();

    // Average over all the face normals
    vector AverageNormal(0, 0, 0);

    forAll(faceNormals, iNormals)
    {
        AverageNormal = AverageNormal + faceNormals[iNormals];
    }

    AverageNormal = AverageNormal / faceNormals.size();

    if (mag(AverageNormal.x() ) < SMALL)
        AverageNormal.x() = 0;

    if (mag(AverageNormal.y()) < SMALL)
        AverageNormal.y() = 0;

    if (mag(AverageNormal.z() ) < SMALL)
        AverageNormal.z() = 0;
    
    for (label nPP = 0; nPP < numOfSubPatches; nPP++)
    {
        forAll(localFaces, i)
        {
            labelList llP(OrigPatch_[i].size() );
            forAll(localFaces[i], j)
            {
                llP[j] = localFaces[i][j] + nPP * ( origPoints.size() );
            }

            face tmpFace(llP);
           
            vector firstNormal = tmpFace.normal( localPoints_ );

            if(mag(firstNormal) == 0)
            {
                Info << "Error " << endl;
            }
            firstNormal /= mag(firstNormal);

            labelList llPtmp = llP;

            //TODO Make this the first step
            if (mag(firstNormal - AverageNormal) > SMALL)
            {  
                forAll(llP, k)
                {
                    llPtmp[k] = llP[llP.size()-1-k];
                }
            }
            
            llP = llPtmp;

            AllFaces.setSize(AllFaces_.size() + 1);
            AllFaces[AllFaces.size() - 1] = face(llP);
        }
        }*/
                          Patch newPatch = new Patch(/*AllFaces,origPoints*/);
    return newPatch;
}

template<class Type, class Patch> 
tmp<Field<Type> > primitiveDomainScalingPatch::interpolateToNeighbour
(
      const tmp<Field<Type> >& sourceField,  
      const Patch& targetPatch
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
   
     scalar lastFieldAdress = 0;
     
     // TODO FB Must be called otherwise the interpolation does not work correctly
      
     if(!masterToPatchPtr_)
       {
           //  the slave 360 degree patch must be built here 
           Patch slavePatch =  buildSlave360Patch(targetPatch,lastFieldAdress );
              masterToPatchPtr_ = new ggiMMInterpolation
              (
                   *primPatch_, // Master patch
                   targetPatch,  // Slave patch 
                                 // This patch must also a 360 degree patch 
                   tensorField(0),
                   tensorField(0),
                   vectorField(0)
              );
          }
//          slaveTensorList_.clear();
//          
//          const List<face>& localFaces = targetPatch.localFaces();
//          
//          const pointField& localPoints = targetPatch.localPoints();
//          forAll(localFaces,i)
//          {   
//              slaveTensorList_.setSize(slaveTensorList_.size() + 1 );
//              
//              point faceCenter = localFaces[i].centre(localPoints);
//                         
//              tensor rotation = tensorFromCCSToCartesian(faceCenter);// makeRotationTensor(theta);
//              
//              slaveTensorList_[i] = rotation;
//          }                   
           
     
     Field<Type> helpField = Field<Type>(sourceField);
     firstField = cloneField(helpField);     
        
     const vectorField& faceNormals = primPatch_->faceNormals(); 

     forAll(faceNormals,i)
     {
        faceNormals[i];        
     }
          
     tmp<Field<Type> > tResult
     (  
         new Field<Type>
         (
             targetPatch.size(),
             pTraits<Type>::zero
         )
     );
            
     Field<Type>& resultField = tResult();
     
     resultField = masterToPatchPtr_->masterToSlave(firstField);
//     forAll(resultField,j)
//     {
//         resultField[j] = 
//             transform
//             (
//                     slaveTensorList_[j] , resultField[j]
//             );
//     }     

          
    return (tResult);
}

template<class Type, class Patch> 
tmp<Field<Type> > primitiveDomainScalingPatch::interpolateToNeighbour
(
       const Field<Type>& sourceField,  
       const Patch& targetPatch
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
    
    Field<Type> helpField = sourceField;
    
    firstField = cloneField(helpField);
    //firstField = cloneField(firstField);
         
     
       
       // TODO FB Must be called otherwise the interpolation does not work correctly
    
    if(!masterToPatchPtr_)
      {
Info << "297" << endl;
             masterToPatchPtr_ = new ggiMMInterpolation
             (
                  *primPatch_,
                  targetPatch,   
                   tensorField(0),
                   tensorField(0),
                   vectorField(0)
             );
         }
         
         // Was used instead one for the radial testcase but has not shown
         // a stabilization for the solving 
         
//          slaveTensorList_.clear();
//          
//          const List<face>& localFaces = targetPatch.localFaces();
//          
//          const pointField& localPoints = targetPatch.localPoints();
//          forAll(localFaces,i)
//          {   
//              slaveTensorList_.setSize(slaveTensorList_.size() + 1 );
//              
//              point faceCenter = localFaces[i].centre(localPoints);
//                         
//              tensor rotation = tensorFromCCSToCartesian(faceCenter);// makeRotationTensor(theta);
//              
//              slaveTensorList_[i] = rotation;
//          }         
    
    
    
//    if(!masterToPatchPtr_)
//    {
//        masterToPatchPtr_ = new ggiMMInterpolation
//       (
//               *primPatch_,// Master 
//               OrigPatch_, // Slave        
//               intersection::VISIBLE,
//               intersection::VECTOR 
//       );
//    }
       // TODO interpolation pointer must be generated and holded
       
       const vectorField& faceNormals = primPatch_->faceNormals(); 

       forAll(faceNormals,i)
       {
           faceNormals[i];        
       }
         
       tmp<Field<Type> > tResult
       (  
           new Field<Type>
           (
               targetPatch.size(),
               pTraits<Type>::zero
           )
       );
         
       Field<Type>& resultField = tResult();
       //Info << "FirstField " << firstField << endl;       
       resultField = masterToPatchPtr_->masterToSlave(firstField);
//     
//       forAll(resultField,j)
//       {
//           resultField[j] = 
//               transform
//               (
//                       slaveTensorList_[j] , resultField[j]
//               );
//       }
//     
       
       return (tResult);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
