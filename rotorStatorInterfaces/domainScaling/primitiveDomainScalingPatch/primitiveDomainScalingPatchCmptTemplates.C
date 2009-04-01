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
#include "diagTensorField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Public Member Templates * * * * * * * * * * * * * //

template<class Type> 
tmp< Field<Type> > primitiveDomainScalingPatch::cloneField
(
        const Field<Type>& sourceField,
        const direction& cmpt,
        const int& rank
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
         
      for (label nPP = 0; nPP < numOfSubPatches; nPP++)
      {
//          // tranform must be used

            //
            forAll(sourceField,i)              
            {
              scalar scale = scalar(1);
              // Was used instead one for the radial testcase but has not shown
              // a stabilization for the solving 
//                    pow(diag(tensorList_[i + nPP * sourceField.size()].T() /*rotationTensor*/).component(cmpt), rank);                  
              firstField[i + nPP * sourceField.size()] = scale * sourceField[i]; 
            }
      }

        
     // Info << sourceField  << endl<< firstField << endl;
      return (tFirst);
}

template<class Type> 
tmp<Field<Type> > primitiveDomainScalingPatch::cloneField
(
     const tmp<Field<Type> >& sourceField,
     const direction& cmpt,
     const int& rank
     
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
     
    for (label nPP = 0; nPP < numOfSubPatches; nPP++)
    {
           forAll(sourceField,i)
           {
               scalar scale = scalar(1);
               // Was used instead one for the radial testcase but has not shown
               // a stabilization for the solving 
                      //  pow(diag(tensorList_[i + nPP * sourceFieldRef.size()].T() /*rotationTensor*/).component(cmpt), rank);                
               firstField[i + nPP * sourceFieldRef.size()]  = scale * sourceFieldRef[i];
           }
    }
   
  //  Info << sourceField  << endl<< firstField << endl;
    
   return (tFirst);
}

template<class Type, class Patch> 
tmp<Field<Type> > primitiveDomainScalingPatch::interpolateToNeighbour
(
      const tmp<Field<Type> >& sourceField,  
      const Patch& targetPatch,
      const direction& cmpt,
      const int& rank
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
      
     if(!masterToPatchPtr_)
       { 
Info << "162" << endl;
	  masterToPatchPtr_ = new ggiMMInterpolation
          (
              *primPatch_, // Master patch
              targetPatch,  // Slave patch            
              tensorField(targetPatch.size(),
	      tensor(1,0,0,0,1,0,0,0,1) ),
	      tensorField( (*primPatch_).size(),
	      tensor(1,0,0,0,1,0,0,0,1) )
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
       
     
     
     Field<Type> helpField = Field<Type>(sourceField);
     
     firstField = cloneField(helpField,cmpt,rank);     
        
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
//         scalar scale =
//                   pow(diag(slaveTensorList_[j]).component(cmpt), rank);           
//         resultField[j] = scale* resultField[j];
//     }

          
    return (tResult);
}

template<class Type, class Patch> 
tmp<Field<Type> > primitiveDomainScalingPatch::interpolateToNeighbour
(
       const Field<Type>& sourceField,  
       const Patch& targetPatch,
       const direction& cmpt,
       const int& rank
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
    
    firstField = cloneField(helpField,cmpt,rank);     
           
    //firstField = cloneField(firstField);
         
    
       
       // TODO FB Must be called otherwise the interpolation does not work correctly
    
    if(!masterToPatchPtr_)
      {
         Info << "257" << endl;
      
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
//         slaveTensorList_.clear();
//         
//         const List<face>& localFaces = targetPatch.localFaces();
//         const pointField& localPoints = targetPatch.localPoints();
//         forAll(localFaces,i)
//         {  
//              slaveTensorList_.setSize(slaveTensorList_.size() + 1 );
//              
//              point faceCenter = localFaces[i].centre(localPoints);
//                         
//              tensor rotation = tensorFromCCSToCartesian(faceCenter);// makeRotationTensor(theta);
//              
//              slaveTensorList_[i] = rotation;
//          }         
      
    
   
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
       
//       forAll(resultField,j)
//       {
//           scalar scale =
//                     pow(diag(slaveTensorList_[j]).component(cmpt), rank);           
//           resultField[j] = scale* resultField[j];
//       }
       
       return (tResult);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
