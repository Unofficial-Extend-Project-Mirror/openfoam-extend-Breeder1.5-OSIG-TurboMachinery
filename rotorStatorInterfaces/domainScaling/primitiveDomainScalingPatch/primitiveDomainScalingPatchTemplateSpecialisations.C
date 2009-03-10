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


tmp< Field<vector> > primitiveDomainScalingPatch::cloneField
(
        const Field<vector>& sourceField
) const
{
      tmp<Field<vector> > tFirst
      (
           new Field<vector>
           (
                AllFaces_.size(),
                pTraits<vector>::zero
           )
       );
      
      label numOfSubPatches = label(360)/label(diffPhi_);
      
      // Check if the geometry is correct for the field
      Field<vector>&  firstField = tFirst();
  
      for (label nPP = 0; nPP < numOfSubPatches; nPP++)
      {
//          // tranform must be used
          forAll(sourceField,i)
          {
             firstField[i + nPP * sourceField.size()] = sourceField[i];
             // Was used for the radial testcase but has not shown
             // a stabilization for the solving 
//             transform
//             (
//               tensorList_[i + nPP * sourceField.size()].T()  /*rotationTensor*/ , sourceField[i]
//             );
          }
      }
        
      return (tFirst);
}


tmp<Field<vector> > primitiveDomainScalingPatch::cloneField
(
     const tmp<Field<vector> >& sourceField
) const
{
    tmp<Field<vector> > tFirst
    (
       new Field<vector>
       (
           AllFaces_.size(),
           pTraits<vector>::zero
       )
    );
    
    const Field<vector>& sourceFieldRef = sourceField();
    
  
    label numOfSubPatches = label(360)/label(diffPhi_);
         
    // Check if the geometry is correct for the field
         
    Field<vector>&  firstField = tFirst();

            
         
    for (label nPP = 0; nPP < numOfSubPatches; nPP++)
    {
//        // tranform must be used
           forAll(sourceFieldRef,i)
           {
               firstField[i + nPP * sourceFieldRef.size()]  =  sourceFieldRef[i];
               // Was used for the radial testcase but has not shown
               // a stabilization for the solving 
//                transform
//                (
//                   tensorList_[i + nPP * sourceFieldRef.size()].T() /*rotationTensor*/ , sourceFieldRef[i]
//                );        
           }
    }
   
    
   return (tFirst);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
