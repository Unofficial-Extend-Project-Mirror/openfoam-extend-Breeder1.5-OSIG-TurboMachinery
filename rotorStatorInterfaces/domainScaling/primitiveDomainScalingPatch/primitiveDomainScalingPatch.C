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
#include "addToRunTimeSelectionTable.H"
#include "SubField.H"
#include "entry.H"
#include "OFstream.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(primitiveDomainScalingPatch, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

primitiveDomainScalingPatch::primitiveDomainScalingPatch
(
   const List<face>& faces, 
   const pointField& points,
   const primitivePatch& OrigPatch,
   const point& CCSCenter,
   const vector& CCSaxis,
   const vector& CCSdirection,
   scalar omega,
   const word& name,
   bool switchBlockMesh
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
deltaTime_(0),
omega_(omega),
timeStep_(0),
startAngle_(0),
masterToPatchPtr_(NULL),
rotationTensor_(NULL),
//relative_(false),
master_(false),
incrementTensor_(NULL)
{    
       
   if(!switchBlockMesh&& !Pstream::parRun())
   {
       Info << "In calculateRadialGeometry" << endl;
       calcRadialGeometry();
   }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

primitiveDomainScalingPatch::~primitiveDomainScalingPatch()
{
    clearAddressing();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void primitiveDomainScalingPatch::clearAddressing()
{
    //TODO check all pointers
    deleteDemandDrivenData(primPatch_);
    deleteDemandDrivenData(masterToPatchPtr_);
}

//bool primitiveDomainScalingPatch::constraintType(const word& pt)
//{
//    return pointPatchField<scalar>::PointPatchConstructorTablePtr_->found(pt);
//}

void primitiveDomainScalingPatch::movePoints(const pointField& p)
{
    primitiveAuxPatch::movePoints(p);
}

void primitiveDomainScalingPatch::setDeltaTime(scalar timeStep, const Time& runTime)
{    
   if( (timeStep_ < timeStep) )
   {       
          deltaTime_ = timeStep - timeStep_;
          timeStep_ = timeStep;
          deleteDemandDrivenData(masterToPatchPtr_);
          calcuteMovedPoints();
          scalar theta = omega_ *360* timeStep_/scalar(60.0);
          incrementTensor_ = new tensor( makeRotationTensor(theta) );
   } 
   
//   scalar numberOfOuterTimeStepsPerRevolution = GREAT;
//
//   if (runTime.controlDict().found("numberOfOuterTimeStepsPerRevolution"))
//   (
//        numberOfOuterTimeStepsPerRevolution = readScalar(runTime.controlDict().lookup("numberOfOuterTimeStepsPerRevolution"))
//   );
//   
//   scalar numRevolutionsPerMin = GREAT;
//
//   if (runTime.controlDict().found("numRevolutionsPerMin"))
//   (
//       numRevolutionsPerMin = readScalar(runTime.controlDict().lookup("numRevolutionsPerMin"))
//   );
// 
//   scalar outerTimeStep = 60/(numRevolutionsPerMin*numberOfOuterTimeStepsPerRevolution);
//   
//   if(omega_< 0)
//       omega_ = -numRevolutionsPerMin *360/60 * outerTimeStep;
//   else
//      omega_ = numRevolutionsPerMin *360/60 * outerTimeStep;
// 
//   dimensionedScalar timeStepNew = runTime.deltaT();
//   
//   scalar diff = (timeStep + timeStepNew.value()) - (timeStep_ + outerTimeStep);
//   if(diff < 0)
//       diff= diff * scalar(-1);
//   
//   if( diff < 0.0000001)
////   if( (timeStep_ < timeStep) )//&& ( omega_ == 0 ) )
//   {       
//      timeStep_ = timeStep + timeStepNew.value();
//      deltaTime_ = outerTimeStep;
//      deleteDemandDrivenData(neighbourToOrig_);
//      deleteDemandDrivenData(masterToPatchPtr_);
//      calcuteMovedPoints();
//      Info << "NextSwitchTimeStep " << timeStep_<< endl; 
//   }
   
}

void primitiveDomainScalingPatch::adoptToOuterTimeStep(scalar angle)
{
      deleteDemandDrivenData(masterToPatchPtr_);
      
      localPoints_ = cylindricalCS_.localPosition( localPoints_ );
      
  //    localPoints_.y() = localPoints_.y() + omega_*deltaT_;

      
      forAll(localPoints_, i)
      {
          localPoints_[i].y() = localPoints_[i].y() + angle;
      }
      
      localPoints_ = cylindricalCS_.globalPosition( localPoints_ );
      
      deleteDemandDrivenData( primPatch_ );
      
      primPatch_ = new primitiveAuxPatch(AllFaces_, localPoints_ );
}

const List<face>& primitiveDomainScalingPatch::GetFaces() const
{
    return ( AllFaces_ );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "calculateDomainScalingGeometry.C"
#   include "primitiveDomainScalingPatchTemplateSpecialisations.C"


// ************************************************************************* //


