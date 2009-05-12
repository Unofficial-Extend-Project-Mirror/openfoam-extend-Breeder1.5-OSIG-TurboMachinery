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

 Class
 PrimitiveRadialPatch class 

 Description
 Calculates the new face geometry by respecting the old clipping of the 
 original patch

 SourceFiles    

 \*---------------------------------------------------------------------------*/
#include "primitiveDomainScalingPatch.H"
#include "face.H"
#include "mathematicalConstants.H"
#include "coordinateSystem.H"
#include "triPointRef.H"

namespace Foam

{

void primitiveDomainScalingPatch::calcRadialGeometry()
{
    labelList MinRB;
    labelList MaxRB;
    labelList MinCB;
    labelList MaxCB;

    Info << "In calculateRadialGeometry" << endl; 
    
    DefineLocalExternalEdges(OrigExtEdges(OrigPatch_), OrigPatch_.localPoints());

    FindCorners();

    Info << "After Find Corners" << endl;
    generateCCSBorderPoints();

    Info << "Before GetBorders" << endl;
    getBorders(MinRB, MaxRB, MinCB, MaxCB); // Should deliver more than here


    labelListList SF;

    Info << "generate PointList " << endl;
    generatePointList(MaxCB, MinCB);

    generateNewFaceList();
    
    const pointField&  test = OrigPatch_.localPoints();    
    const List<face>&   faces = OrigPatch_.localFaces();   
          
    word fileNameOne = "FacesOne" + name_ + ".obj";;
      OFstream os(fileNameOne);
      forAll( faces, i)
      {
         forAll(faces[i], j)
         {
                writeOBJ(
                        os,
                        test[faces[i][j] ] );
        } 
      }
      forAll(faces,i)
      {
                
                os << "f " ; 
                forAll(faces[i],j)
                {
                    os << faces[i][j] + 1 << " " ; 
                                                  
                }
                os << endl;
      }  
    
     word fileName = "faces" + name_;
     OFstream osNew(fileName);
          
     
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
     
     if(startAngle_)
     {   
         localPoints_ = cylindricalCS_.localPosition( localPoints_ );
             
         forAll(localPoints_, i)
         {
             localPoints_[i].y() = localPoints_[i].y() + startAngle_;
         }
         
         localPoints_ = cylindricalCS_.globalPosition( localPoints_ );
     }
     
     if(!primPatch_)
         primPatch_ = new primitiveAuxPatch(AllFaces_, localPoints_ );
}

scalar primitiveDomainScalingPatch::calculateRotationalDivisor(const labelList& MaxCB, const labelList& MinCB)
{
    scalar minPhi = CCSBorderPoints_[BorderEdges_[ MinCB[0] ][0] ].y();
    
    // TODO should be removed 
    forAll(MinCB, i)
    {
        scalar reference = CCSBorderPoints_[BorderEdges_[ MinCB[i] ][0] ].y();
        if( (reference - minPhi) > 10)
        {

            Info  << "primitiveDomainScalingPatch::CalculateRotationalDivisor(const labelList& MinCB, const labelList& MaxCB) const" << endl
                  << "The MinCB line has ambigous phi values. Please check the geometry"<< endl  
                  << " the borders of your domainScalingPatch must lie on a" << endl
                  << " straight line" << endl
                  << endl;            
        }
    }
    
    scalar maxPhi = CCSBorderPoints_[BorderEdges_[ MaxCB[0] ][0] ].y();
    // TODO should be removed    
    forAll(MaxCB, i)
    {
       scalar reference = CCSBorderPoints_[BorderEdges_[ MinCB[i] ][0] ].y();
       if( (reference - minPhi) > 0.01)
       {
            Info  << "primitiveDomainScalingPatch::CalculateRotationalDivisor(const labelList& MinCB, const labelList& MaxCB) const" << endl
                  << "The MaxCB line has ambigous phi values. Please check the geometry"  << endl
                  << " the borders of your domainScalingPatch must lie on a" << endl
                  << " straight line" << endl
                  << endl;
       }
    }
    
    Info << "MaxPhi " << maxPhi << endl;
    Info << "MinPhi " << minPhi << endl;
    
    scalar diffPhi = (maxPhi-minPhi);
    Info << "diffPhi " << diffPhi << endl;
    if(diffPhi == 0)
    {
        Info        << "primitiveDomainScalingPatch::CalculateRotationalDivisor(const labelList& MinCB, const labelList& MaxCB) const" << endl
                    << "The Max phi border and the minium phi border have the same Phi values. Please check the geometry" << endl
                    << endl;
    }
    
    
   // scalar remainder = fmod(360, diffPhi);
    
    Info << "diffPhi " << diffPhi << endl;
  //  Info << "Remainder " << remainder << endl;
//    if( remainder > 0.01  )        
//    {
//       FatalErrorIn("primitiveDomainScalingPatch::CalculateRotationalDivisor(const labelList& MinCB, const labelList& MaxCB) const")
//                 << "The Max phi border and the minium phi border have the same Phi values. Please check the geometry"
//                 << abort(FatalError);
//    }
    if(diffPhi < 0)
        diffPhi = -1* diffPhi;
        
    Info << "diffPhi " << diffPhi << endl;
    

    
    return diffPhi; 
}

void primitiveDomainScalingPatch::generatePointList(const labelList& MaxCB,
        const labelList& MinCB )
{
    diffPhi_ =  calculateRotationalDivisor(MinCB, MaxCB);
   
    label numOfSubPatches = label(360)/label(diffPhi_);
    
    pointField origPoints = cylindricalCS_.localPosition( OrigPatch_.localPoints() );
    
    label origSize = origPoints.size();
   
    localPoints_.setSize(0);
    for (label i = 0; i < label(numOfSubPatches); i++) 
    {        
        forAll(origPoints,j)
        {
            // The cloning and the rotation of the points takes place here
            localPoints_.setSize(localPoints_.size() + 1 );
            localPoints_[origSize*i+j] = origPoints[j];
            localPoints_[origSize*i+j].y() = localPoints_[origSize*i+j].y() + 
            diffPhi_ * i;
        }
        Info << i << endl;
    }
    localPoints_ = cylindricalCS_.globalPosition(localPoints_);
   
    extractNewOrigPatch();
}

void primitiveDomainScalingPatch::extractNewOrigPatch()
{
    pointField globalPForig = OrigPatch_.points();
    
    // Generation of the newOrigPatch
    
    // TODO check if necessary 
    newOrigPatch_.setSize(0);
    newOrigPoints_.setSize(0);
    
    forAll(OrigPatch_,i)
    {
        labelList tmpFlabel(0);
        forAll(OrigPatch_[i],j)
        {
            tmpFlabel.setSize(tmpFlabel.size() + 1);
            tmpFlabel[tmpFlabel.size() - 1] = newOrigPoints_.size();
            
            newOrigPoints_.setSize(newOrigPoints_.size() + 1 );                  
            newOrigPoints_[newOrigPoints_.size() - 1] = 
                globalPForig[ OrigPatch_[i][j] ];
        }
        face tmpFace(tmpFlabel);
        newOrigPatch_.setSize( newOrigPatch_.size() + 1 );
        newOrigPatch_[newOrigPatch_.size() - 1 ] = tmpFace;
    }
}

void primitiveDomainScalingPatch::generateNewFaceList()
{    
    label numOfSubPatches = label(360)/label(diffPhi_);
    
    sliceFaces_.setSize(0);
    AllFaces_.setSize(0);
    
    // Patch faceNormals gives back a normalized normal.
    vectorField faceNormals = OrigPatch_.faceNormals();
    
    pointField origPoints = OrigPatch_.localPoints();
    
    List<face> localFaces = OrigPatch_.localFaces();

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
            // At that time the localPoints_ contain the rotated and cloned points
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

            AllFaces_.setSize(AllFaces_.size() + 1);
            AllFaces_[AllFaces_.size() - 1] = face(llP);
            
        }
        
    }  
}

void primitiveDomainScalingPatch::calcuteMovedPoints()
{
    //This is done with a tensor multiplication because it is easier to track the 
    //rotation with a tensor multiplication
    
    if(omega_ != scalar(0) )
    {
        Info << "In rotation"<< endl;
        scalar theta = omega_ *360* deltaTime_/scalar(60.0);
        
        deleteDemandDrivenData( rotationTensor_ );
        rotationTensor_ = new tensor(makeRotationTensor(theta));
   
        localPoints_ = (*rotationTensor_) & localPoints_;
        deleteDemandDrivenData( primPatch_ );
            
        primPatch_ = new primitiveAuxPatch(AllFaces_, localPoints_ );
    }
    // deletes the interpolation pointer which leads to a rebuild of the 
    // face weights for the next interpolation call
    deleteDemandDrivenData(masterToPatchPtr_);
   
    
    const pointField&  test = OrigPatch_.localPoints();    
    const List<face>&   faces = OrigPatch_.localFaces();   
    word fileNameOne = "Faces" + name_ + ".obj";
        OFstream os(fileNameOne);      
      
     
      forAll( faces, i)
      {
         forAll(faces[i], j)
         {
                writeOBJ(
                        os,
                        test[faces[i][j] ] );
        } 
      }
      forAll(faces,i)
      {
                
                os << "f " ; 
                forAll(faces[i],j)
                {
                    os << faces[i][j] + 1 << " " ; 
                                                  
                }
                os << endl;
      } 
    
      
     word fileName = "faces" + name_;
     OFstream osNew(fileName);
            
       
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
