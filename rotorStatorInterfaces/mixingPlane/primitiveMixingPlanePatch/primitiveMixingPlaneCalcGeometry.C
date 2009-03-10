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
#include "primitiveMixingPlanePatch.H"
#include "face.H"
#include "mathematicalConstants.H"
#include "coordinateSystem.H"
#include "triPointRef.H"

namespace Foam

{

void primitiveMixingPlanePatch::calcRadialGeometry()
{
    bool bdebug = false;
    labelList MinRB;
    labelList MaxRB;
    labelList MinCB;
    labelList MaxCB;

    DefineLocalExternalEdges(OrigExtEdges(OrigPatch_), OrigPatch_.localPoints());

    FindCorners();

    generateCCSBorderPoints();
    Info << "MinCB before gB" << MinCB << endl;
    getBorders(MinRB, MaxRB, MinCB, MaxCB); // Should deliver more than here
     
    Info << "MinCB before gP" << MinCB << endl;
    
    labelListList SF;

    generatePointList(MaxRB, MinCB);

    generateNewFaceList();

    OFstream os("Faces.obj");
   
    if(bdebug)
    {    
        
       // Info << "TEST " << AllFaces_[0].centre(localPoints_) << endl;
        writeOBJ(os, cylindricalCS_.globalPosition( localPoints_[AllFaces_[1][0] ]) );
        writeOBJ(os, cylindricalCS_.globalPosition( localPoints_[AllFaces_[1][1] ]) );
        writeOBJ(os, cylindricalCS_.globalPosition( localPoints_[AllFaces_[1][2] ]) );
        writeOBJ(os, AllFaces_[1].centre(cylindricalCS_.globalPosition( localPoints_) ) );
     
    
        forAll(AllFaces_, i)
        {
            forAll(AllFaces_[i], j)
            {


                writeOBJ(
                        os,
                        localPoints_[AllFaces_[i][j] ] );
            }
            writeOBJ(
                    os,
                    AllFaces_[i].centre(localPoints_ ) );
        }
    }
    
    if(true)
    {
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
    
    generateCCSTransformationTensors();
}


void primitiveMixingPlanePatch::rebuildPatchForEdges
(
        const List<point>& edgePointList
)
{
    bool bdebug = false;
    labelList MinRB;
    labelList MaxRB;
    labelList MinCB;
    labelList MaxCB;
    
    getBorders(MinRB, MaxRB, MinCB, MaxCB);
    
    generateForeignMatchedPointList(edgePointList, MinCB);
        
    generateNewFaceList();
    
    OFstream os("FacesRebuild.obj");
    if(bdebug)
      {    
          
         //Info << "TEST " << AllFaces_[0].centre(localPoints_) << endl;
          writeOBJ(os, cylindricalCS_.globalPosition( localPoints_[AllFaces_[1][0] ]) );
          writeOBJ(os, cylindricalCS_.globalPosition( localPoints_[AllFaces_[1][1] ]) );
          writeOBJ(os, cylindricalCS_.globalPosition( localPoints_[AllFaces_[1][2] ]) );
          writeOBJ(os, AllFaces_[1].centre(cylindricalCS_.globalPosition( localPoints_) ) );
       }
      
      forAll(AllFaces_, i)
      {
          forAll(AllFaces_[i], j)
          {


              writeOBJ(
                      os,
                      localPoints_[AllFaces_[i][j] ] );
          }
          writeOBJ(
                  os,
                  AllFaces_[i].centre(localPoints_ ) );
      }
      
      if(true)
      {
           OFstream osNew("facesRebuild.obj");
           
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
}


List<point> primitiveMixingPlanePatch::calculateForeignMatchedPoints
(
        const List<point>& foreignPoints ,                 
        MachineType mtype
)
{
    List<point> returnPoints = foreignPoints;
    
    forAll(foreignPoints,i)
    {
       returnPoints[i].y() = pointsSide_[0].y();
    }

    return (returnPoints );
}

void primitiveMixingPlanePatch::generatePointList
(
   const labelList& MaxRB, 
   const labelList& MinCB
)
{    
    List<List<point> > orderedPoints;
    bool bdebug = false;
    bool bCSV = false;

    // Finding the starting point for the upper points
    // That is done by finding the matching points between the MaxRB and 
    // MinCB because the      
    
    makeSidePoints(MinCB);
    
    if(bdebug)
    {
        OFstream opS("pointsSide.txt");
        opS << pointsSide_ << endl;
    }
  
    makeMaxRBPoints(MaxRB);

    //      After the pointsSide are generated 
    //      They can be easily transformed by tensor operations
    //      The necessary delta angles are extracted by from the MaxRB 

    OFstream osS("orderedPoints.csv");
    
    OFstream osTest("maxRBallPoints.obj");

    forAll(maxRbPoints_, j)
    {
        writeOBJ(osTest, cylindricalCS_.globalPosition(maxRbPoints_[j]));

    }

    forAll(pointsSide_, j)
    {
        writeOBJ(osTest, cylindricalCS_.globalPosition(pointsSide_[j]));

    }
    // Until here it looks well

    // Now the angles are extracted to a list
    //      osS << "RB Points size " << maxRbPoints_.size() << endl;
    List<scalar> angleList(0);
//    OFstream osAngles("angles.txt");
   
    forAll(maxRbPoints_, i)
    {
        angleList.setSize(angleList.size() +1);
    
//        osAngles << maxRbPoints_[i].y() << endl;
        
        angleList[i] = maxRbPoints_[i].y();
    }
    angleList = mergesort( angleList );
//    osAngles << "Angles Sorted " << endl;
    forAll(angleList,i)
    {
//        osAngles << angleList[i] << endl;
    }
    
    
//    osAngles << "MinCB Y" << endl;
//    forAll(pointsSide_, i)
//    {
//        osAngles << pointsSide_[i].y() << endl;       
//    }



    // Todo are the pointsSide ordered from upside down or from down to up
    OFstream oLast("LastTest.txt");
    // pointsSide is in cylindrical coordinates
    forAll(pointsSide_, i)
    {
        //         osS << "in For pointsSide" << endl;
        orderedPoints.setSize(orderedPoints.size() + 1);

        orderedPoints[i].setSize(orderedPoints[i].size() + 1);
        orderedPoints[i][orderedPoints[i].size() -1 ] = pointsSide_[i];

        oLast << "Side point" << pointsSide_[i] << endl;
        forAll(angleList, j)
        {
            if (j > 0)
            {
                point temp(pointsSide_[i].x(), angleList[j], pointsSide_[i].z());
                oLast <<  temp << endl;
                
                orderedPoints[i].setSize(orderedPoints[i].size() + 1);
                orderedPoints[i][orderedPoints[i].size() -1 ] = temp;

//                osAngles << "J " << j << " Angle " << angleList[j] << endl;
                //        osS << "in if" << endl;
            }
            //      osS << "in For angleList" << endl;
        }
    }
    // osS << angleList.size() << endl;
    OFstream osP("orderedPoints.obj");
    
    
    oLast << CCSBorderPoints_ << endl;

    forAll(orderedPoints, i)
    {
        //          osS << orderedPoints[i].size() << endl;
    }
    
    if (true)
    {
        forAll(orderedPoints, i)
        {   
            forAll(orderedPoints[i], j)
            {
                writeOBJ(osP,
                        cylindricalCS_.globalPosition(orderedPoints[i][j]));
                    osS << orderedPoints[i][j].x() << ", "
                        << orderedPoints[i][j].y() << ", "
                        << orderedPoints[i][j].z() << endl;
            }
        }
    }
    
    writeOBJ(osP,
            cylindricalCS_.origin());
    
    OFstream osCentre("centre.obj");
    writeOBJ(osCentre,
               cylindricalCS_.origin());

    // Generating the pointfield and the labels
    label Counter = -1;
    localPoints_.setSize(0);
    forAll(orderedPoints, i)
    {
        pointLabels_.setSize(i + 1);
        forAll(orderedPoints[i], j)
        {
            Counter++;
            pointLabels_[i].setSize(j + 1);
            pointLabels_[i][j] = Counter;

            localPoints_.setSize(localPoints_.size() + 1);
            localPoints_[localPoints_.size() - 1] = cylindricalCS_.globalPosition(orderedPoints[i][j] );
        }
    }
    
    OFstream of("localPoints.obj");
    forAll(localPoints_,i)
    {
        writeOBJ(of,localPoints_[i]);
        
    }

    OFstream os("test.obj");
    forAll(orderedPoints, i)
    {
        forAll(orderedPoints[i], j)
        {
            writeOBJ(os, cylindricalCS_.globalPosition(orderedPoints[i][j]) );
            if (bCSV)
                Info << cylindricalCS_.globalVector( orderedPoints[i][j] ).x() << ", "
                        << cylindricalCS_.globalVector( orderedPoints[i][j] ).y() << ", "
                        << cylindricalCS_.globalVector( orderedPoints[i][j] ).z() << endl;
        }
    }

    writeOBJ(os, cylindricalCS_.origin() );

    if (bdebug)
        Info << " PointsInternal " << orderedPoints << endl;

    if (bdebug)
        Info << " BorderPoints CS " << CCSBorderPoints_ << endl;
    

}

void primitiveMixingPlanePatch::generateForeignMatchedPointList
(
   const List<point>& matchedSidePoints, 
   const labelList& MinCB
)
{
    List<List<point> > orderedPoints;
    bool bdebug = false;
    bool bCSV = false;

    // Finding the starting point for the upper points
    // That is done by finding the matching points between the MaxRB and 
    // MinCB because the      
    
    OFstream opS("pointsSide.txt");
    
    pointsSide_ = calculateForeignMatchedPoints(matchedSidePoints, mType_);
    
    if(bdebug)
        opS << pointsSide_ << endl;


    OFstream osTest("maxRBallPoints.obj");

    forAll(maxRbPoints_, j)
    {
        writeOBJ(osTest, cylindricalCS_.globalPosition(maxRbPoints_[j]));

    }

    forAll(pointsSide_, j)
    {
        writeOBJ(osTest, cylindricalCS_.globalPosition(pointsSide_[j]));

    }
    // Until here it looks well
    
    // Now the angles are extracted to a list
    //      osS << "RB Points size " << maxRbPoints_.size() << endl;
    List<scalar> angleList(0);
    OFstream osAngles("angles.txt");
    forAll(maxRbPoints_, i)
    {
        angleList.setSize(angleList.size() +1);
        osAngles << maxRbPoints_[i].y() << endl;
        angleList[i] = maxRbPoints_[i].y();
    }
    angleList = mergesort( angleList );
    osAngles << "Angles Sorted " << endl;
    forAll(angleList,i)
    {
        osAngles << angleList[i] << endl;
    }
    
    
    osAngles << "MinCB Y" << endl;
    forAll(pointsSide_, i)
    {
        osAngles << pointsSide_[i].y() << endl;       
    }

    // Sorting the angleList from min angle to max angle
    

    //      osS << angleList << endl; 

    OFstream oLast("LastTest.txt");
    forAll(pointsSide_, i)
    {
        //         osS << "in For pointsSide" << endl;
        orderedPoints.setSize(orderedPoints.size() + 1);

        orderedPoints[i].setSize(orderedPoints[i].size() + 1);
        orderedPoints[i][orderedPoints[i].size() -1 ] = pointsSide_[i];

        oLast << "Side point" << pointsSide_[i] << endl;
        forAll(angleList, j)
        {
            if (j > 0)
            {
                point temp(pointsSide_[i].x(), angleList[j], pointsSide_[i].z());
                oLast <<  temp << endl;
                
                orderedPoints[i].setSize(orderedPoints[i].size() + 1);
                orderedPoints[i][orderedPoints[i].size() -1 ] = temp;

                osAngles << "J " << j << " Angle " << angleList[j] << endl;
                //        osS << "in if" << endl;
            }
            //      osS << "in For angleList" << endl;
        }
    }
    
    // osS << angleList.size() << endl;
    
    OFstream osP("orderedPoints.obj");
    
    
    oLast << CCSBorderPoints_ << endl;

   
    forAll(orderedPoints, i)
    {
        //          osS << orderedPoints[i].size() << endl;
    }
    
    writeOBJ(osP,
            cylindricalCS_.origin());
    
    OFstream osCentre("centre.obj");
    writeOBJ(osCentre,
               cylindricalCS_.origin());

    // Generating the pointfield and the labels
    label Counter = -1;
    
    localPoints_.setSize(0);
    forAll(orderedPoints, i)
    {
        pointLabels_.setSize(i + 1);
        forAll(orderedPoints[i], j)
        {
            Counter++;
            pointLabels_[i].setSize(j + 1);
            pointLabels_[i][j] = Counter;

            localPoints_.setSize(localPoints_.size() + 1);
            localPoints_[localPoints_.size() - 1] = cylindricalCS_.globalPosition(orderedPoints[i][j] );
        }
    }
    
    OFstream of("localPoints.obj");
    forAll(localPoints_,i)
    {
        writeOBJ(of,localPoints_[i]);
        
    }

    OFstream os("test.obj");
    forAll(orderedPoints, i)
    {
        forAll(orderedPoints[i], j)
        {
            writeOBJ(os, cylindricalCS_.globalPosition(orderedPoints[i][j]) );
            if (bCSV)
                Info << cylindricalCS_.globalVector( orderedPoints[i][j] ).x() << ", "
                        << cylindricalCS_.globalVector( orderedPoints[i][j] ).y() << ", "
                        << cylindricalCS_.globalVector( orderedPoints[i][j] ).z() << endl;
        }
    }

    writeOBJ(os, cylindricalCS_.origin() );

    if (bdebug)
        Info << " PointsInternal " << orderedPoints << endl;

    if (bdebug)
        Info << " BorderPoints CS " << CCSBorderPoints_ << endl;
    
   

    // the first point must be the 

    //The index of the MinCB is also used for buiding the points from upside down to the centre

}

void primitiveMixingPlanePatch::generateNewFaceList()
{
    bool bdebug = false;
    label counterfaces = -1;

    sliceFaces_.setSize(0);
    AllFaces_.setSize(0);

    // Patch faceNormals gives back a normalized normal.
    vectorField faceNormals = OrigPatch_.faceNormals();

    // Average over all the face normals
    vector AverageNormal(0, 0, 0);

    forAll(faceNormals, iNormals)
    {
        AverageNormal = AverageNormal + faceNormals[iNormals];
    }

    AverageNormal = AverageNormal / faceNormals.size();

    OFstream os("Mags.txt");

    if (mag( AverageNormal.x() ) < SMALL)
        AverageNormal.x() = 0;

    if ( mag(AverageNormal.y()) < SMALL)
        AverageNormal.y() = 0;

    if ( mag(AverageNormal.z() ) < SMALL)
        AverageNormal.z() = 0;

    
    forAll(pointLabels_, i)
    {

        labelList llP(3);

        if (i == (pointLabels_.size()-1))
            return;

        sliceFaces_.setSize(sliceFaces_.size() + 1);

        if (bdebug)
            Info << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!In if "
                    << pointLabels_[i].size()-1 << endl;

        for (label j = 0; j < pointLabels_[i].size()-1; j++)
        {

            vector FirstNormal = triPointRef
            (
                     localPoints_[pointLabels_[i][j+1]],
                     localPoints_[pointLabels_[i][j]] ,
                    localPoints_[pointLabels_[i+1][j]] 
            ).normal();

            FirstNormal /= mag(FirstNormal);

            if (mag(FirstNormal + AverageNormal) < SMALL)
            {             
                llP[0] = pointLabels_[i][j+1];
                llP[1] = pointLabels_[i][j];
                llP[2] = pointLabels_[i+1][j];
            } else
            {
                llP[1] = pointLabels_[i][j+1];
                llP[0] = pointLabels_[i][j];
                llP[2] = pointLabels_[i+1][j];

            }
            AllFaces_.setSize(AllFaces_.size() + 1);
            AllFaces_[AllFaces_.size() - 1] = face(llP);

            sliceFaces_[i].setSize(sliceFaces_[i].size() + 1);
            sliceFaces_[i][sliceFaces_[i].size() - 1] = AllFaces_.size()-1;

            FirstNormal = triPointRef
            (
                    localPoints_[pointLabels_[i][j+1]],
                    localPoints_[pointLabels_[i+1][j]],
                    localPoints_[pointLabels_[i+1][j+1]] 
            ).normal();


            FirstNormal /= mag(FirstNormal);

            if (mag(FirstNormal + AverageNormal) < SMALL)
            {
                llP[0] = pointLabels_[i][j+1];
                llP[1] = pointLabels_[i+1][j];
                llP[2] = pointLabels_[i+1][j+1];
            } else
            {
                llP[1] = pointLabels_[i][j+1];
                llP[0] = pointLabels_[i+1][j];
                llP[2] = pointLabels_[i+1][j+1];
            }

            
                os << "Mag " << mag(FirstNormal - AverageNormal)
                        << "FirsNormal " << FirstNormal << "AverageNormal "
                        << AverageNormal << endl;

            AllFaces_.setSize(AllFaces_.size() + 1);
            AllFaces_[AllFaces_.size() - 1] = face(llP);

            counterfaces++;
            sliceFaces_[i].setSize(sliceFaces_[i].size() + 1);
            sliceFaces_[i][sliceFaces_[i].size() - 1] = AllFaces_.size()-1;
        }
    }    
}


void primitiveMixingPlanePatch::makeSidePoints( const labelList& MinCB )
{   
    forAll(MinCB, i)
    {
        pointsSide_.setSize(pointsSide_.size() + 1);
//         Info << "MinCB " << MinCB[i] << " Size" << MinCB.size() << endl
//              << "BorderEdges_ " << BorderEdges_[ MinCB[i] ][0] << "Size " 
// 	     << BorderEdges_.size() << "CCSBorderPoints_ " << endl;
//         Info << "Size " << CCSBorderPoints_.size() << " CCSBorderPoints " << CCSBorderPoints_[BorderEdges_[ MinCB[i] ][0] ] <<  endl;
// 	Info << "Size " << pointsSide_.size() << endl;
// 	Info << CCSBorderPoints_[label(BorderEdges_[ label(MinCB[i]) ][0] )] << endl;
// 	Info << "PointsSide "<< pointsSide_.size() - 1 << endl;
//         Info << "CCS Border points " << CCSBorderPoints_ << endl;
//         Info << "MinCB " << MinCB << endl;
// 	Info << "MinCB label " << label(MinCB[i]) << endl;
//         Info << "BorderEdges label " << BorderEdges_[ label(MinCB[i]) ][0]  << endl;
// 	Info << "CCS " << CCSBorderPoints_[label(BorderEdges_[ label(MinCB[i]) ][0] )] << endl;
        pointsSide_[pointsSide_.size() - 1 ] =  CCSBorderPoints_[label(BorderEdges_[ label(MinCB[i]) ][0] )];
        pointsSide_.setSize(pointsSide_.size() + 1);
        pointsSide_[pointsSide_.size() - 1 ] = CCSBorderPoints_[BorderEdges_[ MinCB[i] ][1] ];
    }
    
  
    
    List<point> temppoint(0);
    forAll(pointsSide_, i)
    {
        bool bdouble = false;
        forAll(pointsSide_, j)
        {
            if ( (pointsSide_[i] == pointsSide_[j]) && (i != j))
                bdouble = true;
        }
        if (!bdouble)
        {
            temppoint.setSize(temppoint.size() + 1);
            temppoint[temppoint.size() - 1] = pointsSide_[i];
        }
    }

    temppoint.setSize(0);
    forAll(pointsSide_, i)
    {
        bool bdouble = false;
        forAll(temppoint, j)
        {
            if(mType_ != Axial)
            {
                if ( (temppoint[j].z() == pointsSide_[i].z()))
                      bdouble = true;
            }
            else
            {            
                if ( (temppoint[j].x() == pointsSide_[i].x()))
                    bdouble = true;
            }
        }
        if (!bdouble)
        {
            temppoint.setSize(temppoint.size() + 1);
            temppoint[temppoint.size() - 1] = pointsSide_[i];
        }
    }
    pointsSide_ = temppoint;

    for (int i = 0; i < (pointsSide_.size() ); i++)
    {
        point temp(0, 0, 0);
        for (int j = 1; j < (pointsSide_.size() - i ); j++)
        {
            if(mType_ != Axial)
            {
                if (pointsSide_[j].z() < pointsSide_[j-1].z())
                {
                      temp = pointsSide_[j];
                      pointsSide_[j] = pointsSide_[j-1];
                      pointsSide_[j-1] = temp;
                }
            }
            else
            {
                if (pointsSide_[j].x() < pointsSide_[j-1].x())
                {
                    temp = pointsSide_[j];
                    pointsSide_[j] = pointsSide_[j-1];
                    pointsSide_[j-1] = temp;
                }   
            }
        }
    }    
}

void primitiveMixingPlanePatch::makeMaxRBPoints(const labelList& MaxRB)
{   
    bool bdebug = false;
    
    
    forAll(MaxRB, i)
    {
        maxRbPoints_.setSize(maxRbPoints_.size() + 1);
        maxRbPoints_[maxRbPoints_.size() - 1 ] = CCSBorderPoints_[BorderEdges_[ MaxRB[i] ][0] ];
        maxRbPoints_.setSize(maxRbPoints_.size() + 1);
        maxRbPoints_[maxRbPoints_.size() - 1 ] = CCSBorderPoints_[BorderEdges_[ MaxRB[i] ][1] ];
    }
    
    List<point> temppoint(0);    
    forAll(maxRbPoints_, i)
    {
        bool bdouble = false;
        forAll(temppoint, j)
        {
            if ( (temppoint[j].y() == maxRbPoints_[i].y()))
                bdouble = true;
        }
        if (!bdouble)
        {
            temppoint.setSize(temppoint.size() + 1);
            temppoint[temppoint.size() - 1] = maxRbPoints_[i];
        }
    }
    
    maxRbPoints_ = temppoint;
    
    if(bdebug)
    {
        OFstream os("MaxRBs.txt");
        os << maxRbPoints_ << endl;
        
    }
}

const List<point>& primitiveMixingPlanePatch::getPointsSide() const
{
    return pointsSide_;
}


void primitiveMixingPlanePatch::generateCCSTransformationTensors()
{
    if(tensorListList_.size() <  sliceFaces_.size())
    {
          forAll(sliceFaces_,i)
          {      
              
              tensorListList_.setSize(tensorListList_.size() + 1 );
              
              forAll(sliceFaces_[i], j)
              {                   
                 point faceCenter = AllFaces_[sliceFaces_[i][j]].centre(localPoints_);
               
                 tensor rotation = tensorFromCCSToCartesian(faceCenter);// makeRotationTensor(theta);
                // Info << rotation << endl;
                           
                 tensorListList_[i].setSize(tensorListList_[i].size() +1  );
                 tensorListList_[i][j] = rotation; 
              }
          }
   }
    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
