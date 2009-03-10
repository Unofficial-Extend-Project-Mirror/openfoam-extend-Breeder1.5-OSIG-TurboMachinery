/*---------------------------------------------------------------------------*\
  =========                 |
 /  F ield         | OpenFOAM: The Open Source CFD Toolbox
 /   O peration     |
 /    A nd           | Copyright held by original author
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
#include "primitiveRotationalPatch.H"
#include "face.H"
#include "mathematicalConstants.H"
#include "coordinateSystem.H"
#include "triPointRef.H"

namespace Foam

{

// pure virtual function
void primitiveRotationalPatch::calcRadialGeometry()
{    
}

void primitiveRotationalPatch::getBorders(labelList& MinRB, labelList& MaxRB,
        labelList& MinCB, labelList& MaxCB)
{
    bool bdebug = false;
    //    FindCorners();   

    // Finding which two corners have the highest Radial values
    // Using the bubblesort seems to be appropriate because the sorted list
    // is very small    
    //Corner temp;
    label temp;

    if (bdebug)
        Info << "Begin Function primitiveRotationalPatch::getBorders " << endl;

    labelList labelCorners(Corners_.size());

    forAll(Corners_, i)
    {
        
        if (bdebug)
            Info << Corners_[i].FirstEdge << endl;

        if (Corners_[i].FirstEdge > Corners_[i].SecondEdge)
            labelCorners[i] = Corners_[i].FirstEdge;
        else
            labelCorners[i] = Corners_[i].SecondEdge;
    }

    if (false)
    {
        Info << labelCorners << endl;
        Info << BorderEdges_ << endl;
    }

    labelListList tempOne;

    labelListList Borders(0);

    // now one Corner is used as a starting point the algorithm searches for 
    // minimal label difference between the corner labels. 

    // 080417
    // No more sorting any more. 
    // The idea is to start at the first edge of the
    // first corner and walk through all the following 
    // edges until there is a hit with another corners second edge.

    OFstream oI("j.txt");
    forAll(labelCorners, i)
    {
        bool bhit = false;
        Borders.setSize(Borders.size() + 1);
        int j = labelCorners[i];
        while (!bhit)
        {
            if (j == BorderEdges_.size())
                j = 0;

            forAll(Corners_, k)
            {
                if ( (Corners_[k].FirstEdge == j) && (k != j)
                        && (labelCorners[i] != Corners_[k].FirstEdge ))
                {
                    bhit = true;
                    oI << "k " << k << " j " << j << " i " << labelCorners[i]
                            << endl;
                }

                if ( (Corners_[k].SecondEdge == j) && (k != j)
                        && (labelCorners[i] != Corners_[k].SecondEdge ))
                {
                    bhit = true;
                    oI << "k " << k << " j " << j << " i " << labelCorners[i]
                            << endl;
                }
            }
            Borders[i].setSize(Borders[i].size() +1);
            Borders[i][Borders[i].size()-1] = j;
            j++;
        }
    }

    if (bdebug)
    {
        OFstream os("test.txt");
        OFstream osOBJ("testBorders.obj");
        forAll(Borders, i)
        {
            os << Borders << endl;
            forAll(Borders[i], j)
            {
                writeOBJ(osOBJ,
                        (CCSBorderPoints_[ BorderEdges_[ Borders[i][j] ][0]]));
            }
        }
    }  

    // Now the borders can be sorted for the radial and circumferential 
    // coordinate component. This is done by making a average point 
    // over all the edge points
    List<point> AveragePoints(Borders.size() );
    forAll(Borders, i)
    {
        AveragePoints[i] = point(0,0,0);
        int iCounter = 0;
        forAll(Borders[i], j)
        {
            iCounter++;
            AveragePoints[i]
                    = CCSBorderPoints_[BorderEdges_[ Borders[i][j]][0] ]
                            + AveragePoints[i];
            iCounter++;

            if (false)
            {
                if (i == Borders.size()-1)
                    Info << "CCS For Test"
                            << CCSBorderPoints_[BorderEdges_[ Borders[i][j]][1] ]
                            << endl << " Point Cart "
                            << BorderPoints_[BorderEdges_[ Borders[i][j]][1] ];
            }
            AveragePoints[i]
                    = CCSBorderPoints_[BorderEdges_[ Borders[i][j]][1] ]
                            + AveragePoints[i];
        }
        if (iCounter)
            AveragePoints[i] = AveragePoints[i] / iCounter;
    }
    if (bdebug)
    {
        Info << "AveragePoints " << AveragePoints << endl;

        Info << " CCSBorderPoints " << CCSBorderPoints_ << endl;
    }

    // First Search for the min and max Radial components of the
    // average Points
    labelList avlabels(AveragePoints.size() );
    forAll(avlabels, i)
    {
        avlabels[i] = i;
    }

    temp = -1;

    
    //Sorting for the radial direction
    if(mType_ != Axial)
    {
        for (int i = 0; i < (AveragePoints.size() - 1); i++)
        {
            for (int j = 1; j < AveragePoints.size(); j++)
            {
                if (AveragePoints[avlabels[j] ].z()
                        > AveragePoints[avlabels[j-1] ].z() )
                {
                    temp = avlabels[j];
                    avlabels[j] = avlabels[j-1];
                    avlabels[j-1] = temp;
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < (AveragePoints.size() - 1); i++)
        {
             for (int j = 1; j < AveragePoints.size(); j++)
             {
                 if (AveragePoints[avlabels[j] ].x()
                      > AveragePoints[avlabels[j-1] ].x() )
                 {
                     temp = avlabels[j];
                     avlabels[j] = avlabels[j-1];
                     avlabels[j-1] = temp;
                 }
             }
        }
    }

    MinRB = Borders[avlabels[avlabels.size()-1 ] ];
    MaxRB = Borders[avlabels[0] ];

    
    OFstream oAVLAb("AVLAB.txt");
    forAll(avlabels, i)
    {
        oAVLAb << AveragePoints[avlabels[i]] << endl;
        
    }
    
    labelList circLabels(2);

    circLabels[0] = avlabels[1];
    circLabels[1] = avlabels[2];
    
    OFstream osAV("avlabels.txt");
    forAll(circLabels,j)
    {
        osAV << AveragePoints[circLabels[j] ].y() << endl;
         
    }
    
    if(AveragePoints[circLabels[1] ].y() > AveragePoints[circLabels[0] ].y())
    {
        MaxCB = Borders[circLabels[1]];
        MinCB = Borders[circLabels[0]];
        Info << "MinCB Size " << MinCB << endl;
        Info << Borders[circLabels[0]] << endl;
    }
    else
    {
        MaxCB = Borders[circLabels[0]];
        MinCB = Borders[circLabels[1]];
	Info << "MinCB Size " << MinCB << endl;
        Info << Borders[circLabels[0]] << endl;
//         MinCB.setSize(MinCB.size()+1);
// 	Info << "MinCB Size " << MinCB << endl;
    }
    

    if (bdebug)
    {
        OFstream os("bordersMaxCB.obj");
        forAll(MaxCB, i)
        {
            writeOBJ(os, (BorderPoints_[ BorderEdges_[ MaxCB[i] ][0]]));
            writeOBJ(os, (BorderPoints_[ BorderEdges_[ MaxCB[i] ][1]]));
        }

        OFstream osMinCB("bordersMinCB.obj");
        forAll(MinCB, i)
        {
            if(i != (MinCB.size()-1) )
                writeOBJ(osMinCB, (BorderPoints_[ BorderEdges_[ MinCB[i] ][0]]));
            
            writeOBJ(osMinCB, (BorderPoints_[ BorderEdges_[ MinCB[i] ][1]]));
        }

        OFstream osMaxRB("bordersMaxRB.obj");
        forAll(MaxRB, i)
        {
            writeOBJ(osMaxRB, (BorderPoints_[ BorderEdges_[ MaxRB[i] ][0]]));
            writeOBJ(osMaxRB, (BorderPoints_[ BorderEdges_[ MaxRB[i] ][1]]));
        }

        OFstream osMinRB("bordersMinRB.obj");
        forAll(MinRB, i)
        {
            writeOBJ(osMinRB, (BorderPoints_[ BorderEdges_[ MinRB[i] ][0]]));
            writeOBJ(osMinRB, (BorderPoints_[ BorderEdges_[ MinRB[i] ][1]]));
        }

    }    
    
    
    labelList tempList(0);
    bool changed = false;

 Info << "Before sorting " << MinCB << endl;



    forAll(Corners_,k)
    {
        if( MinCB.size() > 2 )
        {
	  if( MinCB[0] == Corners_[k].FirstEdge  && MinCB[1] == Corners_[k].SecondEdge
                || MinCB[0] == Corners_[k].FirstEdge  && MinCB[1] == Corners_[k].SecondEdge )
	  {
            forAll(MinCB,j)
            {
                if( j > 0 )
                {
                    tempList.setSize(tempList.size() + 1);
                    tempList[j-1] = MinCB[j] ;
                }   
            }
            changed = true;
	  }
        
	  if( MinCB[MinCB.size() - 2] == Corners_[k].FirstEdge  && MinCB[MinCB.size() - 1] == Corners_[k].SecondEdge
	    || MinCB[MinCB.size() - 2] == Corners_[k].FirstEdge  && MinCB[MinCB.size() - 2] == Corners_[k].SecondEdge )
                {
                    forAll(MinCB,j)
                    {
                        if( j < ( MinCB.size() - 1) )
                        {
                            tempList.setSize(tempList.size() + 1);
                            tempList[j] = MinCB[j] ;
                        }
                        
                    }
                    changed = true;
	  }
        }
    }
    
    if(changed)
        MinCB = tempList;

    Info << "MinCB after sorting " << MinCB << endl;
    
    
    // Delete Points of Corners
    if (bdebug)
        Info << "End Function primitiveRotationalPatch::getBorders " << endl;
}

primitiveRotationalPatch::MachineType primitiveRotationalPatch::determineMachine()
{
    bool bdebug = false;
    MachineType mtype = Axial;

    vectorField faceNormals = OrigPatch_.faceNormals();

    // Average over all the face normals
    vector AverageNormal(0, 0, 0);

    forAll(faceNormals, iNormals)
    {
        AverageNormal = AverageNormal + faceNormals[iNormals];
    }

    if(faceNormals.size() > 0)
        AverageNormal = AverageNormal / faceNormals.size();

    if (sqrt(sqr(AverageNormal.x()) ) <= SMALL)
        AverageNormal.x() = 0;

    if (sqrt(sqr(AverageNormal.y()) ) <= SMALL)
        AverageNormal.y() = 0;

    if (sqrt(sqr(AverageNormal.z()) ) <= SMALL)
        AverageNormal.z() = 0;

    if (bdebug)
        Info << AverageNormal << endl;

    AverageNormal = cylindricalCS_.localVector(AverageNormal);

    if (bdebug)
        Info << AverageNormal << endl;
        Info <<    sqrt( AverageNormal.x()*AverageNormal.x() )  << endl;
    
    // Todo must be using abs and eps
    if ( ( sqrt( AverageNormal.x()*AverageNormal.x() )   <= SMALL) && ( 
            sqrt( AverageNormal.z()*AverageNormal.z() ) > SMALL ) )
    {
        mtype = Axial;
        if (bdebug)
            Info << "Axial" << endl;
    } 
    
    if (sqrt( AverageNormal.z()*AverageNormal.z() ) <= SMALL  )
    {
        mtype = Radial;
        if (bdebug)
            Info << "Radial" << endl;
    } 

    if ( ( sqrt( AverageNormal.x()*AverageNormal.x() )   > 0.0001) && ( 
            sqrt( AverageNormal.z()*AverageNormal.z() ) > SMALL ) )
    {
        mtype = Diagonal;
        if (bdebug)
            Info << "Diagonal" << endl;
    }

    return mtype;
}

List<point> primitiveRotationalPatch::calculateForeignMatchedPoints
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

void primitiveRotationalPatch::generateCCSBorderPoints()
{
    CCSBorderPoints_.setSize(BorderPoints_.size());

    CCSBorderPoints_ = cylindricalCS_.localPosition(BorderPoints_);
}

void primitiveRotationalPatch::generatePointList
(
   const labelList& MaxRB, 
   const labelList& MinCB
)
{    
//    List<List<point> > orderedPoints;
//    bool bdebug = false;
//    bool bCSV = false;
//
//    // Finding the starting point for the upper points
//    // That is done by finding the matching points between the MaxRB and 
//    // MinCB because the      
//    
//    makeSidePoints(MinCB);
//    
//    if(bdebug)
//    {
//        OFstream opS("pointsSide.txt");
//        opS << pointsSide_ << endl;
//    }
//  
//    makeMaxRBPoints(MaxRB);
//
//    //      After the pointsSide are generated 
//    //      They can be easily transformed by tensor operations
//    //      The necessary delta angles are extracted by from the MaxRB 
//
//    OFstream osS("orderedPoints.csv");
//    
//    OFstream osTest("maxRBallPoints.obj");
//
//    forAll(maxRbPoints_, j)
//    {
//        writeOBJ(osTest, cylindricalCS_.globalPosition(maxRbPoints_[j]));
//
//    }
//
//    forAll(pointsSide_, j)
//    {
//        writeOBJ(osTest, cylindricalCS_.globalPosition(pointsSide_[j]));
//
//    }
//    // Until here it looks well
//
//    // Now the angles are extracted to a list
//    //      osS << "RB Points size " << maxRbPoints_.size() << endl;
//    List<scalar> angleList(0);
//    OFstream osAngles("angles.txt");
//   
//    forAll(maxRbPoints_, i)
//    {
//        angleList.setSize(angleList.size() +1);
//    
//       osAngles << maxRbPoints_[i].y() << endl;
//        
//        angleList[i] = maxRbPoints_[i].y();
//    }
//    angleList = mergesort( angleList );
////    osAngles << "Angles Sorted " << endl;
//    forAll(angleList,i)
//    {
//        osAngles << angleList[i] << endl;
//    }
//    
//    
////    osAngles << "MinCB Y" << endl;
////    forAll(pointsSide_, i)
////    {
////        osAngles << pointsSide_[i].y() << endl;       
////    }
//
//
//
//    OFstream oLast("LastTest.txt");
//    forAll(pointsSide_, i)
//    {
//        //         osS << "in For pointsSide" << endl;
//        orderedPoints.setSize(orderedPoints.size() + 1);
//
//        orderedPoints[i].setSize(orderedPoints[i].size() + 1);
//        orderedPoints[i][orderedPoints[i].size() -1 ] = pointsSide_[i];
//
//        oLast << "Side point" << pointsSide_[i] << endl;
//        forAll(angleList, j)
//        {
//            if (j > 0)
//            {
//                point temp(pointsSide_[i].x(), angleList[j], pointsSide_[i].z());
//                oLast <<  temp << endl;
//                
//                orderedPoints[i].setSize(orderedPoints[i].size() + 1);
//                orderedPoints[i][orderedPoints[i].size() -1 ] = temp;
//
////                osAngles << "J " << j << " Angle " << angleList[j] << endl;
//                //        osS << "in if" << endl;
//            }
//            //      osS << "in For angleList" << endl;
//        }
//    }
//    // osS << angleList.size() << endl;
//    OFstream osP("orderedPoints.obj");
//    
//    
//    oLast << CCSBorderPoints_ << endl;
//
//    forAll(orderedPoints, i)
//    {
//        //          osS << orderedPoints[i].size() << endl;
//    }
//    
//    if (true)
//    {
//        forAll(orderedPoints, i)
//        {   
//            forAll(orderedPoints[i], j)
//            {
//                writeOBJ(osP,
//                        cylindricalCS_.globalPosition(orderedPoints[i][j]));
//                    osS << orderedPoints[i][j].x() << ", "
//                        << orderedPoints[i][j].y() << ", "
//                        << orderedPoints[i][j].z() << endl;
//            }
//        }
//    }
//    
//    writeOBJ(osP,
//            cylindricalCS_.origin());
//    
//    OFstream osCentre("centre.obj");
//    writeOBJ(osCentre,
//               cylindricalCS_.origin());
//
//    // Generating the pointfield and the labels
//    label Counter = -1;
//    localPoints_.setSize(0);
//    forAll(orderedPoints, i)
//    {
//        pointLabels_.setSize(i + 1);
//        forAll(orderedPoints[i], j)
//        {
//            Counter++;
//            pointLabels_[i].setSize(j + 1);
//            pointLabels_[i][j] = Counter;
//
//            localPoints_.setSize(localPoints_.size() + 1);
//            localPoints_[localPoints_.size() - 1] = cylindricalCS_.globalPosition(orderedPoints[i][j] );
//        }
//    }
//    
//    OFstream of("localPoints.obj");
//    forAll(localPoints_,i)
//    {
//        writeOBJ(of,localPoints_[i]);
//        
//    }
//
//    OFstream os("test.obj");
//    forAll(orderedPoints, i)
//    {
//        forAll(orderedPoints[i], j)
//        {
//            writeOBJ(os, cylindricalCS_.globalPosition(orderedPoints[i][j]) );
//            if (bCSV)
//                Info << cylindricalCS_.globalVector( orderedPoints[i][j] ).x() << ", "
//                        << cylindricalCS_.globalVector( orderedPoints[i][j] ).y() << ", "
//                        << cylindricalCS_.globalVector( orderedPoints[i][j] ).z() << endl;
//        }
//    }
//
//    writeOBJ(os, cylindricalCS_.origin() );
//
//    if (bdebug)
//        Info << " PointsInternal " << orderedPoints << endl;
//
//    if (bdebug)
//        Info << " BorderPoints CS " << CCSBorderPoints_ << endl;
//    

}

void primitiveRotationalPatch::generateNewFaceList()
{
//    bool bdebug = false;
//    label counterfaces = -1;
//
//    sliceFaces_.setSize(0);
//    AllFaces_.setSize(0);
//
//    // Patch faceNormals gives back a normalized normal.
//    vectorField faceNormals = OrigPatch_.faceNormals();
//
//    // Average over all the face normals
//    vector AverageNormal(0, 0, 0);
//
//    forAll(faceNormals, iNormals)
//    {
//        AverageNormal = AverageNormal + faceNormals[iNormals];
//    }
//
//    AverageNormal = AverageNormal / faceNormals.size();
//    
//    
//    OFstream os("Mags.txt");
//
//    if (mag( AverageNormal.x() ) < SMALL)
//        AverageNormal.x() = 0;
//
//    if ( mag(AverageNormal.y()) < SMALL)
//        AverageNormal.y() = 0;
//
//    if ( mag(AverageNormal.z() ) < SMALL)
//        AverageNormal.z() = 0;
//    os << "Original Normal " << name_<< " " << AverageNormal << endl;
//    
//    forAll(pointLabels_, i)
//    {
//
//        labelList llP(3);
//
//        if (i == (pointLabels_.size()-1))
//            return;
//
//        sliceFaces_.setSize(sliceFaces_.size() + 1);
//
//        if (bdebug)
//            Info << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!In if "
//                    << pointLabels_[i].size()-1 << endl;
//
//        for (label j = 0; j < pointLabels_[i].size()-1; j++)
//        {   
//            vector FirstNormal = triPointRef
//            (
//                    localPoints_[pointLabels_[i][j+1]],
//                    localPoints_[pointLabels_[i][j]] ,
//                    localPoints_[pointLabels_[i+1][j]] 
//            ).normal();
//
//            Info << pointLabels_[i][j+1] << " "
//                 << pointLabels_[i][j] << " "
//                 << pointLabels_[i+1][j] << endl; 
//            
//            Info << FirstNormal << endl;
//            if( (mag(FirstNormal) > 0 )|| (mag(FirstNormal) < 0 ) ) 
//                FirstNormal /= mag(FirstNormal);
//           
//            if (mag(FirstNormal + AverageNormal) < SMALL)
//            {             
//                llP[0] = pointLabels_[i][j+1];
//                llP[1] = pointLabels_[i][j];
//                llP[2] = pointLabels_[i+1][j];
//            } else
//            {
//                llP[1] = pointLabels_[i][j+1];
//                llP[0] = pointLabels_[i][j];
//                llP[2] = pointLabels_[i+1][j];
//
//            }
//            AllFaces_.setSize(AllFaces_.size() + 1);
//            AllFaces_[AllFaces_.size() - 1] = face(llP);
//
//            sliceFaces_[i].setSize(sliceFaces_[i].size() + 1);
//            sliceFaces_[i][sliceFaces_[i].size() - 1] = AllFaces_.size()-1;
//
//            FirstNormal = triPointRef
//            (
//                    localPoints_[pointLabels_[i][j+1]],
//                    localPoints_[pointLabels_[i+1][j]],
//                    localPoints_[pointLabels_[i+1][j+1]] 
//            ).normal();
//
//
//            FirstNormal /= mag(FirstNormal);
//
//            if (mag(FirstNormal + AverageNormal) < SMALL)
//            {
//                llP[0] = pointLabels_[i][j+1];
//                llP[1] = pointLabels_[i+1][j];
//                llP[2] = pointLabels_[i+1][j+1];
//            } else
//            {
//                llP[1] = pointLabels_[i][j+1];
//                llP[0] = pointLabels_[i+1][j];
//                llP[2] = pointLabels_[i+1][j+1];
//            }
//
//            
//                os << "Mag " << mag(FirstNormal - AverageNormal)
//                        << " FirsNormal " << FirstNormal << "AverageNormal "
//                        << AverageNormal << endl;
//
//            AllFaces_.setSize(AllFaces_.size() + 1);
//            AllFaces_[AllFaces_.size() - 1] = face(llP);
//
//            counterfaces++;
//            sliceFaces_[i].setSize(sliceFaces_[i].size() + 1);
//            sliceFaces_[i][sliceFaces_[i].size() - 1] = AllFaces_.size()-1;
//        }
//    }    
}

void primitiveRotationalPatch::makeSidePoints( const labelList& MinCB )
{   
    forAll(MinCB, i)
      {
          pointsSide_.setSize(pointsSide_.size() + 1);
          pointsSide_[pointsSide_.size() - 1 ] = CCSBorderPoints_[BorderEdges_[ MinCB[i] ][0] ];
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

void primitiveRotationalPatch::makeMaxRBPoints(const labelList& MaxRB)
{   
//    bool bdebug = false;
//    
//    
//    forAll(MaxRB, i)
//    {
//        maxRbPoints_.setSize(maxRbPoints_.size() + 1);
//        maxRbPoints_[maxRbPoints_.size() - 1 ] = CCSBorderPoints_[BorderEdges_[ MaxRB[i] ][0] ];
//        maxRbPoints_.setSize(maxRbPoints_.size() + 1);
//        maxRbPoints_[maxRbPoints_.size() - 1 ] = CCSBorderPoints_[BorderEdges_[ MaxRB[i] ][1] ];
//    }
//    
//    List<point> temppoint(0);    
//    forAll(maxRbPoints_, i)
//    {
//        bool bdouble = false;
//        forAll(temppoint, j)
//        {
//            if ( (temppoint[j].y() == maxRbPoints_[i].y()))
//                bdouble = true;
//        }
//        if (!bdouble)
//        {
//            temppoint.setSize(temppoint.size() + 1);
//            temppoint[temppoint.size() - 1] = maxRbPoints_[i];
//        }
//    }
//    
//    maxRbPoints_ = temppoint;
//    
//    if(bdebug)
//    {
//        OFstream os("MaxRBs.txt");
//        os << maxRbPoints_ << endl;
//        
//    }
}

tensor primitiveRotationalPatch::makeRotationTensor(const scalar& angle) const
{
    // it has to be assured that the axis is normalized
    vector omega = cylindricalCS_.axis()/mag(cylindricalCS_.axis());
       
    // sin and cos only use radians
    //Info << "Theta" << angle << endl;
    scalar theta = ( angle   ) / 180.0 * mathematicalConstant::pi;            
    //Info << theta << endl;
    tensor omegaHat(tensor::zero);
    //Rodriguez Formula see http://mathworld.wolfram.com/RodriguesRotationFormula.html
    // skew tensor form vector
    omegaHat.replace(tensor::XY,-omega.z() );
    omegaHat.replace(tensor::XZ, omega.y() );
    omegaHat.replace(tensor::YZ,-omega.x() );
    omegaHat = (omegaHat - omegaHat.T());                 
           
//    Info << sin(theta) << endl
//    << omegaHat * (sin(theta)) << endl
//    << (omegaHat & omegaHat)*(1.0 - cos(theta)) << endl;
    tensor rotation = I + omegaHat * (sin(theta)) + (omegaHat & omegaHat)*(1.0 - cos(theta));
  // Info << rotation << endl;       
    return (rotation);     
    
}
       
tensor primitiveRotationalPatch::makeRotationTensor(vector& axis, scalar& angle)
{   
    // it has to be assured that the axis is normalized
    vector omega = axis/mag(axis);
    
    // sin and cos only use radians 
    scalar theta = ( angle   ) / 180.0 * mathematicalConstant::pi;            
 
    tensor omegaHat(tensor::zero);
    //Rodriguez Formula see http://mathworld.wolfram.com/RodriguesRotationFormula.html
    // skew tensor form vector
    omegaHat.replace(tensor::XY,-omega.z() );
    omegaHat.replace(tensor::XZ, omega.y() );
    omegaHat.replace(tensor::YZ,-omega.x() );
    omegaHat = (omegaHat - omegaHat.T());                 
                  
    tensor rotation = I + omegaHat * (sin(theta)) + (omegaHat & omegaHat)*(1.0 - cos(theta));
        
    return (rotation); 
}
       


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
