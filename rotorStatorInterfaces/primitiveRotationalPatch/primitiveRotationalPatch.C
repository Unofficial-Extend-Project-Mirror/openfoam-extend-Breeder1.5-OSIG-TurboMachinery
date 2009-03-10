/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

#include "primitiveRotationalPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "SubField.H"
#include "entry.H"
#include "OFstream.H"
#include "labelList.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(primitiveRotationalPatch, 0);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void primitiveRotationalPatch::calcMeshEdges() const
{
    
    
}

void primitiveRotationalPatch::FindCorners()
{   
    edgeList return_Value(4);
    
    List< Corner > angles;
    
    
    
    edgeList NewOrderedEdges = Reorderedges(); 
    
    BorderEdges_ = NewOrderedEdges;
    
//    Info << NewOrderedEdges.size() << endl;
//    OFstream os("NewEdges.obj");
//    forAll(NewOrderedEdges,iedges)
//    {
//        Info << NewOrderedEdges[iedges] << endl; 
//        writeOBJ(os, NewOrderedEdges[iedges].centre( BorderPoints_ ) );
//    }   
    
    int iCounter = 0;
    forAll(NewOrderedEdges,iEdges)
    {
        
        vector vecFirst = NewOrderedEdges[iEdges].vec(BorderPoints_);
        vector vecSecond;
        
        Corner Storage;
        Storage.FirstEdge = iEdges;
                 
        if(iEdges < NewOrderedEdges.size() - 2 )
        {
            vecSecond = NewOrderedEdges[iEdges+1].vec(BorderPoints_);
            Storage.SecondEdge = iEdges + 1;
        }
        else
        {
            vecSecond = NewOrderedEdges[0].vec(BorderPoints_);
            Storage.SecondEdge = NewOrderedEdges.size() -1;
        }
        // To reduce storage the sorting must be performed directly
        scalar angle = (vecSecond & vecFirst ) // Todo Warning old style Cast
                / (mag(vecSecond) * mag(vecFirst));
       
        Storage.Angle = angle; 
        Storage.ConnectionPoint = -1;
        angles.setSize(iCounter+1);
        angles[iCounter] = Storage;
        iCounter++; // to much why
    }
    
    // Eliminate all the values which have zero angle values
    //Mergesort is working fine
    
    List< Corner > anglesTemp(0);
    iCounter = 0;
    forAll( angles,iangles )
    {
       if( (angles[iangles].Angle * angles[iangles].Angle) != 1 )
       {  
          anglesTemp.setSize(iCounter + 1); 
          anglesTemp[iCounter++] = angles[iangles];          
       }
    }
   
    
    List< Corner > SortedCorners = mergesort(anglesTemp);
    
//    Info << "Border Points size " << BorderPoints_.size() << endl;
    
//    OFstream osC("NewCorners.obj");
//    forAll(SortedCorners,i)
//    {        
//        writeOBJ(osC, NewOrderedEdges[SortedCorners[i].FirstEdge].centre( BorderPoints_ ) );  
//        writeOBJ(osC, NewOrderedEdges[SortedCorners[i].SecondEdge].centre( BorderPoints_ ) );
//        writeOBJ(osC, BorderPoints_[ NewOrderedEdges[SortedCorners[i].FirstEdge][0]] );  
//        writeOBJ(osC, BorderPoints_[NewOrderedEdges[SortedCorners[i].FirstEdge][1]] );
//        writeOBJ(osC, BorderPoints_[ NewOrderedEdges[SortedCorners[i].SecondEdge][0]] );  
//        writeOBJ(osC, BorderPoints_[NewOrderedEdges[SortedCorners[i].SecondEdge][1]] );         
//    }
    
    
    //Until here everything is working perfectly
    
    // Delete doubles
    
    // Debugging must be made 
  
    
   /* 
    for (int i = 0; i < 7; ++i) 
    {        
        int j = 0;
        if(i%2)
        {
            return_Value[i] = BorderEdges_[SortedCorners[j].SecondEdge];
            j++;
        }
        else
        {
            return_Value[i] = BorderEdges_[SortedCorners[j].FirstEdge];
        }
    }*/
        
//    return_Value[0] = BorderEdges_[SortedCorners[0].FirstEdge];
//    return_Value[1] = BorderEdges_[SortedCorners[0].SecondEdge];
//    return_Value[2] = BorderEdges_[SortedCorners[1].FirstEdge];
//    return_Value[3] = BorderEdges_[SortedCorners[1].SecondEdge];
//    return_Value[4] = BorderEdges_[SortedCorners[2].FirstEdge];
//    return_Value[5] = BorderEdges_[SortedCorners[2].SecondEdge];
//    return_Value[6] = BorderEdges_[SortedCorners[3].FirstEdge];
//    return_Value[7] = BorderEdges_[SortedCorners[3].SecondEdge];
        
    
    // Debug only
    
     
//    forAll(SortedCorners, i)
//    {  
//        Info << SortedCorners[i].Angle 
//        << " EdgeFirst " << SortedCorners[i].FirstEdge
//        << " " << NewOrderedEdges[SortedCorners[i].FirstEdge]
//        << " SecondEdge " << SortedCorners[i].SecondEdge << " " << NewOrderedEdges[SortedCorners[i].SecondEdge] << endl;        
//    }
//    
    for (int i = 0; i < 4; i++)
    {
        if(SortedCorners.size() > i)
        {
            Corners_.setSize(i+1);
            Corners_[i] = SortedCorners[i];
            Corners_[i].ConnectionPoint = BorderEdges_[Corners_[i].FirstEdge][0];
        }
    }
    
//    Info << "Corners_;;;;;;;;;;;;;;;;;;;;;;;;:::::::::" << Corners_.size() << endl;
    
    OFstream osSC("SortedCorners.obj");
    forAll(Corners_, i)
    {
//        Info << "Corners_" << Corners_[i].FirstEdge << Corners_[i].SecondEdge << endl
//             << NewOrderedEdges[Corners_[i].FirstEdge].centre(BorderPoints_)  << " "
//             << NewOrderedEdges[Corners_[i].SecondEdge].centre(BorderPoints_)  << endl;
        
        writeOBJ(osSC,
                NewOrderedEdges[Corners_[i].FirstEdge].centre(BorderPoints_) );
        writeOBJ(osSC,
                NewOrderedEdges[Corners_[i].SecondEdge].centre(BorderPoints_) );
        writeOBJ(osSC,
                BorderPoints_[ NewOrderedEdges[Corners_[i].FirstEdge][0]]);
        writeOBJ(osSC,
                BorderPoints_[NewOrderedEdges[Corners_[i].FirstEdge][1]]);
        writeOBJ(osSC,
                BorderPoints_[ NewOrderedEdges[Corners_[i].SecondEdge][0]]);
        writeOBJ(osSC,
                BorderPoints_[NewOrderedEdges[Corners_[i].SecondEdge][1]]);
    }
     
    
    
}

void primitiveRotationalPatch::DefineLocalExternalEdges
(
    const edgeList& extEdges, 
    const pointField& extPoints
)
{
    BorderEdges_.setSize(extEdges.size());
    
    labelList pointsToExt(0);
    
    // Generates the necessary labels for extracting all the
    // points from the Original pointfield to a local pointfield
    forAll(extEdges, iEdges)
    {   
        label start = pointsToExt.size() -1;
        pointsToExt.setSize(extEdges[iEdges].size()+pointsToExt.size());
        pointsToExt[start+1] = extEdges[iEdges][0];
        pointsToExt[start+2] = extEdges[iEdges][1];        
    }
    
    // This checks wether a point is double and replaces that point by -1
    forAll(pointsToExt,ipoints)
    {
        forAll(pointsToExt,iSecond)
        {
            if(pointsToExt[iSecond] == pointsToExt[ipoints] && iSecond != ipoints )
                pointsToExt[iSecond] = -1;
            
        }
    } 
    
    // Removes all the double points from the list
    labelList TempPointsToExt;
    forAll(pointsToExt,ipoints)
    {
        if(pointsToExt[ipoints] > -1)
        {
            label newPos = TempPointsToExt.size();
            TempPointsToExt.setSize(TempPointsToExt.size()+1);
            TempPointsToExt[newPos] = pointsToExt[ipoints];
        }       
    }
    
    // Now the local pointfield is generated
    BorderPoints_.setSize(TempPointsToExt.size());
    forAll(TempPointsToExt,ipoint)
    {
        BorderPoints_[ipoint] = extPoints[ TempPointsToExt[ipoint] ];
    }   
    
    // Now the edges are generated
    forAll(extEdges,iEdge)
    {
        // label point one label point two
        label pointOne = -1;
        label pointTwo = -1;
        forAll(TempPointsToExt, ipoints )
        {
            if(TempPointsToExt[ipoints] == extEdges[iEdge][0] )
                pointOne = ipoints;
            
            if(TempPointsToExt[ipoints] == extEdges[iEdge][1] )
                pointTwo = ipoints;
        }
        
        //
        edge newEdge(pointOne, pointTwo);
        BorderEdges_[iEdge] = newEdge;       
    }   
    
    
}

edgeList primitiveRotationalPatch::Reorderedges()
{
    edgeList ReorderedList(0);
    
    
//    Info << "Size of BorderEdges " << BorderEdges_.size() << endl;
    if (BorderEdges_.size() > 0 )
    {
       // Info << BorderEdges_ << endl;
        ReorderedList.setSize(ReorderedList.size() +1);
        ReorderedList[0] = BorderEdges_[0];   
    }
    
    for (int iCounter = 0; iCounter < BorderEdges_.size()-1; iCounter++) 
    {
        forAll(BorderEdges_, iEdges)
        {
            if (ReorderedList[iCounter][0] == BorderEdges_[iEdges][1] 
                  && ( ReorderedList[iCounter] != BorderEdges_[iEdges]) )
            {
                ReorderedList.setSize(ReorderedList.size() + 1 );
                ReorderedList[iCounter + 1] = BorderEdges_[iEdges];
            }
        }
    }
     
    return (ReorderedList) ;
}

edgeList primitiveRotationalPatch::OrigExtEdges(const primitivePatch& rOrigPatch )
{
    edgeList return_Value;
    edgeList OrigEdges = rOrigPatch.edges();
            
    forAll(OrigEdges, iedges)
    {       
        if ( !rOrigPatch.isInternalEdge(iedges) )
        { 
            label next = return_Value.size();
            return_Value.setSize(return_Value.size() + 1);
            return_Value[next] = OrigEdges[iedges]; 
        }
    }
    
    return (return_Value);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void primitiveRotationalPatch::writeOBJ(Ostream& os, const point& pt)
{
    os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

primitiveRotationalPatch::primitiveRotationalPatch
(
   const List<face>& faces, 
   const pointField& points,
   const primitivePatch& OrigPatch,
   const word& name
)
: 
primitiveAuxPatch 
(
   List<face>(1),
   pointField(1)
),
name_(name),
AllFaces_(0),
localPoints_(0),
OrigPatch_( OrigPatch),
//PointsCorners_(PointsCorners),
CCSBorderPoints_(0),
Corners_(0),
pointsSide_(0),
maxRbPoints_(0),
primPatch_(NULL)
{    
     cylindricalCS_ = cylindricalCS( 
               "ccs", 
               vector(0, 0, 0.005),
               vector(0, 0, 1),  
               vector(1, 0, 0)
      );
     
     mType_ = determineMachine();
     
     
     
     //calcRadialGeometry(); 
}  
   
primitiveRotationalPatch::primitiveRotationalPatch
(       
         const List<face>& faces, 
         const pointField& points,
         const primitivePatch& OrigPatch,
         vector CCSCenter,
         vector CCSaxis,
         vector CCSdirection,
         const word& name
)
: 
primitiveAuxPatch 
(
   faces,
   points
),
name_(name),
AllFaces_(0),
localPoints_(0),
OrigPatch_(OrigPatch),
CCSBorderPoints_(0),
Corners_(0),
pointsSide_(0),
maxRbPoints_(0),
primPatch_(NULL)
{
    Info << "Center"     << CCSCenter <<  endl
         << "Direction " << CCSdirection <<  endl
         << "Axis "      << CCSaxis <<  endl;
    
    cylindricalCS_ = cylindricalCS
    (
            "ccs", 
            CCSCenter, 
            CCSaxis,
            CCSdirection 
    );

    mType_ = determineMachine();
    
   // calcRadialGeometry();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

primitiveRotationalPatch::~primitiveRotationalPatch()
{
    //clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void primitiveRotationalPatch::clearAddressing()
{
    
    
}


void primitiveRotationalPatch::movePoints(const pointField& p)
{
    primitiveAuxPatch::movePoints(p);
}


const cylindricalCS& primitiveRotationalPatch::getCylindricalCoordinateSystem() const
{
    return (cylindricalCS_);    
}

void primitiveRotationalPatch::SetCylindricalCoordinateSystem(const cylindricalCS& ccs )
{
    cylindricalCS_ = cylindricalCS(ccs.name(),ccs.origin(),ccs.axis(),ccs.direction());
    calcRadialGeometry();
}

//const pointField& primitiveRotationalPatch::GetFacePoints() const
//{
//    return (cylindricalCS_.globalPosition( localPoints_ ));
//}
//
//void primitiveRotationalPatch::GetFacePoints(pointField& PF) const
//{
//    PF.setSize(localPoints_.size());
//    PF = cylindricalCS_.globalPosition( localPoints_ );
//}
//
//void primitiveRotationalPatch::GetFacePoints(List<face>& Faces) const
//{
//    Faces.setSize(AllFaces_.size()); 
//    forAll(AllFaces_,i)
//    {
//        Faces[i] = AllFaces_[i];
//    }   
//}

//const List<face>& primitiveRotationalPatch::GetFaces() const
//{
//    return ( AllFaces_ );
//}

vector primitiveRotationalPatch::transformFromCCSToCartesian(vector& ccsVec, point& Position, point& origin ) const
{
  point newPoint = Position-origin;  
    
  scalar sinus   =  newPoint.y()/(Foam::sqrt( scalar(newPoint.y() * newPoint.y() + newPoint.x() * newPoint.x() ) ) ) ; 
  scalar cosinus = newPoint.x()/(Foam::sqrt( scalar(newPoint.y() * newPoint.y() + newPoint.x() * newPoint.x() ) ) ) ;
 // Info << "cosinus" << cosinus << endl;
 // Info << "sinus" << sinus << endl;
  vector newVec = vector(0,0,0);
  newVec.x() = ccsVec.x()*cosinus - ccsVec.y() * sinus;
  newVec.y() = ccsVec.x()*sinus + ccsVec.y() * cosinus;
    
  return  (newVec);  
}

vector primitiveRotationalPatch::transformFromCCSToCartesian( vector& ccsVec, point& Position) const
{   
    
  point newPoint = cylindricalCS_.localPosition(Position);  
  makeRotationTensor(newPoint.y());  
  
  
  scalar sinus   =  newPoint.y()/(Foam::sqrt( scalar(newPoint.y() * newPoint.y() + newPoint.x() * newPoint.x() ) ) ) ; 
  scalar cosinus = newPoint.x()/(Foam::sqrt( scalar(newPoint.y() * newPoint.y() + newPoint.x() * newPoint.x() ) ) ) ;
  Info << "cosinus" << cosinus << endl;
  Info << "sinus" << sinus << endl;
  vector newVec = vector(0,0,0);
  newVec.x() = ccsVec.x()*cosinus - ccsVec.y() * sinus;
  newVec.y() = ccsVec.x()*sinus + ccsVec.y() * cosinus;
    
  return  (newVec);  
}

tensor primitiveRotationalPatch::tensorFromCCSToCartesian( point& Position ) const 
{
    point newPoint = cylindricalCS_.localPosition(Position);  
    
    tensor newTensor = makeRotationTensor(newPoint.y());
     
    return (newTensor );  
}

tensor primitiveRotationalPatch::tensorFromCartesianToCCS( point& Position ) const
{
    point newPoint = Position - cylindricalCS_.origin();  
    
    tensor newTensor = makeRotationTensor(newPoint.y());
    
    return (newTensor.T());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "calcRotationalGeometry.C"
#   include "sortFunctionsRotationalPatch.C"

// ************************************************************************* //


