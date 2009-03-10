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

#include "domainScalingPolyPatch.H"
#include "mathematicalConstants.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "demandDrivenData.H"
#include "polyPatchID.H"
#include "SubField.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(domainScalingPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, domainScalingPolyPatch, word);
  //  addToRunTimeSelectionTable(polyPatch, domainScalingPolyPatch, Istream);
    addToRunTimeSelectionTable(polyPatch, domainScalingPolyPatch, dictionary);
    
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



void Foam::domainScalingPolyPatch::clearOut()
{   
    if(reconFaceCellCentresPtr_)
        deleteDemandDrivenData(reconFaceCellCentresPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components


Foam::domainScalingPolyPatch::domainScalingPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
   /* const word& masterName,
    const word& shadow_Name,    
    const vector& origin,
    const vector& axis,
    const vector& direction*/
)
:
    coupledPolyPatch(name, size, start, index, bm),
    shadow_("master"),       
    master_("slave"),  // Because at that point the geometry of the original patch has not been build
    // this constructor is usually called by the blockMesh therefore the reconstructing of the cellcentres must be 
    // avoided this is done by setting the master equal slave
    origin_(point(0,0,0)),
    axis_(vector(0,0,0)),
    direction_(vector(0,0,0)),
    transformedFaceCellCentresPtr_(NULL),
    masterTransformedPrimitivePointsPtr_(NULL),
    rotPatch_(NULL),  
//    relative_(false),
    reconFaceCellCentresPtr_(NULL),
    transformedPrimitivePatchPtr_(NULL)
    
{        
}

/*// Construct from Istream
Foam::domainScalingPolyPatch::domainScalingPolyPatch
(
    Istream& is,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(is, index, bm),
    shadow__(is),
    shadowIndex_(-1),
    transformPrimitivePatchPtr_(NULL),
    transformFaceCellCentresPtr_(NULL),
    patchToPatchPtr_(NULL),
    reconFaceCellCentresPtr_(NULL)
{}
*/

// Construct from dictionary
Foam::domainScalingPolyPatch::domainScalingPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(name, dict, index, bm),
    shadow_(dict.lookup("shadowPatch")),    
    master_(dict.lookup("hierarchy")),
    origin_(dict.lookup("origin")),
    axis_(dict.lookup("axis")),
    direction_(dict.lookup("direction")),
    shadowIndex_(-1), 
    transformedFaceCellCentresPtr_(NULL),
    masterTransformedPrimitivePointsPtr_(NULL),
    rotPatch_(NULL),
    reconFaceCellCentresPtr_(NULL),
//    relative_(false),
    omega_(readScalar(dict.lookup("omega"))),
    start_(readLabel(dict.lookup("startFace")) ),
    nFaces_( readLabel(dict.lookup("nFaces")))

//    (
//            
//         faceSubList( bm.mesh().allFaces(),
//                      readLabel(dict.lookup("nFaces")),
//                      readLabel(dict.lookup("startFace"))
//                    ),
//          bm.mesh().allPoints(),
//          *this,
//          dict.lookup("origin"),
//          dict.lookup("axis"),
//          dict.lookup("direction"),
//          readScalar(dict.lookup("omega")),
//          name,
//          false
//     )            
{
        Info << "Construct from dictionary " << endl;
        if( dict.found("startAngle") )
        {
             startAngle_ = scalar
             (
                 readScalar
                 ( 
                     dict.lookup("startAngle") 
                 )
              );
        }
        
        
//        if( dict.found("relative") )
//        {
//            relative_ = Switch(
//                                dict.lookup("relative") 
//                              );
//        }        
        
}


//- Construct as copy, resetting the boundary mesh
Foam::domainScalingPolyPatch::domainScalingPolyPatch
(
    const domainScalingPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    shadow_(pp.shadow_),       
    master_(pp.master_),    
    shadowIndex_(-1),
//    relative_(false),    
    transformedFaceCellCentresPtr_(NULL), 
    masterTransformedPrimitivePointsPtr_(NULL),
    reconFaceCellCentresPtr_(NULL),
    transformedPrimitivePatchPtr_(NULL)  
   
   
{
        Info << " Construct as copy, resetting the boundary mesh" << endl;
        rotPatch_= new primitiveDomainScalingPatch(
                       (*this), 
                       (*this).points(),
                       (*this),
                       pp.RotPatch().getCylindricalCoordinateSystem().origin(),
                       pp.RotPatch().getCylindricalCoordinateSystem().axis(),
                       pp.RotPatch().getCylindricalCoordinateSystem().direction(),
                       scalar(0),
                       pp.name(),
                       true
                    );    
     
        
        
}

//- Construct as copy, resetting the face list and boundary mesh data
Foam::domainScalingPolyPatch::domainScalingPolyPatch
(
    const domainScalingPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    shadow_(pp.shadow_),    
    master_(pp.master_),
    origin_(pp.origin_),
    axis_(pp.axis_),
    direction_(pp.direction_),
//    relative_(false),    
    shadowIndex_(-1),       
    transformedFaceCellCentresPtr_(NULL),    
    masterTransformedPrimitivePointsPtr_(NULL),
    reconFaceCellCentresPtr_(NULL)
   
    
{
        
        rotPatch_= new primitiveDomainScalingPatch(
                       (*this), 
                       (*this).points(),
                       (*this),
                       pp.RotPatch().getCylindricalCoordinateSystem().origin(),
                       pp.RotPatch().getCylindricalCoordinateSystem().axis(),
                       pp.RotPatch().getCylindricalCoordinateSystem().direction(),
                       scalar(0),
                       pp.name(),
                       true
                    );           
        
        Info <<  "Construct as copy, resetting the face list and boundary mesh data" << endl;
//   faceList flist(0);
//   pointField pf(0);
//   List<point> pointCorners(0);
//   rotPatch_(flist,pf,*this, pointCorners);
//        
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::domainScalingPolyPatch::~domainScalingPolyPatch()
{
    clearOut();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
const Foam::domainScalingPolyPatch& Foam::domainScalingPolyPatch::shadow() const
{
    return refCast<const domainScalingPolyPatch>(boundaryMesh()[shadowIndex()]);
}

Foam::word Foam::domainScalingPolyPatch::shadowName() const
{
    return ( shadow_ );
}

Foam::label Foam::domainScalingPolyPatch::shadowIndex() const
{
    if (shadowIndex_ == -1)
    {
        polyPatchID shadow(shadow_, boundaryMesh());
        
        if (!shadow.active())
        {
            FatalErrorIn("label domainScalingPolyPatch::shadowIndex() const")
                << "Shadow patch name " << shadow_
                << " not found.  Please check your domainScaling interface definition."
                << abort(FatalError);
        }
        shadowIndex_ = shadow.index();

        if (!isType<domainScalingPolyPatch>(boundaryMesh()[shadowIndex_]))
        {
            FatalErrorIn("label domainScalingPolyPatch::shadow_Index() const")
                << "Shadow of domainScaling patch " << name()
                << " named " << shadowIndex() << " is not a domainScaling." << nl
                << "This is not allowed.  Please check your mesh definition."
                << abort(FatalError);
        }
    }

    return shadowIndex_;
}

const Foam::domainScalingPolyPatch& Foam::domainScalingPolyPatch::shadowPatch() const
{
    return refCast<const domainScalingPolyPatch>(boundaryMesh()[shadowIndex()]);
}

void Foam::domainScalingPolyPatch::calcReconFaceCellCentres() const
{
    Info<< "calcReconFaceCellCentres " << master()<<endl;
//    if (master())
//    {
        // Reconstruct the shadow cell face centres
    
    
    scalar deltaT =  this->boundaryMesh().time().timeOutputValue();
    
    const_cast<primitiveDomainScalingPatch&>(shadowPatch().RotPatch()).setDeltaTime(deltaT, this->boundaryMesh().time());
    
    
    
        primitiveDomainScalingPatch& rotPatch = const_cast<primitiveDomainScalingPatch&>(shadowPatch().RotPatch());
        rotPatch.deleteInterpolationPointers();
        
        const label shadowID = shadowIndex();
        
        vectorField  beforeCellCentres = boundaryMesh()[shadowID].faceCellCentres();
        
        vectorField reconCtrs =        
        shadowPatch().interpolateToShadow
        (
              beforeCellCentres                
        );
        vectorField localCtrs = this->faceCellCentres();

        // Calculate reconstructed centres by eliminating non-orthogonality        
        const vectorField& n = this->faceNormals();
        // localCtrs = localFaceCentres 
        // n
        reconFaceCellCentresPtr_ =                     
            new vectorField((localCtrs + n*(n & (reconCtrs -localCtrs))) );        
//    }
//    else
//    {
//        FatalErrorIn("void domainScalingPolyPatch::calcReconFaceCellCentres() const")
//            << "Attempting to create reconFaceCellCentres on a shadow"
//            << abort(FatalError);
//    }
}

const Foam::vectorField& Foam::domainScalingPolyPatch::reconFaceCellCentres() const
{
    Info<< "reconFaceCellCentres " <<master()<<endl;

    if (!reconFaceCellCentresPtr_)
    {
        calcReconFaceCellCentres();
    }
    
    return *reconFaceCellCentresPtr_;
}

void Foam::domainScalingPolyPatch::initGeometry()
{
    polyPatch::initGeometry();
}

void Foam::domainScalingPolyPatch::calcGeometry()
{
    // Reconstruct the cell face centres
    if (master())
    {
//        reconFaceCellCentres();
    }
}

void Foam::domainScalingPolyPatch::initMovePoints(const pointField& p)
{
    polyPatch::initMovePoints(p);
}


void Foam::domainScalingPolyPatch::movePoints(const pointField& p)
{
    polyPatch::movePoints(p);
    // TODO Here the points of the reconstructed face cell centres must be 
    // rotated if the moving is a rotation
    // if the moving is not a rotation a error message must be thrown.
    // or a translation must be handled
    // it should also be able to handle a combination translation and rotation


}

// TODO When the update should be performed the auxiliary patch should be newly calculated
void Foam::domainScalingPolyPatch::initUpdateMesh()
{
    polyPatch::initUpdateMesh();
}


void Foam::domainScalingPolyPatch::updateMesh()
{
    polyPatch::updateMesh();
    clearOut();
}

//- Initialize ordering (on new mesh)
void Foam::domainScalingPolyPatch::initOrder(const primitivePatch& pp) const
{
    
Info <<  "In init Order " << endl;
}

//- Return new ordering. Ordering is -faceMap: for every face index
//  the new face -rotation:for every new face the clockwise shift
//  of the original face. Return false if nothing changes (faceMap
//  is identity, rotation is 0)
bool Foam::domainScalingPolyPatch::order
(
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    faceMap.setSize(pp.size());
    faceMap = -1;

    rotation.setSize(pp.size());
    rotation = 0;

    return false;
}


// Write
void Foam::domainScalingPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);

     os // << "    nFaces " << size() << token::END_STATEMENT << nl
        //<< "    startFace " << start() << token::END_STATEMENT << nl
        << "    shadowPatch " << shadow_  << token::END_STATEMENT << nl              
        << "    hierarchy " << master_    << token::END_STATEMENT << nl
        << "    origin "    << origin_    << token::END_STATEMENT << nl
        << "    axis "      << axis_      << token::END_STATEMENT << nl
        << "    direction " << direction_ << token::END_STATEMENT << nl
        << "    omega "     << omega_     << token::END_STATEMENT << nl
        << endl;
}

/*
void Foam::domainScalingPolyPatch::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK << nl
        << "    type " << type() << token::END_STATEMENT << nl;
    patchIdentifier::writeDict(os);
    os  << "    nFaces " << size() << token::END_STATEMENT << nl
        << "    startFace " << start() << token::END_STATEMENT << nl
        << "    shadowPatch " << shadow_  << token::END_STATEMENT << nl              
        << "    hierarchy " << master_    << token::END_STATEMENT << nl
        << "    origin "    << origin_    << token::END_STATEMENT << nl
        << "    axis "      << axis_      << token::END_STATEMENT << nl
        << "    direction " << direction_ << token::END_STATEMENT << nl
        << "    omega "     << omega_     << token::END_STATEMENT << nl
        << token::END_BLOCK << endl;
}
*/

void  Foam::domainScalingPolyPatch::makeRotPatch() const
{
    rotPatch_ = new primitiveDomainScalingPatch
                    (
                            (*this),
                            (*this).points(),
                            (*this),
                            origin_,
                            axis_,
                            direction_,
                            omega_,
                            name(),
                            false
                    );
    (*rotPatch_).setStartAngle(startAngle_);
    
}
// ************************************************************************* //
