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

#include "mixingPlanePolyPatch.H"
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
    defineTypeNameAndDebug(mixingPlanePolyPatch, 0);
    
    addToRunTimeSelectionTable(polyPatch, mixingPlanePolyPatch, word);
    // addToRunTimeSelectionTable(polyPatch, mixingPlanePolyPatch, Istream);
    addToRunTimeSelectionTable(polyPatch, mixingPlanePolyPatch, dictionary);
    
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::mixingPlanePolyPatch::mixingPlanePolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
  /*  const word& masterName,
    const word& shadowName,    
    const vector& origin,
    const vector& axis,
    const vector& direction*/
)
:
    polyPatch(name, size, start, index, bm),
    shadow_("slave"),
    master_("master"),
    nameOmega_("Omega"),
    shadowIndex_(-1),
    rotPatch_(NULL)
{
        Info << "In constructor" << endl ;    
}

#if 0
// Construct from Istream
/*Foam::mixingPlanePolyPatch::mixingPlanePolyPatch
(
    Istream& is,
    const label index,
    const polyBoundaryMesh& bm
)
:

polyPatch(is, index, bm),
shadow_(is),
master_(is),
shadowIndex_(-1),
rotPatch_(NULL) 
{   
        
    /*const point& origin(is); // Origin
    const vector& axis(is); // Axis
    const vector& direction(is); // Direction
*/    
/* ;        
        
   } */  
#endif
    
// Construct from dictionary
Foam::mixingPlanePolyPatch::mixingPlanePolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, dict, index, bm),
    shadow_(dict.lookup("shadowPatch")),
    master_(dict.lookup("hierarchy")),  
    nameOmega_("Omega"), 
    shadowIndex_(-1),
    origin_(dict.lookup("origin")),
    axis_(dict.lookup("axis")),
    direction_(dict.lookup("direction")),   
    rotPatch_(NULL)    
{
        Info << "In constructor" << endl;
}

//- Construct as copy, resetting the boundary mesh
Foam::mixingPlanePolyPatch::mixingPlanePolyPatch
(
    const mixingPlanePolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    shadow_(pp.shadow_),
    master_(pp.master_),    
    nameOmega_(pp.nameOmega_),
    shadowIndex_(-1)
{
//        rotPatch_.writePatches(pp.name())     
        Info << "In constructor" << endl;
        origin_ = pp.origin();
        axis_ = pp.axis();
        direction_ = pp.direction();
             
}

//- Construct as copy, resetting the face list and boundary mesh data
Foam::mixingPlanePolyPatch::mixingPlanePolyPatch
(
    const mixingPlanePolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart),
    shadow_(pp.shadow_),
    master_(pp.master_),    
    nameOmega_(pp.master_),
    shadowIndex_(-1)
{  
        Info << "In constructor" << endl;
        origin_ = pp.origin();
        axis_ = pp.axis();
        direction_ = pp.direction();
//        rotPatch_ = new primitiveMixingPlanePatch
//                           (
//                                   (*this),
//                                   (*this).points(),
//                                   (*this),
//                                   pp.RotPatch().getCylindricalCoordinateSystem().origin(),
//                                   pp.RotPatch().getCylindricalCoordinateSystem().axis(),
//                                   pp.RotPatch().getCylindricalCoordinateSystem().direction(),                             
//                                   pp.name(),
//                                   false
//                           );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mixingPlanePolyPatch::~mixingPlanePolyPatch()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::mixingPlanePolyPatch& Foam::mixingPlanePolyPatch::shadow() const
{
    return refCast<const mixingPlanePolyPatch>(boundaryMesh()[shadowIndex()]);
}

Foam::word Foam::mixingPlanePolyPatch::shadowName() const
{
    return ( shadow_ );
}


Foam::label Foam::mixingPlanePolyPatch::shadowIndex() const
{
    if (shadowIndex_ == -1)
    {
        // Grab shadow patch index
        polyPatchID shadow(shadow_, boundaryMesh());
        
        // TODO the does not have class type error must be resolved
        if (!shadow.active())
        {
            FatalErrorIn("label mixingPlanePolyPatch::shadowIndex() const")
                << "Shadow patch name " << shadow_
                << " not found.  Please check your mixingPlane interface definition."
                << abort(FatalError);
        }

        shadowIndex_ = shadow.index();
        
        if (!isType<mixingPlanePolyPatch>(boundaryMesh()[shadowIndex_]))
        {
            FatalErrorIn("label mixingPlanePolyPatch::shadow_Index() const")
                << "Shadow of mixingPlane patch " << name()
                << " named " << shadowIndex() << " is not a mixingPlane." << nl
                << "This is not allowed.  Please check yuor mesh definition."
                << abort(FatalError);
        }
    }

    return shadowIndex_;
}

const Foam::mixingPlanePolyPatch& Foam::mixingPlanePolyPatch::shadowPatch() const
{
    if(shadowIndex_ == -1)
        shadowIndex();
    
    return refCast<const mixingPlanePolyPatch>(boundaryMesh()[shadowIndex()]);
}


void Foam::mixingPlanePolyPatch::adaptToMaster() const
{
    if( !master() )
    {
        if(rotPatch_)
        {
            (*rotPatch_).rebuildPatchForEdges
            (
                    shadowPatch().getPointsSide()
            );
        }
        else
        {
            makeRotPatch();
            (*rotPatch_).rebuildPatchForEdges
            (
                    shadowPatch().getPointsSide()
            );           
            
        }
    }   
}


// Write
void Foam::mixingPlanePolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);

    //os  << nl << shadow_ << endl;
    os  //<< "    nFaces " << size() << token::END_STATEMENT << nl
       // << "    startFace " << start() << token::END_STATEMENT << nl
        << "    shadowPatch " << shadow_ << token::END_STATEMENT << nl              
        << "    hierarchy " << master_ << token::END_STATEMENT << nl
        << "    origin " << (*this).origin() << token::END_STATEMENT << nl
        << "    axis " << (*this).axis() << token::END_STATEMENT << nl
        << "    direction " << (*this).direction() << token::END_STATEMENT << nl
        << endl;
}


/*void Foam::mixingPlanePolyPatch::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK << nl
        << "    type " << type() << token::END_STATEMENT << nl;
    patchIdentifier::writeDict(os);
    
    os  << "    nFaces " << size() << token::END_STATEMENT << nl
        << "    startFace " << start() << token::END_STATEMENT << nl
        << "    shadowPatch " << shadow_ << token::END_STATEMENT << nl              
        << "    hierarchy " << master_ << token::END_STATEMENT << nl
        << "    origin " << (*this).origin() << token::END_STATEMENT << nl
        << "    axis " << (*this).axis() << token::END_STATEMENT << nl
        << "    direction " << (*this).direction() << token::END_STATEMENT << nl
        << token::END_BLOCK << endl;
	}*/

void  Foam::mixingPlanePolyPatch::makeRotPatch() const
{
    rotPatch_ = new primitiveMixingPlanePatch
                    (
                            (*this),
                            (*this).points(),
                            (*this),
                            origin_,
                            axis_,
                            direction_,                           
                            name(),
                            false
                    );
    
}
// ************************************************************************* //

