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

Description
    Generalized grid interface (GGI) patch, providing coupling
    between arbitrary patches which belong to the same fvMesh

Authors
	Franz Blaim

\*---------------------------------------------------------------------------*/

#include "domainScalingFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvBoundaryMesh.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(domainScalingFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, domainScalingFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// Make patch weighting factors
void Foam::domainScalingFvPatch::makeWeights(scalarField& w) const
{   
    
    vectorField n = nf(); //
        scalarField nMe = n & fvPatch::delta();
        // & scalar product
        scalarField nfc =  shadow().interpolateToShadow( shadow().nf() & shadow().delta()); //n & (domainScalingPolyPatch_.reconFaceCellCentres() - Cf());

        w = nfc /( nMe  + nfc);
}

Foam::tmp<Foam::labelField> Foam::domainScalingFvPatch::transfer
(
    const Pstream::commsTypes,
    const unallocLabelList& interfaceData
) const
{
    notImplemented
    (
        "cyclicGgiFvPatchField<Type>::"
        "transfer(const unallocLabelList& interfaceData) const"
    );

    return labelField::null();
}

Foam::tmp<Foam::labelField> Foam::domainScalingFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes,
    const unallocLabelList& iF
) const
{
    return shadow().patchInternalField(iF);
}


// Make patch face - neighbour cell distances
void Foam::domainScalingFvPatch::makeDeltaCoeffs(scalarField& dc) const
{
    scalarField nfc =  shadow().interpolateToShadow( shadow().nf() & shadow().delta()); //n & (domainScalingPolyPatch_.reconFaceCellCentres() - Cf());
      dc = 1.0/(nfc + (nf() & delta()) );
 
  //  Info << "deltas" << ( nf() & delta() ) << endl;
//    Info << "Make Delta coeffs " << dc << endl;
}


Foam::tmp<Foam::vectorField> Foam::domainScalingFvPatch::delta() const
{   
    return Cf() - Cn();
}


const Foam::domainScalingFvPatch& Foam::domainScalingFvPatch::shadow() const
{
    const fvPatch& p = this->boundaryMesh()[domainScalingPolyPatch_.shadowIndex()];

    return refCast<const domainScalingFvPatch>(p);
}

Foam::tmp<Foam::labelField> Foam::domainScalingFvPatch::interfaceInternalField
(
    const unallocLabelList& internalData
) const
{

   notImplemented("interfaceinternalField");
 
   return patchInternalField(internalData);
}


Foam::tmp<Foam::labelField> Foam::domainScalingFvPatch::transfer
(
    const unallocLabelList& interfaceData
) const
{
    notImplemented
    (
        "domainScalingFvPatchField<Type>::"
        "transfer(const unallocLabelList& interfaceData) const"
    );

    return labelField::null();
}


Foam::tmp<Foam::labelField> Foam::domainScalingFvPatch::internalFieldTransfer
(
    const unallocLabelList& iF
) const
{
    notImplemented( "internalFieldTransfer");
    return shadow().patchInternalField(iF);
}

//void Foam::domainScalingFvPatch::makeCorrVecs(vectorField& cv) const
//{
//
//    cv = vector::zero;
//}

// ************************************************************************* //
