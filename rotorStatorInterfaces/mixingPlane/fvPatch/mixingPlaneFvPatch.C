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

#include "mixingPlaneFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "fvBoundaryMesh.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mixingPlaneFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, mixingPlaneFvPatch, polyPatch);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //



// Make patch face - neighbour cell distances
void Foam::mixingPlaneFvPatch::makeDeltaCoeffs(scalarField& dc) const
{   
    dc = 1.0/(nf() & delta());
}

Foam::tmp<Foam::vectorField> Foam::mixingPlaneFvPatch::delta() const
{
     return fvPatch::delta(); // muesste nich eigentlich immer das selbe verwendet werden ?
}


const Foam::mixingPlaneFvPatch& Foam::mixingPlaneFvPatch::shadow() const
{
    const fvPatch& p = this->boundaryMesh()[mixingPlanePolyPatch_.shadowIndex()];

    return refCast<const mixingPlaneFvPatch>(p);
}

Foam::tmp<Foam::labelField> Foam::mixingPlaneFvPatch::interfaceInternalField
(
    const unallocLabelList& internalData
) const
{

   notImplemented("interfaceinternalField");
 
   return patchInternalField(internalData);
}


Foam::tmp<Foam::labelField> Foam::mixingPlaneFvPatch::transfer
(
    const unallocLabelList& interfaceData
) const
{
    notImplemented
    (
        "mixingPlaneFvPatchField<Type>::"
        "transfer(const unallocLabelList& interfaceData) const"
    );

    return labelField::null();
}


Foam::tmp<Foam::labelField> Foam::mixingPlaneFvPatch::internalFieldTransfer
(
    const unallocLabelList& iF
) const
{
    notImplemented( "internalFieldTransfer");
    return shadow().patchInternalField(iF);
}

Foam::tmp<Foam::labelField> Foam::mixingPlaneFvPatch::transfer
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


Foam::tmp<Foam::labelField> Foam::mixingPlaneFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes,
    const unallocLabelList& iF
) const
{
    return shadow().patchInternalField(iF);
}


// ************************************************************************* //
