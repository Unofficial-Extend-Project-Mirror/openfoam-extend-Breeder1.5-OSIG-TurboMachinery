/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

Application
    initializeRotation

Description
    Initializes rotational flow field for pump calculations

Author
    Mikko Auvinen, Aalto University

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  const vectorField& cellCenter = mesh.C();
  vector null = vector::zero;

  forAll( cellCenter , cellId )
  {
    vector cc = ( cellCenter[cellId] - origin );
    scalar r  = mag( rotationAxis ^ cc ); // radius
    if( r <= refRadius )
    {
      U[cellId] = omega ^ cc;
    }
    else
    {
    // Usually, it's best to kill the rotation quickly.
      U[cellId] = (omega ^ cc) * ::pow( (refRadius/r), 12 );
    }
  }

    U.write();
    return(0);
}


// ************************************************************************* //
