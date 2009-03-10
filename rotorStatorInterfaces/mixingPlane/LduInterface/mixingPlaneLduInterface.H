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
    mixingPlaneLduInterface

Description
    An abstract base class for mixingPlane coupled interfaces.

Authors
    Franz Blaim

SourceFiles
    mixingPlaneLduInterface.C

\*---------------------------------------------------------------------------*/

#ifndef mixingPlaneLduInterface_H
#define mixingPlaneLduInterface_H

#include "lduInterface.H"
#include "primitiveFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class mixingPlaneLduInterface Declaration
\*---------------------------------------------------------------------------*/

class mixingPlaneLduInterface
:
    public lduInterface
{

public:
    
    //- Runtime type information
    TypeName("mixingPlaneLduInterface");

    // Constructors

    // Destructor

        virtual ~mixingPlaneLduInterface();


    // Member Functions

        // Access

            //- Return face transformation tensor
           // virtual const tensor& forwardTGGI() const = 0;

            //- Return face reverse transformation tensor
            //virtual const tensor& backwardTGGI() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
