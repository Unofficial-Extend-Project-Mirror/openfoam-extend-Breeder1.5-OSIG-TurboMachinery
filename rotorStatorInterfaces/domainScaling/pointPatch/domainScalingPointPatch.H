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
    domainScalingPointPatch

Description
    Dummy domainScalingPointPatch for post-processing.  No functionality built in!

SourceFiles
    domainScalingPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef domainScalingPointPatch_H
#define domainScalingPointPatch_H

#include "facePointPatch.H"
#include "faceList.H"
#include "domainScalingPolyPatch.H"
#include "coupledFacePointPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class domainScalingPointPatch Declaration
\*---------------------------------------------------------------------------*/

class domainScalingPointPatch
:
    public coupledFacePointPatch
{
    // Private Member Functions

        /**
         * \brief Disallow default construct as copy
         */
        domainScalingPointPatch(const domainScalingPointPatch&);

        /**
         * \brief Disallow default assignment
         */
        void operator=(const domainScalingPointPatch&);

        virtual void initGeometry()
        {   }
        //- Calculate the patch geometry
        virtual void calcGeometry()
        {}

        //- Initialise the patches for moving points
        virtual void initMovePoints(const pointField&)
	{}

        //- Correct patches after moving points
        virtual void movePoints(const pointField&)
        {}

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh()
        {
	    facePointPatch::initUpdateMesh();
        }

        //- Update of the patch topology
        virtual void updateMesh()
        {
	    facePointPatch::updateMesh();	    
	}


public:

    typedef pointBoundaryMesh BoundaryMesh;


    //- Runtime type information
    TypeName(domainScalingPolyPatch::typeName_());


    // Constructors

        /**
         * \brief Construct from components
         */
        domainScalingPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        );


    // Destructor

        virtual ~domainScalingPointPatch();


    // Member Functions

        // Access

            /**
             * \brief Return true because this patch is coupled
             */
            virtual bool coupled() const
            {
                return true;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
