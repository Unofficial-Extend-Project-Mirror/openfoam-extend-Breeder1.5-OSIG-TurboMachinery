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
    mixingPlanePolyPatch

Description
    mixingPlane Patch

Authors
    Franz Blaim

SourceFiles
    mixingPlanePolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef mixingPlanePolyPatch_H
#define mixingPlanePolyPatch_H

#include "polyPatch.H"
#include "primitiveMixingPlanePatch.H"
#include "ggiInterpolation.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class mixingPlanePolyPatch Declaration
\*---------------------------------------------------------------------------*/

class mixingPlanePolyPatch
:
    public polyPatch
{
    // Private data

        //- Shadow patch name        
        const word shadow_;
        
        
        // master name 
        const word master_;
          
        //- Shadow patch index.  Delayed evaluation for construction
        mutable label shadowIndex_;
        
        point origin_;
        
        vector axis_;
        
        vector direction_;
        
        mutable primitiveMixingPlanePatch* rotPatch_;
        
protected:

        // Protected Member functions

            //- Initialise the calculation of the patch geometry
            virtual void initGeometry()
            {
                polyPatch::initGeometry();
            }

            //- Calculate the patch geometry
            virtual void calcGeometry()
            {
            }
            
            //- Initialise the patches for moving points
            virtual void initMovePoints(const pointField& p)
            {
                polyPatch::initMovePoints(p);
            }

            //- Correct patches after moving points
            virtual void movePoints(const pointField& p)
            {
                polyPatch::movePoints(p);
            }

            //- Initialise the update of the patch topology
            virtual void initUpdateMesh()
            {
                polyPatch::initUpdateMesh();
            }


            //- Update of the patch topology
            virtual void updateMesh()
            {
                polyPatch::updateMesh();
            }


public:

    //- Runtime type information
    TypeName("mixingPlane");


    // Constructors

        //- Construct from components
        mixingPlanePolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm
         /*   const word& masterName,
            const word& slaveName,            
            const vector& origin,
            const vector& axis,
            const vector& direction*/
            
        );
        
        //- Construct from istream
        mixingPlanePolyPatch
        (
            Istream& is,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Construct from dictionary
        mixingPlanePolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh&
        );

        //- Construct as copy, resetting the boundary mesh
        mixingPlanePolyPatch
        (
            const mixingPlanePolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct as given the original patch and resetting the
        //  face list and boundary mesh information
        mixingPlanePolyPatch
        (
            const mixingPlanePolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new mixingPlanePolyPatch(*this, bm));
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new mixingPlanePolyPatch
                (
                    refCast<const mixingPlanePolyPatch>(*this),
                    bm,
                    index,
                    newSize,
                    newStart
                )
            );
        }


    // Destructor

        virtual ~mixingPlanePolyPatch();


    // Member functions

        //- Return shadow patch name
        const mixingPlanePolyPatch& shadow() const;

        //- Return shadow patch index
        int shadowIndex() const;

        /**
         * \brief{Is this the master ?.}
         * \date{21-05-2008}
         * The function returns true if this is the master side.
         */
        bool master() const
        {
            if(master_ == "master")
            {
               return true;
            }
            else if(master_ == "slave")
            {
              return false;
            }
            else
            {
               FatalErrorIn("Valid entry for word hierarchy in mixingPlanePolyPatch are master or slave")
                << "Change " << master_ << " to a valid entry"
                << " in patch " << name() << " definition in constant/polyMesh/boundary "
                << abort(FatalError);
                return -1;
            }
        }
        
        /**
         * \brief{Is this the slave ?.}
         * \date{21-05-2008}
         * The function returns true if this is the slave side.
         */
        bool slave() const
        {
            return !master();
        }

        /**
        * \brief{Give back the shadowPatch.}
        * \date{21-05-2008}
        * The function returns true if this is the slave side.
        */
        const mixingPlanePolyPatch& shadowPatch() const;
        
        
        const List<point>& getPointsSide() const
        {   
            if(rotPatch_)
            {
                return ( (*rotPatch_).getPointsSide() );
            }
            else
            {
                makeRotPatch();
                return ( (*rotPatch_).getPointsSide() );
            }
                
        }
        
        
        /**
        * \brief{Is this the slave ?.}
        * \date{21-05-2008}
        * The function returns true if this is the slave side.
        */
        word shadowName() const;
        
        
        void adaptToMaster() const;
      

        //- Write
        virtual void write(Ostream&) const;

        //- Write dictionary
        //virtual void writeDict(Ostream&) const;
        
        
        /**
        * \brief{Is this the slave ?.}
        * \date{01-04-2008}
        * The function returns true if this is the slave side.
        */
        const primitiveMixingPlanePatch& RotPatch() const
        {
            if(rotPatch_)
            {
                return ( *rotPatch_);
            }
            else
            {
                makeRotPatch();
                return ( *rotPatch_);
            }
        }
        
        void makeRotPatch() const;
        
        vector axis()const
        {
            return axis_;
        }
        
        point origin()const
        {
            return origin_;
        }
        
        vector direction()const
        {
            return direction_;
        }
        
//        void writePatches():
            
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
