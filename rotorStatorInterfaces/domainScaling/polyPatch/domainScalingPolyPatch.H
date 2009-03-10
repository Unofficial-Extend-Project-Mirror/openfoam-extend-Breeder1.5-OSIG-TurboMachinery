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
    domainScalingPolyPatch

Description
    domainScaling Patch

Authors
    Franz Blaim

SourceFiles
    domainScalingPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef domainScalingPolyPatch_H
#define domainScalingPolyPatch_H

#include "coupledPolyPatch.H"
#include "primitiveDomainScalingPatch.H"
#include "ggiInterpolation.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class domainScalingPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class domainScalingPolyPatch
:
    public coupledPolyPatch
{
    // Private data

        /**
         * \brief Shadow patch name
         */
        const word shadow_;
        
        /**
         * \brief master name
         */ 
        const word master_;

        /**
         * \brief projectionVectorField
         */
        mutable vectorField slaveToMasterProjectVector_;
                
        mutable vectorField masterToSlaveProjectVector_;
        
        /**
         * \brief Origin of the cylindrical coordinate system
         */
        const vector origin_;
        
        /**
         * \brief Axis of the cylindrical coordinate system
         */
        const vector axis_;
        
        /**
         * \brief Direction of the reference axis for measuring the angle
         */
        const vector direction_;
     
        /**
         * \brief Shadow patch index.  Delayed evaluation for construction
         */
        mutable label shadowIndex_;

        /**
         * \brief Reconstruct masterFaceCellCentres pointField to use for GGinterpolation
         */
        mutable pointField* transformedFaceCellCentresPtr_;
        
        mutable pointField* masterTransformedPrimitivePointsPtr_;             
        
        mutable primitiveDomainScalingPatch* rotPatch_;
        
        /**
         * \brief Reconstructed patch neighbour cell centres
         */
        mutable vectorField* reconFaceCellCentresPtr_;
        
        /**
         * \brief Reconstruct master primitivePatch to use for GGinterpolation
         */
        mutable primitivePatch* transformedPrimitivePatchPtr_;        
         
        /**
         * \brief Stores the angular frequency of for the primitive patch.
         */
        scalar omega_;
        
        /**
         * \brief Stores the start face label for the primitive patch. 
         */
        label start_;
        
        /**
         * \brief Stores the number of faces.
         */
        label nFaces_;
        
        /**
         * \brief The start angle for the interface. Can be used
         * for the frozen rotor approach.
         */
        scalar startAngle_;
        


protected:

    // Protected Member functions

        /**
         * \brief Initialise the calculation of the patch geometry
         */
        virtual void initGeometry();

        /**
         * \brief Calculate the patch geometry
         */
        virtual void calcGeometry(); // Todo here should the geometry of the 
        // auxpatch calculated

        /**
         * \brief Initialise the patches for moving points
         */
        virtual void initMovePoints(const pointField&);

        /**
         * \brief Correct patches after moving points
         */
        virtual void movePoints(const pointField&);

        /**
         * \brief Initialise the update of the patch topology
         */
        virtual void initUpdateMesh();

        /**
         * \brief Update of the patch topology
         */
        virtual void updateMesh();


        
        /**
         * \brief Calculate reconstructed cell centres
         */
        void calcReconFaceCellCentres() const;

        /**
         * \brief Clear out
         */
        void clearOut();

public:

    //- Runtime type information
    TypeName("domainScaling");


    // Constructors

        /**
         * \brief Construct from components
         */
        domainScalingPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm
            /*const word& masterName,
            const word& slaveName,            
            const vector& origin,
            const vector& axis,
            const vector& direction*/
            
        );

        /**
         * \brief Construct from Istream
         */
        domainScalingPolyPatch
        (
            Istream& is,
            const label index,
            const polyBoundaryMesh&
        );

        /**
         * \brief Construct from dictionary
         */
        domainScalingPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh&
        );

        /**
         * \brief Construct as copy, resetting the boundary mesh
         */
        domainScalingPolyPatch
        (
            const domainScalingPolyPatch&,
            const polyBoundaryMesh&
        );

        /**
         * \brief Construct as given the original patch and resetting the
         * face list and boundary mesh information
         */
        domainScalingPolyPatch
        (
            const domainScalingPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        /**
         * \brief Construct and return a clone, resetting the boundary mesh
         */
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new domainScalingPolyPatch(*this, bm));
        }

        /**
         * \brief Construct and return a clone, resetting the face list
         *  and boundary mesh
         */
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
                new domainScalingPolyPatch
                (
                    refCast<const domainScalingPolyPatch>(*this),
                    bm,
                    index,
                    newSize,
                    newStart
                )
            );
        }


    // Destructor

        virtual ~domainScalingPolyPatch();


    // Member functions

        /**
         * \brief Return shadow patch name
         */
        const domainScalingPolyPatch& shadow() const;

        /**
         * \brief Return shadow patch index
         */
        int shadowIndex() const;

        /**
         * \brief Is this the master side?
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
               FatalErrorIn("Valid entry for word hierarchy in domainScalingPolyPatch are master or slave")
                << "Change " << master_ << " to a valid entry"
                << " in patch " << name() << " definition in constant/polyMesh/boundary "
                << abort(FatalError);
                return -1;
            }
        }

        /**
         * \brief Is this the slave side?
         */
        bool slave() const
        {
            return !master();
        }

        /**
         * \brief Return shadow patch
         */
        const domainScalingPolyPatch& shadowPatch() const;
        
        word shadowName() const;
       

        /**
         * \brief Interpolate face field on the shadow side.         
         */
        template<class Type>
        tmp<Field<Type> > interpolateToShadow(const Field<Type>& pf) const
        {
            if(rotPatch_)
            {
                return (*rotPatch_).interpolateToNeighbour(pf, shadowPatch());
            }
            else
            {
                makeRotPatch();  
                return (*rotPatch_).interpolateToNeighbour(pf, shadowPatch());
            }
        }
        
        /**
         * \brief Interpolate face field on the shadow side. 
         */        
        template<class Type>
        tmp<Field<Type> > interpolateToShadow
        (
                const Field<Type>& pf,
                const direction cmpt,
                const int& rank
        ) const
        {
            if(rotPatch_)
            {
                return (*rotPatch_).interpolateToNeighbour
                            (
                                    pf, shadowPatch(),
                                    cmpt,
                                    rank
                            );
            }
            else
            {
                makeRotPatch();  
                return (*rotPatch_).interpolateToNeighbour
                                          (
                                                  pf, shadowPatch(),
                                                  cmpt,
                                                  rank
                                          );
                
            }
        }
        
        /**
         * \brief Interpolate given field on the shadow side.         
         */
        template<class Type> tmp<Field<Type> > interpolateToShadow(
                const tmp<Field<Type> >& tpf) const
        {    

            if(rotPatch_)
            {
                return (*rotPatch_).interpolateToNeighbour(tpf, shadowPatch());
            }
            else
            {
                makeRotPatch(); 
                return (*rotPatch_).interpolateToNeighbour(tpf, shadowPatch());
            }
            
        }

        /**
         * \brief Return reconstructed cell centres
         */
        const vectorField& reconFaceCellCentres() const;

         
        /**
         * \brief Initialize ordering for primitivePatch. Does not
         *  refer to *this (except for name() and type() etc.)
         */
        virtual void initOrder(const primitivePatch&) const;

        /**
         *  Return new ordering for the primitivePatch.
         *  Ordering is -faceMap: for every face
         *  index of the new face -rotation: for every new face the clockwise
         *  shift of the original face. Return false if nothing changes
         *  (faceMap is identity, rotation is 0), true otherwise.
         */
        virtual bool order
        (
            const primitivePatch&,
            labelList& faceMap,
            labelList& rotation
        ) const;

        
        /**
         * \brief Sets the actual time and the actual time step delta.         
         */
        void setDeltaTime(scalar deltaT,const Time& runTime)
        {
            (*rotPatch_).setDeltaTime(deltaT,runTime);                  
        }
        

        /**
         * \brief Write all the parameters to the boundary file.
         */
        virtual void write(Ostream&) const;

        /**
         * \brief Write dictionary
         */
        //virtual void writeDict(Ostream&) const;
        
        /**
         * \brief Can be used to force the interface patch to rotate 
         * about a specified angle.
         */
        void adoptToOuterTimeStep(scalar angle)
        {
            (*rotPatch_).adoptToOuterTimeStep(angle);            
        }        
        
        /**
         * \brief Returns the reference to the primitiveRotationalPatch. 
         * 
         */
        const primitiveDomainScalingPatch& RotPatch() const
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
        
        /**
         * \brief Generates the primitiveRotationalPatch pointer.
         */        
        void makeRotPatch() const;
     
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
