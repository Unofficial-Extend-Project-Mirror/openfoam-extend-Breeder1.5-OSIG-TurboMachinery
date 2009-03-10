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
    domainScalingFvPatch

Description
    The finiteVolumePatch for the domainScaling interface. 
    This handles all the geometric related factor necessary within the 
    finite Volume approach.

Authors
    Franz Blaim

SourceFiles
    ixingPlaneFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef domainScalingFvPatch_H
#define domainScalingFvPatch_H

#include "coupledFvPatch.H"
#include "domainScalingLduInterface.H"
#include "domainScalingPolyPatch.H"
#include "primitiveDomainScalingPatch.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class domainScalingFvPatch Declaration
\*---------------------------------------------------------------------------*/



class domainScalingFvPatch
:
    public coupledFvPatch,
    public domainScalingLduInterface
{
    // Private Data
     
private: 
        /**
         * \brief Reference to the corresponding polyPatch.
         */        
        const domainScalingPolyPatch& domainScalingPolyPatch_;

protected:

        // Protected Member functions
        /**
         * \brief Weights generation for the evaluation of the face values.
         */
        virtual void makeWeights(scalarField&) const;

        /**
         * \brief Make patch face - neighbour cell distances
         */
        virtual void makeDeltaCoeffs(scalarField&) const;

        /**
         * \brief Initialise the patches for moving points
         */
        virtual void initMovePoints()
        {}

        /**
         * \brief Correct patches after moving points
         */
        virtual void movePoints()
        {}



public:

    //- Runtime type information
    TypeName(domainScalingPolyPatch::typeName_());
    
    // Constructors

        /**
         * \brief Construct from components
         */
        domainScalingFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :            
            coupledFvPatch(patch, bm),
            domainScalingPolyPatch_(refCast<const domainScalingPolyPatch>(patch))           
        {
        }


    

        /**
         * \brief Destructor
         */
        virtual ~domainScalingFvPatch()
        {}


    // Member functions

        // Access
            /**
            * \brief Return slave patch index.
            */
            int slaveIndex() const
            {
                return domainScalingPolyPatch_.shadowIndex();
            }

            /**
             * \brief Returns wether this is a master or not.
             */
            bool master() const
            {
                return domainScalingPolyPatch_.master();
            }


            /**
             * \brief Returns wether this is the slave side.
             */
            bool slave() const
            {
                return !master();
            }
            
            
            /**
             * \brief This checks if the slave patch has really setted a master
             */
            bool checkDoubleMaster() const
            {
                if(master())
                {
                   return domainScalingPolyPatch_.shadowPatch().master();
                }
                else
                {
                   return false;
                }
            }
            
            /**
            * \brief Return slave patch.
            */
            const domainScalingFvPatch& shadow() const;

        // Interface transfer functions
            
            /**
            * \brief Return slave patch.
            */
            virtual tmp<labelField> interfaceInternalField
            (
                const unallocLabelList& internalData
            ) const;
       
            
            /**
            * \brief Transfer and return neighbour field.
            */
            virtual tmp<labelField> transfer
            (
                const unallocLabelList& interfaceData
            ) const;

            
            virtual tmp<vectorField> delta() const;
            
            /**
            * \brief Return neighbour field.
            */
            virtual tmp<labelField> internalFieldTransfer
            (
                const unallocLabelList& internalData
            ) const;
            
            //- Transfer and return neighbour field
            virtual tmp<labelField> transfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& interfaceData
            ) const;

            //- Return neighbour field
            virtual tmp<labelField> internalFieldTransfer
            (
                const Pstream::commsTypes commsType,
                const unallocLabelList& internalData
            ) const;


            /**
             * \brief Sets the delta time and the actual runtim within the 
             * domainScalingPolyPatch_ member variable.
             */
            void setDeltaTime(scalar deltaT,const Time& runTime)
            {
               const_cast<domainScalingPolyPatch&>( domainScalingPolyPatch_ ).setDeltaTime(deltaT,runTime);                  
            }
            
            /**
             * \brief Return the face to Cells adressing label list.
             */
            virtual const unallocLabelList& faceCells() const
            {
                return ( coupledFvPatch::faceCells() );
            }
            

            /**
             * \brief Interpolate from the original patch to the shadow patch, 
             * without knowing about the rank and the components.
             * \param tpf temporary field object, which should be interpolated to the shadow patch.
             */
            template<class Type> 
            tmp<Field<Type> > interpolateToShadow(const Field<Type>& pf) const
            {
                return domainScalingPolyPatch_.interpolateToShadow(pf);
            }

            /**
             * \brief  Interpolate from the original patch to the shadow patch, 
             * without knowing about the rank and the components.
             * \param tpf temporary field object, which should be interpolated to the shadow patch. 
             */
            template<class Type> 
            tmp<Field<Type> > interpolateToShadow(const tmp<Field<Type> >& tpf) const
            {
                 return domainScalingPolyPatch_.interpolateToShadow(tpf);
            }
            
            /**
             * \brief  Interpolate from the original patch to the shadow patch.
             * \param pf field reference, which should be interpolated to the shadow patch.
             * \param cmpt component of the actual tensor.
             * \param rank The rank of the actual tensor field. 
             */
            template<class Type> 
            tmp<Field<Type> > interpolateToShadow
            (
                    const Field<Type>& pf,
                    const direction& cmpt,
                    const int& rank
            ) const
            {
                return domainScalingPolyPatch_.interpolateToShadow(pf,cmpt,rank);
            }
            
            /**
             * \brief Interpolate from the original patch to the shadow patch.
             * \param tpf temporary field object which should be interpolated to the shadow patch.
             * \param cmpt component of the actual tensor.
             * \param rank The rank of the actual tensor field.
             */
            template<class Type> 
            tmp<Field<Type> > interpolateToShadow
            (
                    const tmp<Field<Type> >& tpf,
                    const direction& cmpt,
                    const int& rank
            ) const
            {
                 return domainScalingPolyPatch_.interpolateToShadow(tpf,cmpt,rank);
            }
            
            /**
             * \brief Returns constant reference on the domain scaling poly patch.
             */
            const domainScalingPolyPatch& GetdomainScalingPolyPatch() const
            {
                return (domainScalingPolyPatch_);
            }
            
            /**
             * \brief Returns a constant reference to the primitiveDomainScalingPatch.
             */
            const primitiveDomainScalingPatch& RotPatch() const
            {
                const primitiveDomainScalingPatch& rotPatch = domainScalingPolyPatch_.RotPatch();
                return ( rotPatch );                
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
