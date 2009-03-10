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
    domainScalingFvPatchField

Description
    The fvPatchField class for the domainScaling patch

Authors
    Franz Blaim

SourceFiles
    domainScalingFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef domainScalingFvPatchField_H
#define domainScalingFvPatchField_H

#include "coupledFvPatchField.H"
#include "domainScalingLduInterfaceField.H"
#include "domainScalingFvPatch.H"
#include "transform.H"
#include "scalarField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class domainScalingFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class domainScalingFvPatchField
:
    public domainScalingLduInterfaceField,
    public coupledFvPatchField<Type>
{
    // Private data

       /**
        * \brief The constant reference on the domainScalingFvPatch,
        * it is mainly used to retrieve the reference on the primitiveDomainScalingPatch.
        * This reference is used to call the interpolation functions.
        */
       const domainScalingFvPatch& domainScalingPatch_; 
     
       /**
        * \brief The name of the actual field, it is mainly used 
        * to retrieve the correct shadow patch field.
        */
       word fieldName_; 
        
       
       
      
public:

    /**
     * \brief Runtime type information
     */
    TypeName(domainScalingFvPatch::typeName_());


    // Constructors

        /**
         * \brief Construct from patch and internal field
         */
        domainScalingFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        /**
         * \brief Construct from patch, internal field and dictionary
         */        
        domainScalingFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        /**
         * \brief Construct by mapping given domainScalingFvPatchField onto a new patch
         */
        domainScalingFvPatchField
        (
            const domainScalingFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        /**
         * \brief Construct as copy setting internal field reference
         */
        domainScalingFvPatchField
        (
            const domainScalingFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        /**
         * \brief Construct and return a clone
         */
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new domainScalingFvPatchField<Type>(*this)
            );
        }

        /**
         * \brief Construct and return a clone setting internal field reference
         */
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new domainScalingFvPatchField<Type>(*this, iF)
            );
        }
        
        
        /**
         * \brief
         */
        const domainScalingFvPatchField<Type>& shadowPatchField() const;
        
        /**
         *  \brief Returns the rank of the corresponding tensor field.
         */
        virtual int rank() const
        {
           return pTraits<Type>::rank;
        }
        
      

    // Member functions

        // Evaluation functions

            /**
             * \brief Return neighbour field given internal cell data
             */
            virtual tmp<Field<Type> > patchNeighbourField() const;
           
            /**
             * \brief Evaluate the patch field
             */
            virtual void evaluate
            (
                const Pstream::commsTypes commsType
            );
            
          
            /**
             * \brief Initialize the update of the matrix.
             * In this function nothing important is done, but it 
             * has to be implemented.
             */
            virtual void initInterfaceMatrixUpdate
            (
                const scalarField& psiInternal,
                scalarField&,
                const lduMatrix&,
                const scalarField&,
                const direction,
                const bool bufferdTransfer
            ) const;
            
            /**
            * \brief Update result field based on interface functionality 
            */
/*            virtual void updateInterfaceMatrix
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix&,
                const scalarField& coeffs,
                const direction cmpt
            ) const;*/

            virtual void updateInterfaceMatrix
            (
                const scalarField& psiInternal,
                scalarField& result,
                const lduMatrix&,
                const scalarField& coeffs,
                const direction cmpt,
                const Pstream::commsTypes commsType
            ) const;
            
            /**
             * \brief Returns true if the corresponding polyPatch is a master.
             */
            virtual bool master() const
            {
                return domainScalingPatch_.master();
            }
            
            /**
             * \brief Returns true if this is a slave
             */
            virtual bool slave() const
            {
                return domainScalingPatch_.slave();
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "domainScalingFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
