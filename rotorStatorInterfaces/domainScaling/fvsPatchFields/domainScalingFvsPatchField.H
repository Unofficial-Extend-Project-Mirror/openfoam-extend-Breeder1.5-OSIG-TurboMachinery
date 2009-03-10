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
    domainScalingFvsPatchField

Description
    domainScaling Patch Field

Authors
    Franz Blaim

SourceFiles
    domainScalingFvsPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef domainScalingFvsPatchField_H
#define domainScalingFvsPatchField_H

#include "coupledFvsPatchField.H"
#include "domainScalingFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class domainScalingFvsPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class domainScalingFvsPatchField
:
    public coupledFvsPatchField<Type>
{
    // Private data
        
        /**
         * \brief Constant reference on the domainScalingFvPatch.
         */
        const domainScalingFvPatch& domainScalingPatch_;


public:

    /**
     * \brief Runtime type information
     */
    TypeName(domainScalingPolyPatch::typeName_());


    // Constructors

        /**
         * \brief Construct from patch and internal field
         */
        domainScalingFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&
        );

        /**
         * \brief Construct from patch, internal field and dictionary
         */
        domainScalingFvsPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const dictionary&
        );

        /**
         * \brief Construct by mapping given domainScalingFvsPatchField onto a new patch
         */
        domainScalingFvsPatchField
        (
            const domainScalingFvsPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, surfaceMesh>&,
            const fvPatchFieldMapper&
        );

        /**
         * \brief Construct as copy setting the internal field reference.
         */
        domainScalingFvsPatchField
        (
            const domainScalingFvsPatchField<Type>&,
            const DimensionedField<Type, surfaceMesh>&
        );

        /**
         * \brief Construct and return a clone.
         */
        virtual tmp<fvsPatchField<Type> > clone() const
        {
            return tmp<fvsPatchField<Type> >
            (
                new domainScalingFvsPatchField<Type>(*this)
            );
        }

        /**
         * \brief Construct and return a clone by setting internal field reference
         */
        virtual tmp<fvsPatchField<Type> > clone
        (
            const DimensionedField<Type, surfaceMesh>& iF
        ) const
        {
            return tmp<fvsPatchField<Type> >
            (
                new domainScalingFvsPatchField<Type>(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "domainScalingFvsPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
