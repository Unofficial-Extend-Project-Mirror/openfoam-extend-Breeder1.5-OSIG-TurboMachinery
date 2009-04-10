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
    mixingPlaneFvPatchField

Description
    Mixing Plane Patch Field

Authors
    Franz Blaim

SourceFiles
    mixingPlaneFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef mixingPlaneFvPatchField_H
#define mixingPlaneFvPatchField_H

#include "fvPatchField.H"
#include "mixingPlaneFvPatch.H"
#include "transform.H"
#include "scalarField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class mixingPlaneFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class mixingPlaneFvPatchField
:
   public fvPatchField<Type>
{
    // Private data
       
       const mixingPlaneFvPatch& mixingPlanePatch_; 
        
       /**
        * \brief{ stores the name of the field }
        * \date{20-08-2008}
        */
       mutable word fieldName_;
       
       /**
        * \brief{ defines if the field is a extensive state variable }
        * \date{20-08-2008}
        */
       bool bextensive_;
       
       /**
        * \brief{ stores the name of the density variable }
        * \date{20-08-2008}
        * This name is used to retrieve the density surface scalar field which is used
        * to generate a mass weighted average.
        */
       word rhoName_;
       
       /**
        * \brief{ defines if the boundary is handles as zerogradient condition }
        * \date{20-08-2008}
        */
       bool bzeroGradient_;
       
       /**
        * \brief{ defines if the geometry of the underlying patch has been adopted correctly }
        * \date{20-08-2008}
        */
       bool adopted_;
       
       bool operatorTrigger_;
       
       
      
public:

    //- Runtime type information
    TypeName("mixingPlaneFvPatchField");


    // Constructors

        //- Construct from patch and internal field
        mixingPlaneFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary        
        mixingPlaneFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mixingPlaneFvPatchField onto a new patch
        mixingPlaneFvPatchField
        (
            const mixingPlaneFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        mixingPlaneFvPatchField
        (
            const mixingPlaneFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new mixingPlaneFvPatchField<Type>(*this)
            );
        }

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new mixingPlaneFvPatchField<Type>(*this, iF)
            );
        }
        

    // Member functions

        // derived functions   
        
            /* \brief{This function is used to generate the circumferential average }
             * Important function otherwise instabilities will occur
             */
            virtual bool fixesValue() const
            {
                 if (bzeroGradient_)
                 {
                    return false;
                 }
                 else
                 {
                    return true;
                 }
             }
        
            /**
             * \brief{This function is used to generate the circumferential average }
             * \date{20-08-2008}
             */      
             //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType
            );
            
            /**
             * \brief{This function is used to generate the circumferential average }
             * \date{20-08-2008}
             * TODO FB check for what they are used theoretical background
             */  
            virtual tmp<Field<Type> > valueInternalCoeffs
            (
                const tmp<scalarField>& w
            ) const;
            
            /**
             * \brief{This function is used to generate the circumferential average }
             * \date{20-08-2008}
             * TODO FB  check for what they are used theoretical background
             */  
            virtual tmp<Field<Type> > valueBoundaryCoeffs
            ( 
                const tmp<scalarField>& w 
            ) const;
            
            
            /**
             * \brief{This function is used to generate the circumferential average }
             * \date{20-08-2008}
             * TODO FB check for what they are used theoretical background
             */  
            virtual tmp<Field<Type> > gradientInternalCoeffs() const;
            
            
            /**
             * \brief{This function is used to generate the circumferential average }
             * \date{20-08-2008}
             * TODO FB check for what they are used theoretical background
             */  
            virtual tmp<Field<Type> > gradientBoundaryCoeffs() const;
            
            
            void write(Ostream& os) const;
            
            
            /**
             * \brief{This function is used to generate the circumferential average }
             * \date{20-08-2008}
             */      
          
            
            virtual tmp<Field<Type> > snGrad() const;
            
            virtual bool master() const
            {
                return mixingPlanePatch_.master();
            }
            
            /**
             * \brief{ Is this the slave side? }
             */
            virtual bool slave() const
            {
                return mixingPlanePatch_.slave();
            }
            
            // TODO  it is necessary but for what must it be public  ??
            void operatorUpdate();
            
       // Special functions
            
            const mixingPlaneFvPatchField<Type>& shadowPatchField() const;
            
            //TODO check if one is obsolete
            template<class Average>
            tmp<Field<Average> > makeCircumferentialAverage
            (
               const Field<Average>& iField, 
               const List<Type>& sliceAvg
            ) const;           
                      
            template<class Average>
            tmp<Field<Average> > makeCircumferentialAverage
            (
                const Field<Average>& iField, 
                List<Type>& sliceAvg, 
                bool bgenerate=true 
            ) const;            
                  
            // TODO FB when is that function called necessary ??
            bool zeroGradient() const
            {
                return bzeroGradient_;
            }
            
            // TODO FB Is that function necessary 
            word rhoName() const
            {
                 return rhoName_;
            }
            
            // TODO FB Is that function necessary
            bool extensive() const
            {
                 return bextensive_;
            }
            
            // TODO FB is that necessary
            bool triggered() const
            {
                return operatorTrigger_;
                
            }
            
            
            /**
             * \brief{To fix the value this operator has been disabled}
             * \date{20-08-2008}
             */
            virtual void operator=(const UList<Type>& ul);          

            /**
             * \brief{To fix the value this operator has been disabled}
             * \date{20-08-2008}
             */
            virtual void operator=(const fvPatchField<Type>& ptf);
            
            /**
             * \brief{To fix the value this operator has been disabled}
             * \date{20-08-2008}
             */
            virtual void operator+=(const fvPatchField<Type>& ptf);
            
            /**
             * \brief{To fix the value this operator has been disabled}
             * \date{20-08-2008}
             */
            virtual void operator-=(const fvPatchField<Type>& ptf);
            
            /**
             * \brief{To fix the value this operator has been disabled}
             * \date{20-08-2008}
             */  
            virtual void operator*=(const fvPatchField<scalar>& ptf);
            
            /**
             * \brief{To fix the value this operator has been disabled}
             * \date{20-08-2008}
             */  
            virtual void operator/=(const fvPatchField<scalar>& ptf);
                        
            /**
             * \brief{To fix the value this operator has been disabled}
             * \date{20-08-2008}
             */           
            virtual void operator+=(const Field<Type>& tf);
            
            /**
             * \brief{To fix the value this operator has been disabled}
             * \date{20-08-2008}
             */  
            virtual void operator-=(const Field<Type>& tf);
            
            /**
             * \brief{To fix the value this operator has been disabled}
             * \date{20-08-2008}
             */  
            virtual void operator*=(const Field<scalar>&);
            
            /**
             * \brief{To fix the value this operator has been disabled}
             * \date{20-08-2008}
             */
            virtual void operator/=(const Field<scalar>& tf);
            
            /**
             * \brief{To fix the value this operator has been disabled}
             * \date{20-08-2008}
             */
            virtual void operator=(const Type& t);
            
            /**
             * \brief{To fix the value this operator has been disabled}
             * \date{20-08-2008}
             */
            virtual void operator+=(const Type& t);
            
            /**
             * \brief{To fix the value this operator has been disabled}
             * \date{20-08-2008}
             */
            virtual void operator-=(const Type& t);
            
            /**
             * \brief{To fix the value this operator has been disabled}
             * \date{20-08-2008}
             */
            virtual void operator*=(const scalar s);
            
            /**
             * \brief{To fix the value this operator has been disabled}
             * \date{20-08-2008}
             */
            virtual void operator/=(const scalar s);
            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "mixingPlaneFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
