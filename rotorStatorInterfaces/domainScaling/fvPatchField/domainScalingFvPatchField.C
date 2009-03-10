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

Authors
    Franz Blaim

\*---------------------------------------------------------------------------*/

#include "domainScalingFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "primitiveDomainScalingPatch.H"
#include "coupledFvPatchFieldsFwd.H"
#include "OFstream.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
domainScalingFvPatchField<Type>::domainScalingFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    domainScalingLduInterfaceField(),
    coupledFvPatchField<Type>(p, iF),
    domainScalingPatch_(refCast<const domainScalingFvPatch>(p)),
    fieldName_(iF.name())
{
        
}

template<class Type>
domainScalingFvPatchField<Type>::domainScalingFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    domainScalingLduInterfaceField(),
    coupledFvPatchField<Type>(p, iF, dict),
    domainScalingPatch_(refCast<const domainScalingFvPatch>(p)),
    fieldName_(iF.name())
{
    if (!isType<domainScalingFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "domainScalingFvPatchField<Type>::domainScalingFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not domainScalingFvPatchField type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }

    // Same as in the cyclic interface
    this->evaluate(Pstream::blocking);
}

template<class Type>
domainScalingFvPatchField<Type>::domainScalingFvPatchField
(
    const domainScalingFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    domainScalingLduInterfaceField(),
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    domainScalingPatch_(refCast<const domainScalingFvPatch>(p)),
    fieldName_("")
{
    if (!isType<domainScalingFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "domainScalingFvPatchField<Type>::domainScalingFvPatchField\n"
            "(\n"
            "    const domainScalingFvPatchField<Type>& ptf,\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const fvPatchFieldMapper& mapper\n"
            ")\n"
        )   << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}

template<class Type>
domainScalingFvPatchField<Type>::domainScalingFvPatchField
(
    const domainScalingFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    domainScalingLduInterfaceField(),
    coupledFvPatchField<Type>(ptf, iF),
    domainScalingPatch_(refCast<const domainScalingFvPatch>(ptf.patch())),
    fieldName_(iF.name())   
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class Type>
const domainScalingFvPatchField<Type>&  domainScalingFvPatchField<Type>::shadowPatchField() const
{
     fieldName_ = this->dimensionedInternalField().name();
     
     Info << "fieldName_ " << fieldName_ << endl;
     
      typedef GeometricField<Type, fvPatchField, volMesh> GeoField;
        
      const domainScalingFvPatchField<Type>& returnField =
      refCast<const domainScalingFvPatchField<Type> >
      (
             domainScalingPatch_.shadow().lookupPatchField<GeoField,Type>(fieldName_)
      );        
       
     return ( returnField );
        
}

template<class Type>
tmp<Field<Type> > domainScalingFvPatchField<Type>::patchNeighbourField() const
{     
    const domainScalingFvPatch& shadowPatch = domainScalingPatch_.shadow();
    
    tmp<Field<Type> > tresult
    (
            new Field<Type>
            (
                    shadowPatch.size(),
                    pTraits<Type>::zero
            )
    );
    
    Field<Type>&  firstField = tresult() ;
    
    const Field<Type>& iField =  this->internalField();

    const unallocLabelList& sfc = shadowPatch.faceCells();
    
    forAll(firstField, i)
    {
        firstField[i] = iField[sfc[i]];
    }
    
  //  Info << firstField << endl;
    firstField = shadowPatch.interpolateToShadow(firstField);
   // Info << firstField << endl;
    return (tresult);
    
}

template<class Type>
void domainScalingFvPatchField<Type>::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField&,
    const lduMatrix&,
    const scalarField&,
    const direction,
    const bool bufferdTransfer
) const
{
    if(Pstream::parRun())
    {
    
      labelList globalIndex;

      //this should be global size
      scalarField globalField(domainScalingPatch_.shadow().size(),0.0);

    
      scalarField neighbourField = mag(patchNeighbourField());
          
      forAll(neighbourField , facei)
      {
          label globalface = globalIndex[facei];
          globalField[globalface] = neighbourField[facei];
      }

      for(
               int slave = 0;
               slave <= Pstream::lastSlave();
                slave++
         )
      {
         
         //OPstream toProc(slave);

         //toProc << globalField;
             
      }
   } 
    
   // Set delta time to calculate the right position
   const fvMesh& mesh = this->dimensionedInternalField().mesh();
   scalar deltaT = mesh.time().timeOutputValue();//.value();
  
   const_cast<domainScalingFvPatch&>(domainScalingPatch_).setDeltaTime(deltaT,mesh.time());
    
}

// Return matrix product for coupled boundary
template<class Type>
void domainScalingFvPatchField<Type>::updateInterfaceMatrix
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    if(Pstream::parRun())
    {
        labelList globalIndex;

        List<scalarField> globalFields(Pstream::nProcs()); // nProcs is number of processors

        for(
                int slave = 0;
                slave <= Pstream::lastSlave();
                slave++
            )
         {
           
             //IPstream fromProc(slave);
             //fromProc >> globalFields[slave];
           
         }

        scalarField globalField(globalFields[0].size(),0.0);
 
        forAll(globalField, facei)
        {
           forAll(globalFields, procI)
           {
              globalField[facei] += globalFields[procI][facei];
           }
        }
        
        labelListList faceAddr;
        scalarListList faceWeights;
        
        if(domainScalingPatch_.master())
        {
           // faceAddr = read masterfaceaddressing
            
           // faceWeights = read masterfaceweights
        }
        else
        {
           // faceAddr = read slavefaceaddressing
            
           // faceWeights = read slavefaceweights
        }

        scalarField pnf(this->size(),0.0);
        
        forAll(pnf,facei)
        {
           label glI(globalIndex[facei]);

           forAll(faceAddr[glI], neighFacei)
           {
               scalar faceValue = globalField[faceAddr[glI][neighFacei]];

               scalar weight = faceWeights[glI][neighFacei];

               pnf[facei] +=  faceValue * weight; 
           }
        }

        // Transform according to the transformation tensor
       //  transformCoupleField(pnf, cmpt);

        // Multiply the field by coefficients and add into the result
        const unallocLabelList& fc = domainScalingPatch_.faceCells();

        forAll(fc, elemI)
        {
            result[fc[elemI]] -= coeffs[elemI]*pnf[elemI];
        }
    }
    else
    {
        // Get shadow face-cells and assemble shadow field
        const unallocLabelList& sfc = domainScalingPatch_.shadow().faceCells();

        scalarField sField(sfc.size());

        forAll (sField, i)
        {
            sField[i] = psiInternal[sfc[i]];
        }

        scalarField pnf = domainScalingPatch_.shadow().interpolateToShadow(sField,cmpt,this->rank());
        
        // Multiply the field by coefficients and add into the result
        const unallocLabelList& fc = domainScalingPatch_.faceCells();

        forAll(fc, elemI)
        {
            result[fc[elemI]] -= coeffs[elemI]*pnf[elemI];
        }
    }
}

template<class Type>
void domainScalingFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    Field<Type>::operator=
    (
       this->patch().weights()*this->patchInternalField()
       + (1.0 - this->patch().weights())*this->patchNeighbourField()
    );    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
