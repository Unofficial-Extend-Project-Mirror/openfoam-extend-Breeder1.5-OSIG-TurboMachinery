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

#include "mixingPlaneFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "primitiveMixingPlanePatch.H"
#include "OFstream.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
mixingPlaneFvPatchField<Type>::mixingPlaneFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    mixingPlanePatch_(refCast<const mixingPlaneFvPatch>(p)),
    fieldName_(iF.name()),
    bzeroGradient_(false),
    adopted_(false)
{
  bextensive_ = false;  
 // Info << "In fucking constructor!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"  << endl;
}

template<class Type>
mixingPlaneFvPatchField<Type>::mixingPlaneFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict),
    mixingPlanePatch_(refCast<const mixingPlaneFvPatch>(p)),
    fieldName_(""),
    bzeroGradient_(false),
    adopted_(false),
    operatorTrigger_(false)
{
    if (!isType<mixingPlaneFvPatch>(p))
    {
        FatalIOErrorIn
        (
            "mixingPlaneFvPatchField<Type>::mixingPlaneFvPatchField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<Type, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict
        )   << "patch " << this->patch().index() << " not mixingPlane type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
    
    // checking if extensive
  
    
    if( dict.found("extensive") )
    {
        Switch extensive
        (
            dict.lookup("extensive")        
        );
        if(extensive)
            bextensive_ = true;            
    }
    else
    {
        bextensive_ = false;        
    }
    
   /* Info << fieldName_ 
         << "bextensive " << bextensive_ << endl;*/
    
    
    if( dict.found("rhoName") )
    {
      
        word rhoName( dict.lookup("rhoName") );
        rhoName_ = rhoName;
        Info << "RhoName " <<  rhoName_ << endl;
    }
    else
    {   
        rhoName_= word("false");
        Info << "RhoName " <<  rhoName_ << endl;
    }    
   
    
    word btest(dict.lookup("zeroGradient"));
    
    if( btest == "true" )
    {
        bzeroGradient_ = true;
        fvPatchField<Type>::operator=(this->patchInternalField());
    }
    else
    {
        bzeroGradient_ = false;
        
    }
    
    
//    Info << "In constructor "<< endl;
//    fvPatchField<Type>::operator=(this->patchInternalField()); // Dangerous
    
    
}

template<class Type>
mixingPlaneFvPatchField<Type>::mixingPlaneFvPatchField
(
    const mixingPlaneFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
    mixingPlanePatch_(refCast<const mixingPlaneFvPatch>(p)),
    fieldName_(""),
    bextensive_(ptf.extensive()),
    rhoName_(ptf.rhoName()),    
    bzeroGradient_(ptf.zeroGradient()),
    operatorTrigger_(false)
{
    if (!isType<mixingPlaneFvPatch>(this->patch()))
    {
        FatalErrorIn
        (
            "mixingPlaneFvPatchField<Type>::mixingPlaneFvPatchField\n"
            "(\n"
            "    const mixingPlaneFvPatchField<Type>& ptf,\n"
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
   // fvPatchField<Type>::operator= (this->patchInternalField());
   // Info << "In Copy constructor "<< endl;
}

template<class Type>
mixingPlaneFvPatchField<Type>::mixingPlaneFvPatchField
(
    const mixingPlaneFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    mixingPlanePatch_(refCast<const mixingPlaneFvPatch>(ptf.patch())),
    fieldName_(iF.name()),
    bextensive_(ptf.extensive()),
    rhoName_(ptf.rhoName()),    
    bzeroGradient_(ptf.zeroGradient()),
    operatorTrigger_(false)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class Type>
const mixingPlaneFvPatchField<Type>&  mixingPlaneFvPatchField<Type>::shadowPatchField() const
{
    
    fieldName_ = this->dimensionedInternalField().name();
    typedef GeometricField<Type, fvPatchField, volMesh> GeoField;
    
    const mixingPlaneFvPatchField<Type>& returnField =
    refCast<const mixingPlaneFvPatchField<Type> >
    (
            mixingPlanePatch_.shadow().lookupPatchField<GeoField,Type>(fieldName_)
    );        
    
    return ( returnField );
    
}   

template<class Type>
tmp<Field<Type> > mixingPlaneFvPatchField<Type>::snGrad() const
{    
    if ( bzeroGradient_ )
    {
    
        return tmp<Field<Type> >
        (
             new Field<Type>(this->size(), pTraits<Type>::zero)
        );

    }
    else
    {
        return (fvPatchField<Type>::snGrad());
    }
}

template<class Type>
void mixingPlaneFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{   
   
    if(!adopted_)
    {
        mixingPlanePatch_.adaptToMaster();
        adopted_ = true;
    }
    
    if ( bzeroGradient_ )
    {
        if (!this->updated())
        {
            this->updateCoeffs();
        }
        //Info << "In zeroGradient "<< fieldName_ << endl;

        fvPatchField<Type>::operator==(this->patchInternalField());
        fvPatchField<Type>::evaluate(Pstream::blocking);
    }
    else
    {
        List<Type> sliceAvg(1);
        Field<Type> helpField = shadowPatchField().makeCircumferentialAverage
                              ( this->internalField(), sliceAvg  );  
        
       
        const List<Type>& rsliceAvg(sliceAvg);
//        Info << "SliceAvg " << rsliceAvg << endl;
        Field<Type> firstField =
        this->makeCircumferentialAverage<Type>( this->internalField(), rsliceAvg  );

//        Info << firstField << endl;
       
        
        fvPatchField<Type>::operator=
        (
                firstField
        );
    }
}

template<class Type>
tmp<Field<Type> > mixingPlaneFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    if( bzeroGradient_ )
    {
        return tmp<Field<Type> >
          (
              new Field<Type>(this->size(), pTraits<Type>::one)
          );
    }
    else
    {
        return tmp<Field<Type> >
        (
            new Field<Type>(this->size(), pTraits<Type>::zero)
        );
    }
}

template<class Type>
tmp<Field<Type> > mixingPlaneFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    if( bzeroGradient_ )
    {
        return tmp<Field<Type> >
        (
            new Field<Type>(this->size(), pTraits<Type>::zero)
        );
    }
    else
    {
        
//        Info << *this << endl;
        return *this;
    }
    
}

template<class Type>
tmp<Field<Type> > mixingPlaneFvPatchField<Type>::gradientInternalCoeffs() const
{
    if( bzeroGradient_ )
    { 
        return tmp<Field<Type> >
        (
                new Field<Type>(this->size(), pTraits<Type>::zero)
        );
    }
    else
    {
        return -pTraits<Type>::one * this->patch().deltaCoeffs();
    }
}

template<class Type>
tmp<Field<Type> > mixingPlaneFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    if( bzeroGradient_  )
    {
        return tmp<Field<Type> >
        (
           new Field<Type>(this->size(), pTraits<Type>::zero)
        );
    }
    else
    {
        return this->patch().deltaCoeffs() * (*this);
    }
}

template<class Type> 
template<class Average>
tmp<Field<Average> > mixingPlaneFvPatchField<Type>::makeCircumferentialAverage
(
   const Field<Average>& iField,
   List<Type>& sliceAvg,
   bool bgenerate
) const
{
    
    bool bdebug = false;
   
    const unallocLabelList& fc = mixingPlanePatch_.faceCells();
   
//    OFstream oField("intField.txt",ios_base::app);
    
    Field<Average> mainField(fc.size(), pTraits<Average>::zero);
 
    
    typedef GeometricField<vector, fvPatchField, volMesh> GeoField;   
    const vectorField& USurfP = 
    refCast<const Field<vector> >
    (
        mixingPlanePatch_.lookupPatchField<GeoField,vector>("U")
    );  
    vectorField massField( USurfP.size() );
    
      
  // Info << *this << fieldName_ << endl;
  // Info << "RhoName " << rhoName_ << endl;
    if( rhoName_=="false" || rhoName_=="")        
    {            
        scalar srho(1.225);
        
        forAll(USurfP, i)
        {
             massField[i] = srho*USurfP[i];        
        }       
    }
    else
    {    
        
        const scalarField& rhoSurf = 
            refCast<const Field< scalar > >
        (
            mixingPlanePatch_.lookupPatchField<GeoField,scalar>("p") 
            //rhoName_)
        );
                
        if(USurfP.size() != rhoSurf.size())
        {
           // Info << "Error " << endl;
        }
        
     
        forAll(USurfP, i)
        {
             massField[i] = rhoSurf[i]*USurfP[i];        
        }
    }  

    forAll(mainField, i)
    {
        mainField[i] = iField[fc[i]];
    } 
    
    tmp<Field<Average> > tFirst
    (
            new Field<Average>
            (
                    mainField.size(),
                    pTraits<Average>::zero
            )
    );

    Field<Average> firstField = tFirst();


    const mixingPlanePolyPatch& rpolyPatch =
    mixingPlanePatch_.GetmixingPlanePolyPatch();

    const primitiveMixingPlanePatch& primRotPatch = rpolyPatch.RotPatch();
    
    
//    
  
//    Info <<  USurfP << endl;
    
    if(bextensive_ && bgenerate )
    {       
//        Info << *this << fieldName_ 
//             << "bextensive " << bextensive_ << endl;
        tFirst = const_cast<primitiveMixingPlanePatch&>(primRotPatch).MakeCircumferentialAverage
                            (  
                                    mainField,    
                                    mixingPlanePatch_.GetmixingPlanePolyPatch(),
                                    sliceAvg,
                                    massField
                             );
    }
    else
    {
    
    if( bgenerate )
    {
    
        
        tFirst = const_cast<primitiveMixingPlanePatch&>(primRotPatch).MakeCircumferentialAverage
                    (  
                            mainField,    
                            mixingPlanePatch_.GetmixingPlanePolyPatch(),
                            sliceAvg
                     );
    
    }
    else
    {
        tFirst = const_cast<primitiveMixingPlanePatch&>(primRotPatch).MakeCircumferentialAverage
                        (                                
                           mixingPlanePatch_.GetmixingPlanePolyPatch(),
                           sliceAvg
                         );   
        
    }
    }

    if(bdebug)
        Info << "In Neighbourfield Master " << fieldName_ << endl;    
    
    return (tFirst);
}

template<class Type> 
template<class Average>
tmp<Field<Average> > mixingPlaneFvPatchField<Type>::makeCircumferentialAverage
(
   const Field<Average>& iField,
   const List<Type>& sliceAvg

) const
{
    
   // bool bdebug = false;
   
    const unallocLabelList& fc = mixingPlanePatch_.faceCells();
   
//    OFstream oField("intField.txt",ios_base::app);
    
    Field<Average> mainField(fc.size(), pTraits<Average>::zero);

    forAll(mainField, i)
    {
        mainField[i] = iField[fc[i]];
    }
    
    
    tmp<Field<Average> > tFirst
    (
            new Field<Average>
            (
                    mainField.size(),
                    pTraits<Average>::zero
            )
    );

    Field<Average> firstField = tFirst();

    const mixingPlanePolyPatch& rpolyPatch =
    mixingPlanePatch_.GetmixingPlanePolyPatch();

    const primitiveMixingPlanePatch& primRotPatch = rpolyPatch.RotPatch();
    
 
   tFirst = const_cast<primitiveMixingPlanePatch&>(primRotPatch).MakeCircumferentialAverage
            (                                
                mixingPlanePatch_.GetmixingPlanePolyPatch(),
                sliceAvg
            );   
      
    
    return (tFirst);
}
    
template<class Type>
void mixingPlaneFvPatchField<Type>::operatorUpdate()
{
    operatorTrigger_ = false;
    this->evaluate(Pstream::blocking);  
}

template<class Type>
void mixingPlaneFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    
    if( bzeroGradient_ )
    {
        os.writeKeyword("zeroGradient") << "true"
            << token::END_STATEMENT << nl;
    }
    else
    {
        os.writeKeyword("zeroGradient") << "false"
           << token::END_STATEMENT << nl;
    }
    
    this->writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "mixingPlaneFvPatchFieldOperators.C"

// ************************************************************************* //
