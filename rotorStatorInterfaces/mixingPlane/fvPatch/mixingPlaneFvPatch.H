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
    mixingPlaneFvPatch

Description
    The finiteVolumePatch for the mixing plane

Authors
    Franz Blaim

SourceFiles
    ixingPlaneFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef MixingPlaneFvPatch_H
#define MixingPlaneFvPatch_H

#include "coupledFvPatch.H"
#include "mixingPlanePolyPatch.H"
#include "primitiveMixingPlanePatch.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class MixingPlaneFvPatch Declaration
\*---------------------------------------------------------------------------*/
 

class mixingPlaneFvPatch
:
    public fvPatch
{
    // Private Data
        
      
public:
        struct Storage 
        { 
              scalarField Coeffs;
              scalarField Field;
              int iCounter;
        };
        
private: 
        
    
        //mutable Storage storage_;
        
        //- Reference to polyPatch        
        const mixingPlanePolyPatch& mixingPlanePolyPatch_;

protected:

        // Protected Member functions

        //- Make patch face - neighbour cell distances
        virtual void makeDeltaCoeffs(scalarField&) const;
       

        //- Make patch face non-orthogonality correction vectors
        virtual void makeCorrVecs(vectorField&) const
        {}

        //- Initialise the patches for moving points
        virtual void initMovePoints()
        {}

        //- Correct patches after moving points
        virtual void movePoints()
        {}



public:

    //- Runtime type information
    TypeName(mixingPlanePolyPatch::typeName_());
    
      // Constructors

        //- Construct from components
        mixingPlaneFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            fvPatch(patch, bm),
            mixingPlanePolyPatch_(refCast<const mixingPlanePolyPatch>(patch))           
        {
            Info << "fvPatch" << endl;
        }


    // Destructor

        virtual ~mixingPlaneFvPatch()
        {}


    // Member functions

        // Access
            /**
            * \brief{ Return slave patch index.}
            */
            int slaveIndex() const
            {
                return mixingPlanePolyPatch_.shadowIndex();
            }

            /**
             * \brief{ Returns wether this is a master or not.}
             */
            bool master() const
            {
                return mixingPlanePolyPatch_.master();
            }


            /**
             * \brief{ Returns wether this is the slave side.}
             */
            bool slave() const
            {
                return !master();
            }
            
            
            // This checks if the slave patch has really setted a master
            bool checkDoubleMaster() const
            {
                if(master())
                {
                   return mixingPlanePolyPatch_.shadowPatch().master();
                }
                else
                {
                   return false;
                }
            }
    
        word getNameOmega() const
        {
            return mixingPlanePolyPatch_.getNameOmega();
        }

            /**
            * \brief{ Return slave patch.}
            */
            const mixingPlaneFvPatch& shadow() const;

        // Interface transfer functions
            
            /**
            * \brief{ Return slave patch.}
            */
            virtual tmp<labelField> interfaceInternalField
            (
                const unallocLabelList& internalData
            ) const;
            
            /**
            * \brief{ Transfer and return neighbour field.}
            */
            virtual tmp<labelField> transfer
            (
                const unallocLabelList& interfaceData
            ) const;

            /**
            * \brief{ Calculate the delta distance for this patch.}
            */
            virtual tmp<vectorField> delta() const;
            
            /**
            * \brief{ Return neighbour field.}
            * 
            */
            virtual tmp<labelField> internalFieldTransfer
            (
                const unallocLabelList& internalData
            ) const;
               
            
            /**
             * \brief{returns ref on mixingPlanePolyPatch.}
             * \date{22-05-2008}
             * The function returns the mixingPlanePolyPatch member variable.
             */
            const mixingPlanePolyPatch& GetmixingPlanePolyPatch() const
            {
                return (mixingPlanePolyPatch_);
            }
            
            /**
             * \brief{returns ref on rotational patch.}
             * \date{21-05-2008}
             * The function returns the Rotational primitive patch .
             */
            const primitiveRotationalPatch& RotPatch() const
            {
                const primitiveRotationalPatch& rotPatch = mixingPlanePolyPatch_.RotPatch();
                return ( rotPatch );                
            }
            
            void adaptToMaster() const
            {   
                mixingPlanePolyPatch_.adaptToMaster();
            }

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
            
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
