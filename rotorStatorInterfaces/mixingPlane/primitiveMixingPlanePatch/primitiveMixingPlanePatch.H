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
    primitiveRadialPatch

Description
    A patch is a list of labels, which address the faces in the global face
    list. Based on the global faces, the patch can calculate its own edge.
    Patch also contains all addressing between the faces.

SourceFiles
    primitiveRadialPatch.C    

\*---------------------------------------------------------------------------*/

#ifndef primitiveMixingPlanePatch_H
#define primitiveMixingPlanePatch_H

#include "fvPatchFields.H"
#include "patchIdentifier.H"
#include "primitiveRotationalPatch.H"
#include "primitivePatch.H"
#include "primitiveAuxPatch.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "labelList.H"
#include "PrimitivePatch.H"
#include "face.H"
#include "pointField.H"
#include "cylindricalCS.H"
#include "ggiInterpolation.H"
#include "GGIInterpolation.H"




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                           Class polyPatch Declaration
\*---------------------------------------------------------------------------*/

class primitiveMixingPlanePatch
:
    public primitiveRotationalPatch
{
    public:        
        typedef GGIInterpolation
        <
             primitiveAuxPatch,
             primitivePatch
        >    ggiMMInterpolation;

    private:    
    // Private Member Functions

        ggiMMInterpolation* masterToPatchPtr_;
        
        /**
         * \brief{Devides the given patch into radialy spaced slices .}
         * \date{01-04-2008}
         * The function calls the GetRadialClipping the GetCircumferentialClipping
         */
        virtual void calcRadialGeometry(); //
        
        
        List< List<tensor> > tensorListList_; 
        
        void generatePointList
        (   
            const labelList& MaxRB,
            const labelList& MinCB            
        );  
        
        void generateForeignMatchedPointList
        (
           const List<point>& matchedSidePoints, 
           const labelList& MinCB
        );
        
        List<point> calculateForeignMatchedPoints
        (
                const List<point>& foreignPoints ,                 
                MachineType mtype
        );
        
        
        /**
         * \brief{Generates the new FaceList grouped in radially spaced slices .}
         * \date{01-04-2008}
         * \param{CC is the circumferentialclipping list}
         * \param{RC is the radialclipping list}
         * \param{SF is the Slice to Face mapping list}
         * The faces are ordered by slices, every slice has a label and that label 
         * is used as a key to get the ending number of the face.
         * slice with label 0 starts at face zero and ends a face f,
         * than the next slice starts with f+1 and ends with ff ...
         */
        virtual void generateNewFaceList
        (          
        );
        
        void makeMaxRBPoints(const labelList& MaxRB);
        
        void makeSidePoints( const labelList& MinCB );
        
        
        
        
    
    // Protected Member Functions


//        //- Return the polyMesh allPoints
//        const pointField& allPoints() const;

        // The primitiveMixingPlanePatch geometry initialisation is called by polyBoundaryMesh
        friend class polyBoundaryMesh;

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry()
        {}

        //- Calculate the patch geometry
        virtual void calcGeometry()
        {}

        //- Initialise the patches for moving points
        virtual void initMovePoints(const pointField&)
        {}

        //- Correct patches after moving points
        virtual void movePoints(const pointField& p);
        
        template<class Type > 
        Field<Type> Average
        (
             const Field<Type>& tpf
        )const ;
             
        template<class Type > 
        Field<Type> Average
        (
            const Field<Type>& tpf,
            List<Type>& sliceAvg
        )const ;
             
        template<class Type > 
        Field<Type> Average
        (
            const List<Type>& sliceAvg
        )const;
        
        
        template<class Type > 
        Field<Type> Average
        (
                  const Field<Type>& tpf,
                  List<Type>& sliceAvg,
                  scalarField massField
        )const; 
        
        Field<vector> Average
        (
             const Field<vector>& tpf,
             List<vector>& sliceAvg,
             scalarField massField
        )const; 
        
        Field<vector> Average
        (
            const Field<vector>& tpf,
            List<vector>& sliceAvg
        )const ;        
        
        Field<vector> Average
        (
            const List<vector>& sliceAvg
        );
        
public:

    //- Runtime type information
    // TODO Must be enabled when it is included as a dynamic library
    TypeName("primitiveMixingPlanePatch"); 


   
    
    // Constructors

        //- Construct from components
        primitiveMixingPlanePatch
        (
            const List<face>& faces, 
	        const pointField& points,
            const primitivePatch& OrigPatch,
            const point& CCSCenter,
            const vector& CCSaxis,
            const vector& CCSdirection,
            const word& name,
            bool blockMeshDict
        );
        
        //- Construct as copy
        //TODO Copyconstructor must be correctly implemented
        primitiveMixingPlanePatch(const primitiveMixingPlanePatch&);
  
        
        const List<face>& GetFaces() const;        
        
        const pointField& GetFacePoints() const;
             
        void getFacePoints(pointField& PF) const;
       
        void getFacePoints(List<face>& Faces) const;
        
    // Selectors

        //- Return a point er to a new patch created on freestore from
        //  components
        static autoPtr<primitiveMixingPlanePatch> New
        (
            const List<face>& faces,
            const pointField& points
        );
        
        // Destructor

        virtual ~primitiveMixingPlanePatch();
        
        
        const List<point>& getPointsSide() const;
       
        void rebuildPatchForEdges( const List<point>& edgePointList );
        
        template<class Type, class Patch> 
        tmp<Field<Type> > MakeCircumferentialAverage
        (
               const Field<Type>& masterField, 
               const Patch& masterPatch,
               List<Type>& sliceAvg,
               const vectorField& masterFieldMass
        ) ;
             
        
        // Member functions
        template<class Type, class Patch>
        tmp<Field<Type> > MakeCircumferentialAverage
        (
               const Field<Type>& masterField, 
               const Patch& masterPatch,
               List<Type>& sliceAvg
        );
               
        template<class Type, class Patch> 
        tmp<Field<Type> > MakeCircumferentialAverage
        (              
             const Patch& masterPatch,             
             const List<Type>& sliceAvg                 
        );
                                
                
                
        template<class Type, class Source, class Target> 
        tmp<Field<Type> > InterpolateBetweenPatches
        (
              const Field<Type>& tpf, 
              const Source& sourcePatch,
              const Target& targetPatch
        ) const;
        
        void generateCCSTransformationTensors();
        
        
        /**
         *\brief{ Advices the cylindrical coordinate system (CCS)}
         * Just a setter for the cylindrical coordinate system (CCS) 
         * which is used to generate the radialy spaced faces.
         * A cylindrical coordinate system can be generated by:
         * 1. including the cylindricalCS.H file
         * 2. constructing a cylindrical coordinate system
         *    Example :
         *    cylindricalCS ccs
         *    ( "ccs", // Name of cylindrical coordinate system
         *       vector(0, 1, 0), //center point of ccs
         *       vector(0, 0, 1), //axis of ccs 
         *       vector(1, 0, 0) //direction of the coordinatesystem 
         *                       // defining the angle 
         *     );
         *   Note : Bove the base axis for cylindrical angle and 
         *          the axis of the ccs must be perpendicular on each 
         *          other 
         *          
         */
         void setCylindricalCoordinateSystem(const cylindricalCS& ccs );         
    
         
         void resetInterpolation()
         {
            deleteDemandDrivenData(masterToPatchPtr_);
         }

         
         void writePatches( const word& patchName );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#      include "primitiveMixingPlanePatchTemplates.C"
#endif

#endif

// ************************************************************************* //
