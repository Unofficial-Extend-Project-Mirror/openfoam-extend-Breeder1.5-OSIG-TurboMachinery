/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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

#ifndef PrimitiveRotationalPatch_H
#define PrimitiveRotationalPatch_H

#include "patchIdentifier.H"
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




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam

{


//
//typedef PrimitivePatch<face, List, const pointField&>
//        primitiveAuxPatch;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                           Class polyPatch Declaration
\*---------------------------------------------------------------------------*/

class primitiveRotationalPatch
:
//    public patchIdentifier,
    public primitiveAuxPatch // TODO erase that inheritance
{
    protected:
        
        typedef struct
        {
            public:
            label FirstEdge;
            label SecondEdge;
            scalar Angle;
            label ConnectionPoint;
        }Corner; // Lists are constricted for use in FindCorners.
        
        /**
         *  \brief{defining the different turbomachine types  }
         *  Because the machine type  on the geometric
         * conditions of the interface it i  
         */       
        enum MachineType
        {
            Axial,
            Radial,
            Diagonal,
        };
        

        //- Start label of this patch in the polyMesh face list
        label start_;
        
        word name_;        
      

        //- Face-cell addressing
        mutable labelList::subList* faceCellsPtr_;

        //- Global edge addressing
        mutable labelList* mePtr_;

        /**
         *  \brief{All faces of the interpolation sector  }
         *  \todo{ rename it to originalInterpolationFaces }   
         */ 
        List<face> AllFaces_;

        /**
         *  \brief{points of the interpolation sector  }
         *  \todo{ rename it to originallocalPoints_ }         
         */         
        pointField localPoints_;
        
        /**
         *  \brief{reference to ther original primitive patch }         
         */ 
        const primitivePatch& OrigPatch_;

        /**
         *  \brief{the border points in cylindrical coordinates  }
         */         
        pointField CCSBorderPoints_;

        /**
         *  \brief{the corners of the patch}
         *  Those corners are determined by the FindCorners method.
         */         
        List<Corner> Corners_;

        /**
         *  \brief{points of the original patch}
         */         
        pointField AllPoints_;

        /**
         *  \brief{list of edges storing the border edges of the original patch}
         */
        edgeList BorderEdges_;
        
        /**
         *  \brief{list of edges storing the border edges of the original patch}
         */
        pointField BorderPoints_;
        
        /**
         *  \brief{list of edges storing the border edges of the original patch}
         */
        cylindricalCS cylindricalCS_;

        /**
         *  \brief{list of edges storing the border edges of the original patch}
         */        
        List<List<label> > pointLabels_;
        
        /**
         *  \brief{list of edges storing the border edges of the original patch}
         */
        List<List<label> > sliceFaces_;
        
        /**
         *  \brief{storing the type of the turbomachine radial diagonal or axial}
         *  This is determined in the method determineMachine
         */
        MachineType mType_;

        /**
         *  \brief{list of edges storing the border edges of the original patch}
         */        
        pointField pointsSide_;
        
        /**
         *  \brief{}
         * \todo{ move it to the mixingplane}
         */        
        pointField maxRbPoints_;
        
        /**
         *  \brief{storing the interpolation patch}
         *  \todo{this should be moved to the using class}
         */        
        primitiveAuxPatch*  primPatch_; 
        
        
            


    // Private Member Functions

        //- Calculate labels of mesh edges
        void calcMeshEdges() const; 
        
        /**
         * \brief{Devides the given patch into radialy spaced slices .}
         * \date{01-04-2008}
         * The function calls the GetRadialClipping the GetCircumferentialClipping
         * \todo{all the necessary steps which define the data for generating 
         *       the points and faces should be placed in here}
         */
        virtual void calcRadialGeometry(); //
        
        /**
         * \brief{Detects the two borders in radial direction.}
         * \date{01-04-2008}
         * \param <MinRB> {Border with the minimal radial values}
         * \param <MaxRB> {Border with maximal radial values }
         * \param <MinCB> {Border with minimal circumferential values}
         * \param <MaxCB> {Border with maximal circumferential values }
         *            
         *            MaxRB
         *         __---^--___
         *         |    |   |
         *   MinCB  |  R|  | MaxCB
         *           |  | |
         *           ------
         *              | MinRB
         *              x
         * It has to be noted that MinCB and MaxCB are dependent on the angle
         * measuring direction.
         * 
         * In case of the radial machine there is no MaxRB and MinRB than those
         * carry the maximal and minimal borders in z direction.
         * \todo{split it up into pieces too long}
         */
        void getBorders 
        (
          labelList& MinRB, 
          labelList& MaxRB,
          labelList& MinCB, 
          labelList& MaxCB
        );
        
        /**
         * \brief{generates the points for the interpolation sector.}
         * \date{12-09-2008}
         * \param{CC is the circumferentialclipping list}
         * \param{RC is the radialclipping list}
         * \param{SF is the Slice to Face mapping list}
         * The faces are ordered by slices, every slice has a label and that label 
         * is used as a key to get the ending number of the face.
         * slice with label 0 starts at face zero and ends a face f,
         * than the next slice starts with f+1 and ends with ff ...
         */
        void generatePointList
        (   
            const labelList& MaxRB,
            const labelList& MinCB            
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
        
        /**
         * \brief{Finds the corners of the original Patch}
         * \date{01-04-2008}
         * \Author{Franz Blaim}
         * The output is a edgelist wich delivers the pairs 
         * of edges defining a corner.
         * E.G. edgeList[2] is first edge of second corner
         * and edgeList[3] is the second edge of the second corner.
         * How it works :
   
         * 
         */
        void FindCorners();
        /*
        * \f$ \alpha = \frac{ \bold{v}_{1} \cdot w_{1} + }{}*\f$  
           * \alpha=\frac{\mathbf{v}_{1}\cdot\mathbf{w}_{1}+
           * \mathbf{v}_{2}\cdot\mathbf{w}_{2}+\mathbf{v}_{3}
           * \cdot\mathbf{w}_{3}}{\mathbf{\left|v\right|\cdot}
           * \left|\mathbf{w}\right|}*/     
        
        /**
         * \brief{Reorder the local edges}
         * This function generates a edgelist so that 
         * every edge is directly followed by its neighbour.        
         */
        edgeList Reorderedges();
        
        /**
         * \brief{Generates the external edges from the origininal Patch}
         * \param <rOrigPatch> {reference to the OriginalPatch}
         * \todo{check what is done in that}
         */
        void DefineLocalExternalEdges
        (
            const edgeList& extEdges, 
            const pointField& extPoints
        );
        
        /**
         * \brief{Generates the external edges from the origininal Patch}
         * \param <rOrigPatch> {reference to the OriginalPatch}
         */
        edgeList OrigExtEdges(const primitivePatch& rOrigPatch);
        
        /**
         * \brief{Generates the border points in cylindrical coordinatess}         
         */        
        void generateCCSBorderPoints();
        
        
        /**
         * \brief{rebuilds the rotational patch by using foreign edges as border.}
         * \date{22-05-2008}
         * The target of this function is to rebuild the rotational patch by using 
         * the distances between foreign points. So the original border will be conserved.
         * But within the original border the spacing between the points are adjusted too 
         * be equal with the foreign spacing.
         * A case discrimination has to be made between the three geometrical cases.
         * diagonal and radial: The z coordinate spacing has to be used for 
         *                      generting the new edge points.
         * axial              : The r coordinate has to be used for generating 
         *                      the new edge points  
         */
        List<point> calculateForeignMatchedBorder
        (
            const List<point>& edgePointList ,
            const labelList& MinCB 
            
        );
        
        
        MachineType determineMachine();
        
        
        /* * * * * * * * * * * * * * * * * * * * * * * * * * * */
        /*       Definition of the sort helper functions       */
        /* * * * * * * * * * * * * * * * * * * * * * * * * * * */
        
        /**
         * \brief{rebuilds the rotational patch by using foreign edges as border.}
         * \date{22-05-2008}
         * The target of this function is to rebuild the rotational patch by using 
         * the distances between foreign points. So the original border will be conserved.
         * But within the original border the spacing between the points are adjusted too 
         * be equal with the foreign spacing.
         * A case discrimination has to be made between the three geometrical cases.
         * diagonal and radial: The z coordinate spacing has to be used for 
         *                      generting the new edge points.
         * axial              : The r coordinate has to be used for generating 
         *                      the new edge points  
         */
         List< Corner > mergesort
         (
                 const List< primitiveRotationalPatch::Corner >& a              
         );
        
         
         List< Corner> merge
         (
                 List< primitiveRotationalPatch::Corner>& a, 
                 List< primitiveRotationalPatch::Corner>& b
         );

         void removeFirst
         (
                 List< primitiveRotationalPatch::Corner>& rT
         );
         
         
         List< scalar > mergesort(const List< scalar >& a);
         
         List< scalar > merge(List< scalar >& a, List< scalar >& b);
         
         void removeFirst(List< scalar>& rT);
         
         
         
         /**
          * \brief{rebuilds the rotational patch by using foreign edges as border.}
          * \date{22-05-2008}
          * The target of this function is to rebuild the rotational patch by using 
          * the distances between foreign points. So the original border will be conserved.
          * But within the original border the spacing between the points are adjusted too 
          * be equal with the foreign spacing.
          * A case discrimination has to be made between the three geometrical cases.
          * diagonal and radial: The z coordinate spacing has to be used for 
          *                      generting the new edge points.
          * axial              : The r coordinate has to be used for generating 
          *                      the new edge points  
          */
         
         List<point> calculateForeignMatchedPoints
         (
                 const List<point>& foreignPoints ,                 
                 MachineType mtype
         );
         
         
         virtual void makeSidePoints( const labelList& MinCB );
         
         virtual void makeMaxRBPoints(const labelList& MaxRB);       
         

    // Protected static data

        //- Relative tolerance (for geometric matching). Is factor of
        //  maximum edge length per face.
         /**
          * \todo{necessary}
          */
        static const scalar matchTol_;
    
    
    // Protected Member Functions
        // The primitiveRotationalPatch geometry initialisation is called by polyBoundaryMesh
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
        
        
        void writeOBJ(Ostream& os, const point& pt);
        
        static void writeOBJ
        (
             Ostream& os, 
             const point& p0, 
             const point& p1,
             label& vertI
        );
        
        virtual void clearAddressing();
        // Auxiliary patch analysis functionality

            //- Write selected points in OBJ format


public:

    //- Runtime type information
    // TODO Must be enabled when it is included as a dynamic library
    TypeName("PrimitiveRadialpatch"); 

    const primitivePatch& getOrigPatch()
    {
        return OrigPatch_;
    }

    // Constructors

        //- Construct from components
//        primitiveRotationalPatch
//        (
//            const List<face>& faces, // Schwierig da nicht ganz sauber muss von der abgeleiteten Klasse abhängen
//            const pointField& points,
//            const word& name
//        );
        
        
        /**
         * \brief{Constructor with additional clipping points}
         * \param<OrigPatch> {Original Patch for getting the bounding edges for this patch}
         * \param<PointsCorners> {additional clipping points used for specifying the corners of the patch}
         *         
         */        
        primitiveRotationalPatch
        (       
                const List<face>& faces, 
                const pointField& points,
                const primitivePatch& OrigPatch,
                const word& name
        );
        
        
        /**
         * \brief{Constructor with additional clipping points}
         * \param<OrigPatch> {Original Patch for getting the bounding edges for this patch}
         * \param<PointsCorners> {additional clipping points used for specifying the corners of the patch}
         *         
         */        
        primitiveRotationalPatch
        (       
              const List<face>& faces, 
              const pointField& points,
              const primitivePatch& OrigPatch,
              vector CCSCenter,
              vector CCSaxis,
              vector CCSdirection,
              const word& name
        );

        //- Construct as copy
        primitiveRotationalPatch(const primitiveRotationalPatch&);
   
        
    // Selectors

        //- Return a point er to a new patch created on freestore from
        //  components
        static autoPtr<primitiveRotationalPatch> New
        (
            const List<face>& faces,
            const pointField& points
        );
        
        // Destructor

        virtual ~primitiveRotationalPatch();
        
        // Access methods        
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
         void SetCylindricalCoordinateSystem(const cylindricalCS& ccs );
         
         /**
          *\brief{Returns a reference to the cylindrical coordinate System =}
          * The values of the reference cannnot be changed but can be accessed 
          * with the following functions : 
          *    name() returns a reference on the name of the CCS.
          *    origin() returns a reference on the vector defining the 
          *             origin of the CCS
          *    axis() returns a reference on the vector defining the axis.
          *    direction() returns a reference on the vector defining the 
          *                direction.   
          */
         const cylindricalCS& getCylindricalCoordinateSystem() const;
         
   
        // Member methods 
        /**
         * \brief{Generates a tensor for coordinate rotation}
         * \date{28-07-2008}
         * Usage the returned tensor can be muliplied by a dot product
         * with a vector or point. And this vector or point will be rotated 
         * around the cylindrical coordinate system axis with the given angle.
         * Note: In OF the & is the operator for dot products. 
         */
        tensor makeRotationTensor(const scalar& angle) const;

        /**
         * \brief{Generates a tensor for coordinate rotation}
         * \date{28-07-2008}
         * Usage the returned tensor can be muliplied by a dot product
         * with a vector or point. And this vector or point will be rotated 
         * around a given axis with the given angle.
         * Note: In OF the & is the operator for dot products. 
         */        
        tensor makeRotationTensor(vector& axis, scalar& angle);
        
        /**
         * \brief{tensor for transformation from cylindrical to rectangular coordinates}
         * \date{28-07-2008}
         * \param<Position>{the position of the face to which the field variable which 
         *                  should be transformed belongs}
         * \param<ccsVec>{vector in cylinder coordinates which must be transformed in 
         *                rectangular coordinates}
         */                
        vector transformFromCCSToCartesian( vector& ccsVec, point& Position, point& origin ) const;

        /**
         * \brief{tensor for transformation from cylindrical to rectangular coordinates}
         * \date{28-07-2008}
         * \param<Position>{the position of the face to which the field variable which 
         *                  should be transformed belongs}
         * \param<ccsVec>{vector in cylinder coordinates which must be transformed in 
         *                rectangular coordinates}
         */                        
        vector transformFromCCSToCartesian( vector& ccsVec, point& Position) const;
        
        

        /**
         * \brief{tensor for transformation from cylindrical to rectangular coordinates}
         * \date{28-07-2008}
         * \param<Position>{the position of the face to which the field variable which 
         *                  should be transformed belongs}
         */                        
        tensor tensorFromCCSToCartesian( point& Position ) const;
        

        /**
         * \brief{tensor for transformation from cylindrical to rectangular coordinates}
         * \date{28-07-2008}
         * \param<Position>{the position of the face to which the field variable which 
         *                  should be transformed belongs}
         */                        
        tensor tensorFromCartesianToCCS( point& Position ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

/// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
