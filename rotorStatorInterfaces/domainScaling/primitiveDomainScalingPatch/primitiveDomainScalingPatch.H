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
    primitiveDomainScalingPatch

Description
    This patch is used to generate and update the interpolation disk which is 
    used to provide the adressing for the fully 360 degree spanning field which is
    generated also within that patch.
    This patch can also be used to directly interpolate between the original patch 
    and the interpolation disk. 
    Therefore the main logic of the domainScaling interface is implemented within this 
    class. 

SourceFiles
    primitiveRadialPatch.C    

\*---------------------------------------------------------------------------*/

#ifndef primitiveDomainScalingPatch_H
#define primitiveDomainScalingPatch_H

#include "patchIdentifier.H"
#include "primitivePatch.H"
#include "primitiveRotationalPatch.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "labelList.H"
#include "PrimitivePatch.H"
#include "face.H"
#include "pointField.H"
//#include "cylindricalCS.H"
#include "ggiInterpolation.H"
#include "GGIInterpolation.H"
#include "Time.H"




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

class primitiveDomainScalingPatch
:
    public primitiveRotationalPatch
{       
        /**
         *  \brief Typdef for the intepolation pointers used in this class 
         * \date 31-08-2008
         * This is a template instantiation of the GGIInterpolation class
         * which is specialy suited for the use with a custom designed primitiveAuxPatch
         * which is a typdef in the primitiveRotationalPatch.
         */        
        typedef GGIInterpolation
        <
              primitiveAuxPatch,
              primitivePatch
        >    ggiMMInterpolation;       
       
        /**
         *  \brief This pointField contains the points for the interpolation disk faces 
         */
        pointField newOrigPoints_;
        
        
        /**
         *  \brief The facelist containing the interpolation disk faces  
         */        
        faceList newOrigPatch_;
        
        /**
         *  \brief  Storing the opening angle of the underlying patch 
         *  This variable is used to store the angle difference between the 
         *  corner with the maximum angle and the corner with the minimum angle.
         *  This is done within the calculateRotationalDivisor by getting the MinCB phi values
         *  and the MaxCB values
         * TODO FB calculation must be improved   
         */
        scalar diffPhi_; 
        
        /**
         * \brief  the timestep difference
         *          
         * Stores the difference between the actual timestep and 
         * the old timestep this difference is used to rotate the 
         * points of the patch
         * \date 31-08-2008
         */        
        scalar deltaTime_;
        
        /**
         * \brief The angular velocity in rpm
         *          
         * This member variable stores the angular velocity 
         * which rotates the interpolation disk.
         * \date 31-08-2008 
         */      
        scalar omega_;
        
        /**
         * \brief The timestep difference
         *   
         * Is used to store the actual timestep, which is used 
         * to detect if the patch must be rotated or not.
         * \date 31-08-2008
         */         
        scalar timeStep_;
        
        /**
         *  \brief The starting angle in degrees
         * 
         * This member variable is used to set the initial  
         * angular position for the interpolation disk.
         * In combination with the omega set to zero this leads to 
         * a frozen rotor approach
         * \date 31-08-2008 
         */              
        scalar startAngle_;
        
        Switch master_;
        
        
        /**
         * \brief Stores the interpolation patch to patch pointer 
         * \date 31-08-2008
         * The master means the original patch and the Patch means
         * the interpolation disk patch.
         */        
        ggiMMInterpolation* masterToPatchPtr_;
        
        /**
         * \brief This rotation tensor  is used to rotate the points in case of a moving interpolation disk 
         * \date 31-08-2008
         */        
        tensor* rotationTensor_;
        
        /**
         *  \brief This tensor is used to increment the rotation of the field
         * after every timestep. Was a experimental idea for the alternative 
         * handling of cylindrical shaped patch geometries, 
         * \date 31-08-2008
         */                
        tensor* incrementTensor_;
        
        /**
         * \brief This list of tensors is used to transform the field values
         *  at the master patch into the cylindrical coordinate frame.
         * \date 31-08-2008                  
         */        
        List<tensor> tensorList_;
        
        /**
         * \brief This list of tensors is used to transform the field after 
         *  interpolation to the slave. 
         * \date 17-09-2008
         * 
         */        
        List<tensor>  slaveTensorList_;
        

    
        /**
         * \brief{Devides the given patch into radialy spaced slices .}
         * \date{01-04-2008}
         * The function calls the GetRadialClipping the GetCircumferentialClipping
         */
        virtual void calcRadialGeometry(); //
        
        
        /**
         * \brief{Devides the given patch into radialy spaced slices .}
         * \date{03-06-2008}
         * The method extracts the points of the original patch 
         * to easily use the labels of the newOrignal patch which 
         * has been generated.
         */
        void extractNewOrigPatch();
        
        /**
         * \brief{method inherited from the primitiveRotationalPatch.}
         * \date{03-06-2008}
         * This method is used to generate the points of the interpolation disk
         * after the constructor is called
         */              
        void generatePointList
        (
                const labelList& MaxCB,
                const labelList& MinCB 
        );
        
        
        /**
         * \brief{Generates the new faceList.}
         * \date{31-08-2008}
         * This method generated the new faces of the interpolation disk patch.
         */
        void generateNewFaceList
        (          
        );        
        
        /**
         * \brief{calculates the opening angle of the patch .}
         * \date{31-08-2008}
         * \param<MinCB>{The border containing the minimum angle}
         * \param<MaxCB>{The border containing the maximum angle}
         * This method calculates the difference between the maximum circumferential angle
         * and the minimum circumferential angle line
         * TODO FB check if the value is independent from a angle averaging.
         *      Check if always the maximum angle is retrieved better
         *      a check must be introduced to check wether the maximum angle remains a
         *      constant because the whole thing should not change its angle value otherwise it will be 
         *      hard to generate a 360 degree disk 
         */        
        scalar calculateRotationalDivisor
        (
                const labelList& MaxCB, 
                const labelList& MinCB
        );
        
        
        /**
         * \brief{Moves the points if the timestep has been increased .}
         * \date{31-08-2008}
         * This method calculates the rotationTensor_ member variable
         * by using the deltaTime_ amd the omega_ member variable 
         * It has to be noted that the retrieving function for the deltaTime_
         * must be called before the points can be moved
         */ 
        void calcuteMovedPoints();
        
        /**
         * \brief{Used to distribute the values other the interpoaltion disk .}
         * \date{31-08-2008}
         * \param{sourceField the original field which will be distributed}
         * This function is used for dealing with the interpolation of field
         * values
         */  
        tmp<Field<vector> > cloneField
        (
           const Field<vector>& sourceField           
        ) const;
        
        /**
         * \brief{Used to distribute the values other the interpoaltion disk .}
         * \date{31-08-2008}
         * \param{sourceField the original field which will be distributed}
         * This function is used for dealing with the interpolation of field
         * values
         */    
        tmp<Field<vector> > cloneField
        (
             const tmp<Field<vector> >& sourceField             
        ) const;       
        
        /**
         * \brief{Used to distribute the values other the interpoaltion disk .}
         * \date{31-08-2008}
         * \param{sourceField the original field which will be distributed}
         * This function is used for dealing with the interpolation of field
         * values
         */       
        template<class Type> 
        tmp<Field<Type> > cloneField
        (
            const Field<Type>& sourceField
        ) const;


    /**
     * \brief{Used to build a rotational symmetric patch for the slave targetpatch}
     * \date{11-04-2009}
     * \param{targetPatch the original targetPatch}
     * \param{lastFieldAdress scalar which marks the end of the field}
     * This function is used to generate 360 degree for the target Patch.
     */
    template<class Patch>
    Patch& buildSlave360Patch(const Patch& targetPatch,scalar& lastFieldAdress);

        /**
         * \brief{Used to distribute the values other the interpoaltion disk .}
         * \date{31-08-2008}
         * \param{sourceField the original field which will be distributed}
         * This function is used for dealing with the interpolation of field
         * values
         */      
        template<class Type> 
        tmp<Field<Type> > cloneField
        (
           const tmp<Field<Type> >& sourceField
        ) const; 
        
        /**
         * \brief{Used to distribute the values other the interpoaltion disk .}
         * \date{31-08-2008}
         * \param{sourceField the original field which will be distributed}
         * \param{cmpt the direction it is necessary to rotate the scalar field correctly}
         * \param{rank of the original field necessary to rotate the scalar field correctly}
         * This function is used for dealing with the interpolation of field
         * values
         */              
        template<class Type> 
        tmp<Field<Type> > cloneField
        (
           const Field<Type>& sourceField,
           const direction& cmpt,
           const int& rank
        ) const;
        
        /**
         * \brief{Used to distribute the values other the interpoaltion disk .}
         * \date{31-08-2008}
         * \param{sourceField the original field which will be distributed}
         * \param{cmpt the direction it is necessary to rotate the scalar field correctly}
         * \param{rank of the original field necessary to rotate the scalar field correctly}
         * This function is used for dealing with the interpolation of field
         * values
         */                                
        template<class Type> 
        tmp<Field<Type> > cloneField
        (
           const tmp<Field<Type> >& sourceField,
           const direction& cmpt,
           const int& rank
           
        ) const;         
        
        
   
         
protected:

    // Protected static data

        //- Relative tolerance (for geometric matching). Is factor of
        //  maximum edge length per face.
        static const scalar matchTol_;
    
    
    // Protected Member Functions

        // The primitiveDomainScalingPatch geometry initialisation is called by polyBoundaryMesh
        friend class polyBoundaryMesh;

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry()
        {}

        //- Calculate the patch geometry
        virtual void calcGeometry()
        {}

        virtual void clearAddressing();
        
        void CalcRadialGeometry();
        
        
        //- Initialise the patches for moving points
        virtual void initMovePoints(const pointField&)
        {}

        //- Correct patches after moving points
        virtual void movePoints(const pointField& p);


public:

    //- Runtime type information
    // TODO Must be enabled when it is included as a dynamic library
    TypeName("primitiveDomainScalingPatch"); 


    // Constructors

     

        
        
        /**
         * \brief{Constructor with additional clipping points}
         * \param<facesh> {}
         * \param<points> {}
         * \param<OrigPatch> {the reference to the original patch}
         * \param<CCSCenter> {the origin of the cylindrical coordinate system}
         * \param<CCSaxis> {the axis for the cylindrical coordinate system}
         * \param<CCSdirection> {the direction for the cylindrical coordinate system} 
         * \param<omega> {the angular velocity of the virtual patch}
         * \param<name> {the name of the patch}
         * \param<switchBlockMesh> {This switch is used to disable the calculation of the geometry in case of being called from the blockMesh executable}
         */        
        primitiveDomainScalingPatch
        (       
                const List<face>& faces, 
                const pointField& points,
                const primitivePatch& OrigPatch,
                const point& CCSCenter,
                const vector& CCSaxis,
                const vector& CCSdirection,
                scalar omega,
                const word& name,
                bool switchBlockMesh
        );

        //- Construct as copy
        primitiveDomainScalingPatch(const primitiveDomainScalingPatch&);
        
        
        
        /**
         * \brief{Used to set the timestep and the actual runtime .}
         * \date{01-04-2008}
         * \param{timestep the actual timestep width}
         * \param{runTime reference to the actual time object}         
         * This function is used for dealing with the interpolation of field
         * values
         */                      
        void setDeltaTime(scalar timeStep, const Time& runTime);
       // Construct from dictionary
        
        /**
         * \brief{Used to set the timestep and the actual runtime .}
         * \date{01-04-2008}
         * \param{timestep the actual timestep width}
         * \param{runTime reference to the actual time object}         
         * This function is used for dealing with the interpolation of field
         * values
         */                              
        //TODO FB check if that is necessary
        const List<face>& GetFaces() const;        
        
        //TODO FB when is that called ??
        const pointField& GetFacePoints() const;
        
        /**
         * \brief{Used to set the timestep and the actual runtime .}
         * \date{01-04-2008}
         * \param{timestep the actual timestep width}
         * \param{runTime reference to the actual time object}         
         * This function is used for dealing with the interpolation of field
         * values
         */            
        template<class Type, class Patch> 
        tmp<Field<Type> > interpolateToNeighbour
        (
               const Field<Type>& sourceField,  
               const Patch& targetPatch
        );
        
        /**
         * \brief{Used to set the timestep and the actual runtime .}
         * \date{01-04-2008}
         * \param{timestep the actual timestep width}
         * \param{runTime reference to the actual time object}         
         * This function is used for dealing with the interpolation of field
         * values
         */            
        template<class Type, class Patch> 
        tmp<Field<Type> > interpolateToNeighbour
        (
             const tmp<Field<Type> >& sourceField,  
             const Patch& targetPatch
        );
        
        /**
         * \brief{Used to set the timestep and the actual runtime .}
         * \date{01-04-2008}
         * \param{timestep the actual timestep width}
         * \param{runTime reference to the actual time object}         
         * This function is used for dealing with the interpolation of field
         * values
         */            
        template<class Type, class Patch> 
        tmp<Field<Type> > interpolateToNeighbour
        (
               const Field<Type>& sourceField,  
               const Patch& targetPatch,
               const direction& cmpt,
               const int& rank
        );
        
        /**
         * \brief{Used to set the timestep and the actual runtime .}
         * \date{01-04-2008}
         * \param{timestep the actual timestep width}
         * \param{runTime reference to the actual time object}         
         * This function is used for dealing with the interpolation of field
         * values
         */            
        template<class Type, class Patch> 
        tmp<Field<Type> > interpolateToNeighbour
        (
             const tmp<Field<Type> >& sourceField,  
             const Patch& targetPatch,
             const direction& cmpt,
             const int& rank
        );

        
        
        // TODO FB when is that called
        void GetFacePoints(pointField& PF) const;
        // TODO FB when is that called
        void GetFacePoints(List<face>& Faces) const;
        
    // Selectors

        //- Return a point er to a new patch created on freestore from
        //  components
        static autoPtr<primitiveDomainScalingPatch> New
        (
            const List<face>& faces,
            const pointField& points
        );
        
        // Destructor
        virtual ~primitiveDomainScalingPatch();
        
        // Member functions
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
            
         void adoptToOuterTimeStep(scalar angle);
         
         // Access
         /**
          * \brief{Used to get the omega value .}
          * \date{01-04-2008}
          */               
         const scalar& omega()const
         {
             return (omega_);
         }
         
         /**
          * \brief{Used to delete the interpolation pointer from the outside.}
          * \date{01-04-2008}
          */         
         void deleteInterpolationPointers()
         {   
             deleteDemandDrivenData(masterToPatchPtr_);
         }
            
         /**
          * \brief{Used to set the start angle.}
          * \date{01-04-2008}
          */         
         void setStartAngle(scalar startAngle)
         {
             startAngle_ = startAngle;
             calcRadialGeometry();
         }
      
         
//         void setRelative(Switch relative)
//         {
//             relative_ = relative;
//         }
         
         void setMaster(Switch master)
         {
             master_ = master;
         }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#      include "primitiveDomainScalingPatchTemplates.C"
#      include "primitiveDomainScalingPatchCmptTemplates.C"
#endif

#endif

// ************************************************************************* //
