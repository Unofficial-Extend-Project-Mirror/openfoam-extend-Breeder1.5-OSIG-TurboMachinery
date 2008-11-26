// Parametrized test case for the ERCOFTAC diffuser.

//         Created by Omar Bounous 

//Run using:
//m4 -P blockMeshDict.m4 > blockMeshDict

//m4 definitions:
m4_changecom(//)m4_changequote([,])
m4_define(calc, [m4_esyscmd(perl -e 'use Math::Trig; printf ($1)')])
m4_define(VCOUNT, 0)
m4_define(vlabel, [[// ]Vertex $1 = VCOUNT m4_define($1, VCOUNT)m4_define([VCOUNT], m4_incr(VCOUNT))])

//Mathematical constants:
m4_define(pi, 3.1415926536)

//Geometry
m4_define(openingAngle, 10.0)
m4_define(wedgeAngle, 5.0)
m4_define(diffuserLength, 0.51)
m4_define(extensionLength, 0.59)
m4_define(rIn,0.13)
m4_define(rOut, calc(rIn+diffuserLength*tan(openingAngle*pi/180.0)))

//Grid points (integers!):
m4_define(rNumberOfCells, 25)
m4_define(xABnumberOfCells, 10)
m4_define(xBCnumberOfCells, 4)
m4_define(xCDnumberOfCells, 6)
m4_define(xDEnumberOfCells, 30)
m4_define(xEFnumberOfCells, 10)
m4_define(rGrading, 0.2)

//Plane A:
m4_define(xA, -0.50)
m4_define(rA, rIn)

//Plane B:
m4_define(xB, -0.10)
m4_define(rB, rIn)

//Plane C:
m4_define(xC, -0.025)
m4_define(rC, rIn)

//Plane D:
m4_define(xD, 0)
m4_define(rD, rIn)

//Plane E:
m4_define(xE, diffuserLength)
m4_define(rE, rOut)

//Plane F:
m4_define(xF, calc(diffuserLength+extensionLength))
m4_define(rF, rOut)

/*---------------------------------------------------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  1.5                                   |
|   \\  /    A nd           | Web:      http://www.openfoam.org               |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/

FoamFile
{
    version         2.0;
    format          ascii;

    root            "";
    case            "";
    instance        "";
    local           "";

    class           dictionary;
    object          blockMeshDict;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 1;

vertices
(
//Plane A:
(xA 0 0) vlabel(A0)
(xA calc(rA*cos((wedgeAngle/2)*pi/180.0)) -calc(rA*sin((wedgeAngle/2)*pi/180.0))) vlabel(A1)
(xA calc(rA*cos((wedgeAngle/2)*pi/180.0)) calc(rA*sin((wedgeAngle/2)*pi/180.0))) vlabel(A2)

//Plane B:
(xB 0 0) vlabel(B0)
(xB calc(rB*cos((wedgeAngle/2)*pi/180.0)) -calc(rB*sin((wedgeAngle/2)*pi/180.0))) vlabel(B1)
(xB calc(rB*cos((wedgeAngle/2)*pi/180.0)) calc(rB*sin((wedgeAngle/2)*pi/180.0))) vlabel(B2)

//Plane C:
(xC 0 0) vlabel(C0)
(xC calc(rC*cos((wedgeAngle/2)*pi/180.0)) -calc(rC*sin((wedgeAngle/2)*pi/180.0))) vlabel(C1)
(xC calc(rC*cos((wedgeAngle/2)*pi/180.0)) calc(rC*sin((wedgeAngle/2)*pi/180.0))) vlabel(C2)

//Plane D:
(xD 0 0) vlabel(D0)
(xD calc(rD*cos((wedgeAngle/2)*pi/180.0)) -calc(rD*sin((wedgeAngle/2)*pi/180.0))) vlabel(D1)
(xD calc(rD*cos((wedgeAngle/2)*pi/180.0)) calc(rD*sin((wedgeAngle/2)*pi/180.0))) vlabel(D2)

//Plane E:
(xE 0 0) vlabel(E0)
(xE calc(rE*cos((wedgeAngle/2)*pi/180.0)) -calc(rE*sin((wedgeAngle/2)*pi/180.0))) vlabel(E1)
(xE calc(rE*cos((wedgeAngle/2)*pi/180.0)) calc(rE*sin((wedgeAngle/2)*pi/180.0))) vlabel(E2)

//Plane F:
(xF 0 0) vlabel(F0)
(xF calc(rF*cos((wedgeAngle/2)*pi/180.0)) -calc(rF*sin((wedgeAngle/2)*pi/180.0))) vlabel(F1)
(xF calc(rF*cos((wedgeAngle/2)*pi/180.0)) calc(rF*sin((wedgeAngle/2)*pi/180.0))) vlabel(F2)
);

// Defining blocks:
blocks
(
    //Blocks between plane A and plane B:
    // block0 
    hex (A0 B0 B1 A1 A0 B0 B2 A2) AB
    (xABnumberOfCells rNumberOfCells 1) 
    simpleGrading (1 rGrading 1)

    //Blocks between plane B and plane C:
    // block0
    hex (B0 C0 C1 B1 B0 C0 C2 B2) BC
    (xBCnumberOfCells rNumberOfCells 1)
    simpleGrading (1 rGrading 1)
    
    //Blocks between plane C and plane D:
    // block0
    hex (C0 D0 D1 C1 C0 D0 D2 C2) CD
    (xCDnumberOfCells rNumberOfCells 1)
    simpleGrading (1 rGrading 1)
    
    //Blocks between plane D and plane E:
    // block0
    hex (D0 E0 E1 D1 D0 E0 E2 D2) DE
    (xDEnumberOfCells rNumberOfCells 1)
    simpleGrading (1 rGrading 1)
    
    //Blocks between plane E and plane F:
    // block0
    hex (E0 F0 F1 E1 E0 F0 F2 E2) EF
    (xEFnumberOfCells rNumberOfCells 1)
    simpleGrading (1 rGrading 1)
    
);

edges
(
    //Plane A:
    line A0 A1 
    arc A1 A2 (xA rA 0)
    line A2 A0

    //Plane B:
    line B0 B1
    arc B1 B2 (xB rB 0)
    line B2 B0

    //Plane C:
    line C0 C1
    arc C1 C2 (xC rC 0)
    line C2 C0

    //Plane D:
    line D0 D1
    arc D1 D2 (xD rD 0)
    line D2 D0

    //Plane E:
    line E0 E1
    arc E1 E2 (xE rE 0)
    line E2 E0

    //Plane F:
    line F0 F1
    arc F1 F2 (xF rF 0)
    line F2 F0

);

// Defining patches:
patches
(
    symmetryPlane axis
    (
        (A0 B0 B0 A0)
        (B0 C0 C0 B0)
        (C0 D0 D0 C0)
        (D0 E0 E0 D0)
        (E0 F0 F0 E0)
    )
    patch inlet
    (
       (A0 A2 A1 A0)
    )
    patch outlet
    (
       (F0 F1 F2 F0)
    )
    wall wallProlongation
    (
      (E1 E2 F2 F1)
    )
    wall wallDiffuser
    (
      (D1 D2 E2 E1)
    )
    wall statSwirlWall
    (
      (B1 B2 C2 C1)
      (C1 C2 D2 D1)
    )
    wall rotSwirlWall
    (
      (A1 A2 B2 B1)
    )
    wedge back
    (
      (A0 A1 B1 B0)
      (B0 B1 C1 C0)
      (C0 C1 D1 D0)
      (D0 D1 E1 E0)
      (E0 E1 F1 F0)
    )
    wedge front
    (
     (A2 A0 B0 B2)
     (B2 B0 C0 C2)
     (C2 C0 D0 D2)
     (D2 D0 E0 E2)
     (E2 E0 F0 F2)
     )
);

mergePatchPairs 
(
);

// ************************************************************************* //
