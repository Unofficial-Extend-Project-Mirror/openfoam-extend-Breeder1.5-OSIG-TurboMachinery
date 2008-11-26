// Parametrized test case for the ERCOFTAC diffusor.

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
m4_define(diffuserLength, 0.51)
m4_define(rIn, 0.13)
m4_define(rOut, calc(rIn+diffuserLength*tan(openingAngle*pi/180.0)))

//Grid points (integers!):
m4_define(rNumberOfCells, 15)
m4_define(tNumberOfCells, 20)
m4_define(zABnumberOfCells, 6)
m4_define(zBCnumberOfCells, 30)
m4_define(rGrading, 5)

//Plane A:
m4_define(zA, -0.025)
m4_define(rA, rIn)
m4_define(rRelA, 0.7)
m4_define(rRelAc, 0.8)

//Plane B:
m4_define(zB, 0)
m4_define(rB, rIn)
m4_define(rRelB, 0.7)
m4_define(rRelBc, 0.8)

//Plane C:
m4_define(zC, diffuserLength)
m4_define(rC, rOut)
m4_define(rRelC, 0.7)
m4_define(rRelCc, 0.8)

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
(calc(rRelA*rA*cos(pi/4)) -calc(rRelA*rA*sin(pi/4)) zA) vlabel(A0)
(calc(rRelA*rA*cos(pi/4)) calc(rRelA*rA*sin(pi/4)) zA) vlabel(A1)
(calc(-rRelA*rA*cos(pi/4)) calc(rRelA*rA*sin(pi/4)) zA) vlabel(A2)
(calc(-rRelA*rA*cos(pi/4)) -calc(rRelA*rA*sin(pi/4)) zA) vlabel(A3)
(calc(rA*cos(pi/4)) -calc(rA*sin(pi/4)) zA) vlabel(A4)
(calc(rA*cos(pi/4)) calc(rA*sin(pi/4)) zA) vlabel(A5)
(calc(-rA*cos(pi/4)) calc(rA*sin(pi/4)) zA) vlabel(A6)
(calc(-rA*cos(pi/4)) -calc(rA*sin(pi/4)) zA) vlabel(A7)

//Plane B:
(calc(rRelB*rB*cos(pi/4)) -calc(rRelB*rB*sin(pi/4)) zB) vlabel(B0)
(calc(rRelB*rB*cos(pi/4)) calc(rRelB*rB*sin(pi/4)) zB) vlabel(B1)
(calc(-rRelB*rB*cos(pi/4)) calc(rRelB*rB*sin(pi/4)) zB) vlabel(B2)
(calc(-rRelB*rB*cos(pi/4)) -calc(rRelB*rB*sin(pi/4)) zB) vlabel(B3)
(calc(rB*cos(pi/4)) -calc(rB*sin(pi/4)) zB) vlabel(B4)
(calc(rB*cos(pi/4)) calc(rB*sin(pi/4)) zB) vlabel(B5)
(calc(-rB*cos(pi/4)) calc(rB*sin(pi/4)) zB) vlabel(B6)
(calc(-rB*cos(pi/4)) -calc(rB*sin(pi/4)) zB) vlabel(B7)

//Plane C:
(calc(rRelC*rC*cos(pi/4)) -calc(rRelC*rC*sin(pi/4)) zC) vlabel(C0)
(calc(rRelC*rC*cos(pi/4)) calc(rRelC*rC*sin(pi/4)) zC) vlabel(C1)
(calc(-rRelC*rC*cos(pi/4)) calc(rRelC*rC*sin(pi/4)) zC) vlabel(C2)
(calc(-rRelC*rC*cos(pi/4)) -calc(rRelC*rC*sin(pi/4)) zC) vlabel(C3)
(calc(rC*cos(pi/4)) -calc(rC*sin(pi/4)) zC) vlabel(C4)
(calc(rC*cos(pi/4)) calc(rC*sin(pi/4)) zC) vlabel(C5)
(calc(-rC*cos(pi/4)) calc(rC*sin(pi/4)) zC) vlabel(C6)
(calc(-rC*cos(pi/4)) -calc(rC*sin(pi/4)) zC) vlabel(C7)
);

// Defining blocks:
blocks
(

    //Blocks between plane A and plane B:
    // block0 - positive x O-grid block
    hex (A5 A1 A0 A4 B5 B1 B0 B4 ) AB
    (rNumberOfCells tNumberOfCells zABnumberOfCells)
    simpleGrading (rGrading 1 1)
    // block1 - positive y O-grid block
    hex (A6 A2 A1 A5 B6 B2 B1 B5 ) AB
    (rNumberOfCells tNumberOfCells zABnumberOfCells)
    simpleGrading (rGrading 1 1)
    // block2 - negative x O-grid block
    hex (A7 A3 A2 A6 B7 B3 B2 B6 ) AB
    (rNumberOfCells tNumberOfCells zABnumberOfCells)
    simpleGrading (rGrading 1 1)
    // block3 - negative y O-grid block
    hex (A4 A0 A3 A7 B4 B0 B3 B7 ) AB
    (rNumberOfCells tNumberOfCells zABnumberOfCells)
    simpleGrading (rGrading 1 1)
    // block4 - central O-grid block
    hex (A0 A1 A2 A3 B0 B1 B2 B3 ) AB
    (tNumberOfCells tNumberOfCells zABnumberOfCells)
    simpleGrading (1 1 1)

    //Blocks between plane B and plane C:
    // block0 - positive x O-grid block
    hex (B5 B1 B0 B4 C5 C1 C0 C4 ) BC
    (rNumberOfCells tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading 1 1)
    // block1 - positive y O-grid block
    hex (B6 B2 B1 B5 C6 C2 C1 C5 ) BC
    (rNumberOfCells tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading 1 1)
    // block2 - negative x O-grid block
    hex (B7 B3 B2 B6 C7 C3 C2 C6 ) BC
    (rNumberOfCells tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading 1 1)
    // block3 - negative y O-grid block
    hex (B4 B0 B3 B7 C4 C0 C3 C7 ) BC
    (rNumberOfCells tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading 1 1)
    // block4 - central O-grid block
    hex (B0 B1 B2 B3 C0 C1 C2 C3 ) BC
    (tNumberOfCells tNumberOfCells zBCnumberOfCells)
    simpleGrading (1 1 1)

);

edges
(
    //Plane A:
    arc A0 A1 (calc(rRelAc*rRelA*rA) 0 zA)
    arc A1 A2 (0 calc(rRelAc*rRelA*rA) zA)
    arc A2 A3 (-calc(rRelAc*rRelA*rA) 0 zA)
    arc A3 A0 (0 -calc(rRelAc*rRelA*rA) zA)
    arc A4 A5 (rA 0 zA)
    arc A5 A6 (0 rA zA)
    arc A6 A7 (-rA 0 zA)
    arc A7 A4 (0 -rA zA)

    //Plane B:
    arc B0 B1 (calc(rRelBc*rRelB*rB) 0 zB)
    arc B1 B2 (0 calc(rRelBc*rRelB*rB) zB)
    arc B2 B3 (-calc(rRelBc*rRelB*rB) 0 zB)
    arc B3 B0 (0 -calc(rRelBc*rRelB*rB) zB)
    arc B4 B5 (rB 0 zB)
    arc B5 B6 (0 rB zB)
    arc B6 B7 (-rB 0 zB)
    arc B7 B4 (0 -rB zB)

    //Plane C:
    arc C0 C1 (calc(rRelCc*rRelC*rC) 0 zC)
    arc C1 C2 (0 calc(rRelCc*rRelC*rC) zC)
    arc C2 C3 (-calc(rRelCc*rRelC*rC) 0 zC)
    arc C3 C0 (0 -calc(rRelCc*rRelC*rC) zC)
    arc C4 C5 (rC 0 zC)
    arc C5 C6 (0 rC zC)
    arc C6 C7 (-rC 0 zC)
    arc C7 C4 (0 -rC zC)

);

// Defining patches:
patches
(
    patch inlet
    (
       (A1 A5 A4 A0)
       (A2 A6 A5 A1)
       (A3 A7 A6 A2)
       (A0 A4 A7 A3)
       (A3 A2 A1 A0)
    )
    patch outlet
    (
       (C0 C4 C5 C1)
       (C1 C5 C6 C2)
       (C2 C6 C7 C3)
       (C3 C7 C4 C0)
       (C0 C1 C2 C3)
    )
    wall wallDiffuser
    (
       (A4 A5 B5 B4)
       (A5 A6 B6 B5)
       (A6 A7 B7 B6)
       (A7 A4 B4 B7)

       (B4 B5 C5 C4)
       (B5 B6 C6 C5)
       (B6 B7 C7 C6)
       (B7 B4 C4 C7)
    )
);

mergePatchPairs 
(
);

// ************************************************************************* //
