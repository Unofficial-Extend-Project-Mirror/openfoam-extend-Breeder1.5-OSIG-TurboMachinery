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
m4_define(rIn,0.13)
m4_define(rOut, calc(rIn+diffuserLength*tan(openingAngle*pi/180.0)))
m4_define(rDump, calc(3*rOut))
m4_define(dumpLength, calc(10*rOut))
m4_define(outletPipeLength, calc(6*rOut))

//Grid points (integers!):
m4_define(rNumberOfCells1st, 15)  // in the first O-grid
m4_define(rNumberOfCells2nd, 15)  // in the second O-grid
m4_define(rNumberOfCells3rd, 15)  // in the third O-grid
m4_define(tNumberOfCells, 20)
m4_define(zABnumberOfCells, 10)
m4_define(zBCnumberOfCells, 4)
m4_define(zCDnumberOfCells, 6)
m4_define(zDEnumberOfCells, 30)
m4_define(zEFnumberOfCells, 10)
m4_define(zFGnumberOfCells, 10)
m4_define(zGHnumberOfCells, 10)

m4_define(rGrading1, 5)
m4_define(rGrading2, 0.2)
m4_define(zGrading1, 5)
m4_define(zGrading2, 0.2)

//Plane A:
m4_define(zA, -0.50)
m4_define(rA, rIn)
m4_define(rRelA, 0.7)
m4_define(rRelAc, 0.8)

//Plane B:
m4_define(zB, -0.10)
m4_define(rB, rIn)
m4_define(rRelB, 0.7)
m4_define(rRelBc, 0.8)

//Plane C:
m4_define(zC, -0.025)
m4_define(rC, rIn)
m4_define(rRelC, 0.7)
m4_define(rRelCc, 0.8)

//Plane D:
m4_define(zD, 0)
m4_define(rD, rIn)
m4_define(rRelD, 0.7)
m4_define(rRelDc, 0.8)

//Plane E:
m4_define(zE, diffuserLength)
m4_define(rE, rOut)
m4_define(rHalfDumpE, calc(2*rOut))
m4_define(rDumpE, rDump)
m4_define(rRelE, 0.7)
m4_define(rRelEc, 0.8)

//Plane F:
m4_define(zF, calc(diffuserLength+dumpLength/2))
m4_define(rF, rOut)
m4_define(rHalfDumpF, calc(2*rOut))
m4_define(rDumpF, rDump)
m4_define(rRelF, 0.7)
m4_define(rRelFc, 0.8)

//Plane G:
m4_define(zG, calc(diffuserLength+dumpLength))
m4_define(rG, rOut)
m4_define(rHalfDumpG, calc(2*rOut))
m4_define(rDumpG, rDump)
m4_define(rRelG, 0.7)
m4_define(rRelGc, 0.8)

//Plane H:
m4_define(zH, calc(diffuserLength+dumpLength+outletPipeLength))
m4_define(rH, rOut)
m4_define(rRelH, 0.7)
m4_define(rRelHc, 0.8)

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

//Plane D:
(calc(rRelD*rD*cos(pi/4)) -calc(rRelD*rD*sin(pi/4)) zD) vlabel(D0)
(calc(rRelD*rD*cos(pi/4)) calc(rRelD*rD*sin(pi/4)) zD) vlabel(D1)
(calc(-rRelD*rD*cos(pi/4)) calc(rRelD*rD*sin(pi/4)) zD) vlabel(D2)
(calc(-rRelD*rD*cos(pi/4)) -calc(rRelD*rD*sin(pi/4)) zD) vlabel(D3)
(calc(rD*cos(pi/4)) -calc(rD*sin(pi/4)) zD) vlabel(D4)
(calc(rD*cos(pi/4)) calc(rD*sin(pi/4)) zD) vlabel(D5)
(calc(-rD*cos(pi/4)) calc(rD*sin(pi/4)) zD) vlabel(D6)
(calc(-rD*cos(pi/4)) -calc(rD*sin(pi/4)) zD) vlabel(D7)

//Plane E:
(calc(rRelE*rE*cos(pi/4)) -calc(rRelE*rE*sin(pi/4)) zE) vlabel(E0)
(calc(rRelE*rE*cos(pi/4)) calc(rRelE*rE*sin(pi/4)) zE) vlabel(E1)
(calc(-rRelE*rE*cos(pi/4)) calc(rRelE*rE*sin(pi/4)) zE) vlabel(E2)
(calc(-rRelE*rE*cos(pi/4)) -calc(rRelE*rE*sin(pi/4)) zE) vlabel(E3)
(calc(rE*cos(pi/4)) -calc(rE*sin(pi/4)) zE) vlabel(E4)
(calc(rE*cos(pi/4)) calc(rE*sin(pi/4)) zE) vlabel(E5)
(calc(-rE*cos(pi/4)) calc(rE*sin(pi/4)) zE) vlabel(E6)
(calc(-rE*cos(pi/4)) -calc(rE*sin(pi/4)) zE) vlabel(E7)
(calc(rHalfDumpE*cos(pi/4)) -calc(rHalfDumpE*sin(pi/4)) zE) vlabel(E8)
(calc(rHalfDumpE*cos(pi/4)) calc(rHalfDumpE*sin(pi/4)) zE) vlabel(E9)
(calc(-rHalfDumpE*cos(pi/4)) calc(rHalfDumpE*sin(pi/4)) zE) vlabel(E10)
(calc(-rHalfDumpE*cos(pi/4)) -calc(rHalfDumpE*sin(pi/4)) zE) vlabel(E11)
(calc(rDumpE*cos(pi/4)) -calc(rDumpE*sin(pi/4)) zE) vlabel(E12)
(calc(rDumpE*cos(pi/4)) calc(rDumpE*sin(pi/4)) zE) vlabel(E13)
(calc(-rDumpE*cos(pi/4)) calc(rDumpE*sin(pi/4)) zE) vlabel(E14)
(calc(-rDumpE*cos(pi/4)) -calc(rDumpE*sin(pi/4)) zE) vlabel(E15)

//Plane F:
(calc(rRelF*rF*cos(pi/4)) -calc(rRelF*rF*sin(pi/4)) zF) vlabel(F0)
(calc(rRelF*rF*cos(pi/4)) calc(rRelF*rF*sin(pi/4)) zF) vlabel(F1)
(calc(-rRelF*rF*cos(pi/4)) calc(rRelF*rF*sin(pi/4)) zF) vlabel(F2)
(calc(-rRelF*rF*cos(pi/4)) -calc(rRelF*rF*sin(pi/4)) zF) vlabel(F3)
(calc(rF*cos(pi/4)) -calc(rF*sin(pi/4)) zF) vlabel(F4)
(calc(rF*cos(pi/4)) calc(rF*sin(pi/4)) zF) vlabel(F5)
(calc(-rF*cos(pi/4)) calc(rF*sin(pi/4)) zF) vlabel(F6)
(calc(-rF*cos(pi/4)) -calc(rF*sin(pi/4)) zF) vlabel(F7)
(calc(rHalfDumpF*cos(pi/4)) -calc(rHalfDumpF*sin(pi/4)) zF) vlabel(F8)
(calc(rHalfDumpF*cos(pi/4)) calc(rHalfDumpF*sin(pi/4)) zF) vlabel(F9)
(calc(-rHalfDumpF*cos(pi/4)) calc(rHalfDumpF*sin(pi/4)) zF) vlabel(F10)
(calc(-rHalfDumpF*cos(pi/4)) -calc(rHalfDumpF*sin(pi/4)) zF) vlabel(F11)
(calc(rDumpF*cos(pi/4)) -calc(rDumpF*sin(pi/4)) zF) vlabel(F12)
(calc(rDumpF*cos(pi/4)) calc(rDumpF*sin(pi/4)) zF) vlabel(F13)
(calc(-rDumpF*cos(pi/4)) calc(rDumpF*sin(pi/4)) zF) vlabel(F14)
(calc(-rDumpF*cos(pi/4)) -calc(rDumpF*sin(pi/4)) zF) vlabel(F15)

//Plane G:
(calc(rRelG*rG*cos(pi/4)) -calc(rRelG*rG*sin(pi/4)) zG) vlabel(G0)
(calc(rRelG*rG*cos(pi/4)) calc(rRelG*rG*sin(pi/4)) zG) vlabel(G1)
(calc(-rRelG*rG*cos(pi/4)) calc(rRelG*rG*sin(pi/4)) zG) vlabel(G2)
(calc(-rRelG*rG*cos(pi/4)) -calc(rRelG*rG*sin(pi/4)) zG) vlabel(G3)
(calc(rG*cos(pi/4)) -calc(rG*sin(pi/4)) zG) vlabel(G4)
(calc(rG*cos(pi/4)) calc(rG*sin(pi/4)) zG) vlabel(G5)
(calc(-rG*cos(pi/4)) calc(rG*sin(pi/4)) zG) vlabel(G6)
(calc(-rG*cos(pi/4)) -calc(rG*sin(pi/4)) zG) vlabel(G7)
(calc(rHalfDumpG*cos(pi/4)) -calc(rHalfDumpG*sin(pi/4)) zG) vlabel(G8)
(calc(rHalfDumpG*cos(pi/4)) calc(rHalfDumpG*sin(pi/4)) zG) vlabel(G9)
(calc(-rHalfDumpG*cos(pi/4)) calc(rHalfDumpG*sin(pi/4)) zG) vlabel(G10)
(calc(-rHalfDumpG*cos(pi/4)) -calc(rHalfDumpG*sin(pi/4)) zG) vlabel(G11)
(calc(rDumpG*cos(pi/4)) -calc(rDumpG*sin(pi/4)) zG) vlabel(G12)
(calc(rDumpG*cos(pi/4)) calc(rDumpG*sin(pi/4)) zG) vlabel(G13)
(calc(-rDumpG*cos(pi/4)) calc(rDumpG*sin(pi/4)) zG) vlabel(G14)
(calc(-rDumpG*cos(pi/4)) -calc(rDumpG*sin(pi/4)) zG) vlabel(G15)

//Plane H:
(calc(rRelH*rH*cos(pi/4)) -calc(rRelH*rH*sin(pi/4)) zH) vlabel(H0)
(calc(rRelH*rH*cos(pi/4)) calc(rRelH*rH*sin(pi/4)) zH) vlabel(H1)
(calc(-rRelH*rH*cos(pi/4)) calc(rRelH*rH*sin(pi/4)) zH) vlabel(H2)
(calc(-rRelH*rH*cos(pi/4)) -calc(rRelH*rH*sin(pi/4)) zH) vlabel(H3)
(calc(rH*cos(pi/4)) -calc(rH*sin(pi/4)) zH) vlabel(H4)
(calc(rH*cos(pi/4)) calc(rH*sin(pi/4)) zH) vlabel(H5)
(calc(-rH*cos(pi/4)) calc(rH*sin(pi/4)) zH) vlabel(H6)
(calc(-rH*cos(pi/4)) -calc(rH*sin(pi/4)) zH) vlabel(H7)
);

// Defining blocks:
blocks
(
    //Blocks between plane A and plane B:
    // block0 - positive x O-grid block
    hex (A5 A1 A0 A4 B5 B1 B0 B4 ) AB
    (rNumberOfCells1st tNumberOfCells zABnumberOfCells)
    simpleGrading (rGrading1 1 1)
    // block1 - positive y O-grid block
    hex (A6 A2 A1 A5 B6 B2 B1 B5 ) AB
    (rNumberOfCells1st tNumberOfCells zABnumberOfCells)
    simpleGrading (rGrading1 1 1)
    // block2 - negative x O-grid block
    hex (A7 A3 A2 A6 B7 B3 B2 B6 ) AB
    (rNumberOfCells1st tNumberOfCells zABnumberOfCells)
    simpleGrading (rGrading1 1 1)
    // block3 - negative y O-grid block
    hex (A4 A0 A3 A7 B4 B0 B3 B7 ) AB
    (rNumberOfCells1st tNumberOfCells zABnumberOfCells)
    simpleGrading (rGrading1 1 1)
    // block4 - central O-grid block
    hex (A0 A1 A2 A3 B0 B1 B2 B3 ) AB
    (tNumberOfCells tNumberOfCells zABnumberOfCells)
    simpleGrading (1 1 1)

    //Blocks between plane B and plane C:
    // block0 - positive x O-grid block
    hex (B5 B1 B0 B4 C5 C1 C0 C4 ) BC
    (rNumberOfCells1st tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading1 1 1)
    // block1 - positive y O-grid block
    hex (B6 B2 B1 B5 C6 C2 C1 C5 ) BC
    (rNumberOfCells1st tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading1 1 1)
    // block2 - negative x O-grid block
    hex (B7 B3 B2 B6 C7 C3 C2 C6 ) BC
    (rNumberOfCells1st tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading1 1 1)
    // block3 - negative y O-grid block
    hex (B4 B0 B3 B7 C4 C0 C3 C7 ) BC
    (rNumberOfCells1st tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading1 1 1)
    // block4 - central O-grid block
    hex (B0 B1 B2 B3 C0 C1 C2 C3 ) BC
    (tNumberOfCells tNumberOfCells zBCnumberOfCells)
    simpleGrading (1 1 1)

    //Blocks between plane C and plane D:
    // block0 - positive x O-grid block
    hex (C5 C1 C0 C4 D5 D1 D0 D4 ) CD
    (rNumberOfCells1st tNumberOfCells zCDnumberOfCells)
    simpleGrading (rGrading1 1 1)
    // block1 - positive y O-grid block
    hex (C6 C2 C1 C5 D6 D2 D1 D5 ) CD
    (rNumberOfCells1st tNumberOfCells zCDnumberOfCells)
    simpleGrading (rGrading1 1 1)
    // block2 - negative x O-grid block
    hex (C7 C3 C2 C6 D7 D3 D2 D6 ) CD
    (rNumberOfCells1st tNumberOfCells zCDnumberOfCells)
    simpleGrading (rGrading1 1 1)
    // block3 - negative y O-grid block
    hex (C4 C0 C3 C7 D4 D0 D3 D7 ) CD
    (rNumberOfCells1st tNumberOfCells zCDnumberOfCells)
    simpleGrading (rGrading1 1 1)
    // block4 - central O-grid block
    hex (C0 C1 C2 C3 D0 D1 D2 D3 ) CD
    (tNumberOfCells tNumberOfCells zCDnumberOfCells)
    simpleGrading (1 1 1)

    //Blocks between plane D and plane E:
    // block0 - positive x O-grid block
    hex (D5 D1 D0 D4 E5 E1 E0 E4 ) DE
    (rNumberOfCells1st tNumberOfCells zDEnumberOfCells)
    simpleGrading (rGrading1 1 1)
    // block1 - positive y O-grid block
    hex (D6 D2 D1 D5 E6 E2 E1 E5 ) DE
    (rNumberOfCells1st tNumberOfCells zDEnumberOfCells)
    simpleGrading (rGrading1 1 1)
    // block2 - negative x O-grid block
    hex (D7 D3 D2 D6 E7 E3 E2 E6 ) DE
    (rNumberOfCells1st tNumberOfCells zDEnumberOfCells)
    simpleGrading (rGrading1 1 1)
    // block3 - negative y O-grid block
    hex (D4 D0 D3 D7 E4 E0 E3 E7 ) DE
    (rNumberOfCells1st tNumberOfCells zDEnumberOfCells)
    simpleGrading (rGrading1 1 1)
    // block4 - central O-grid block
    hex (D0 D1 D2 D3 E0 E1 E2 E3 ) DE
    (tNumberOfCells tNumberOfCells zDEnumberOfCells)
    simpleGrading (1 1 1)

    //Blocks between plane E and plane F:
    // block0 - positive x O-grid block 1st belt
    hex (E5 E1 E0 E4 F5 F1 F0 F4 ) EF
    (rNumberOfCells1st tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading1 1 zGrading1)
    // block1 - positive y O-grid block 1st belt
    hex (E6 E2 E1 E5 F6 F2 F1 F5 ) EF
    (rNumberOfCells1st tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading1 1 zGrading1)
    // block2 - negative x O-grid block 1st belt
    hex (E7 E3 E2 E6 F7 F3 F2 F6 ) EF
    (rNumberOfCells1st tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading1 1 zGrading1)
    // block3 - negative y O-grid block 1st belt
    hex (E4 E0 E3 E7 F4 F0 F3 F7 ) EF
    (rNumberOfCells1st tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading1 1 zGrading1)
    // block4 - central O-grid block 
    hex (E0 E1 E2 E3 F0 F1 F2 F3 ) EF
    (tNumberOfCells tNumberOfCells zEFnumberOfCells)
    simpleGrading (1 1 zGrading1)
    // block5 - positive x O-grid block 2nd belt 
    hex (E9 E5 E4 E8 F9 F5 F4 F8 ) EF
    (rNumberOfCells2nd tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading2 1 zGrading1)
    // block6 - positive y O-grid block 2nd belt
    hex (E10 E6 E5 E9 F10 F6 F5 F9 ) EF
    (rNumberOfCells2nd tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading2 1 zGrading1)
    // block7 - negative x O-grid block 2nd belt
    hex (E11 E7 E6 E10 F11 F7 F6 F10 ) EF
    (rNumberOfCells2nd tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading2 1 zGrading1)
    // block8 - negative y O-grid block 2nd belt
    hex (E8 E4 E7 E11 F8 F4 F7 F11 ) EF
    (rNumberOfCells2nd tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading2 1 zGrading1)
    // block9 - positive x O-grid block 3rd belt
    hex (E13 E9 E8 E12 F13 F9 F8 F12 ) EF
    (rNumberOfCells3rd tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading1 1 zGrading1)
    // block10 - positive y O-grid block 3rd belt
    hex (E14 E10 E9 E13 F14 F10 F9 F13 ) EF
    (rNumberOfCells3rd tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading1 1 zGrading1)
    // block11 - negative x O-grid block 3rd belt
    hex (E15 E11 E10 E14 F15 F11 F10 F14 ) EF
    (rNumberOfCells3rd tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading1 1 zGrading1)
    // block12 - negative y O-grid block 3rd belt
    hex (E12 E8 E11 E15 F12 F8 F11 F15 ) EF
    (rNumberOfCells3rd tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading1 1 zGrading1)

//Blocks between plane F and plane G:
    // block0 - positive x O-grid block 1st belt
    hex (F5 F1 F0 F4 G5 G1 G0 G4 ) FG
    (rNumberOfCells1st tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading1 1 zGrading2)
    // block1 - positive y O-grid block 1st belt
    hex (F6 F2 F1 F5 G6 G2 G1 G5 ) FG
    (rNumberOfCells1st tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading1 1 zGrading2)
    // block2 - negative x O-grid block 1st belt
    hex (F7 F3 F2 F6 G7 G3 G2 G6 ) FG
    (rNumberOfCells1st tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading1 1 zGrading2)
    // block3 - negative y O-grid block 1st belt
    hex (F4 F0 F3 F7 G4 G0 G3 G7 ) FG
    (rNumberOfCells1st tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading1 1 zGrading2)
    // block4 - central O-grid block 
    hex (F0 F1 F2 F3 G0 G1 G2 G3 ) FG
    (tNumberOfCells tNumberOfCells zFGnumberOfCells)
    simpleGrading (1 1 zGrading2)
    // block5 - positive x O-grid block 2nd belt 
    hex (F9 F5 F4 F8 G9 G5 G4 G8 ) FG
    (rNumberOfCells2nd tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading2 1 zGrading2)
    // block6 - positive y O-grid block 2nd belt
    hex (F10 F6 F5 F9 G10 G6 G5 G9 ) FG
    (rNumberOfCells2nd tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading2 1 zGrading2)
    // block7 - negative x O-grid block 2nd belt
    hex (F11 F7 F6 F10 G11 G7 G6 G10 ) FG
    (rNumberOfCells2nd tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading2 1 zGrading2)
    // block8 - negative y O-grid block 2nd belt
    hex (F8 F4 F7 F11 G8 G4 G7 G11 ) FG
    (rNumberOfCells2nd tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading2 1 zGrading2)
    // block9 - positive x O-grid block 3rd belt
    hex (F13 F9 F8 F12 G13 G9 G8 G12 ) FG
    (rNumberOfCells3rd tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading1 1 zGrading2)
    // block10 - positive y O-grid block 3rd belt
    hex (F14 F10 F9 F13 G14 G10 G9 G13 ) FG
    (rNumberOfCells3rd tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading1 1 zGrading2)
    // block11 - negative x O-grid block 3rd belt
    hex (F15 F11 F10 F14 G15 G11 G10 G14 ) FG
    (rNumberOfCells3rd tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading1 1 zGrading2)
    // block12 - negative y O-grid block 3rd belt
    hex (F12 F8 F11 F15 G12 G8 G11 G15 ) FG
    (rNumberOfCells3rd tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading1 1 zGrading2)

//Blocks between plane G and plane H:
    // block0 - positive x O-grid block 1st belt
    hex (G5 G1 G0 G4 H5 H1 H0 H4 ) GH
    (rNumberOfCells1st tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading1 1 zGrading1)
    // block1 - positive y O-grid block 1st belt
    hex (G6 G2 G1 G5 H6 H2 H1 H5 ) GH
    (rNumberOfCells1st tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading1 1 zGrading1)
    // block2 - negative x O-grid block 1st belt
    hex (G7 G3 G2 G6 H7 H3 H2 H6 ) GH
    (rNumberOfCells1st tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading1 1 zGrading1)
    // block3 - negative y O-grid block 1st belt
    hex (G4 G0 G3 G7 H4 H0 H3 H7 ) GH
    (rNumberOfCells1st tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading1 1 zGrading1)
    // block4 - central O-grid block 
    hex (G0 G1 G2 G3 H0 H1 H2 H3 ) GH
    (tNumberOfCells tNumberOfCells zFGnumberOfCells)
    simpleGrading (1 1 zGrading1)
   
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

    //Plane D:
    arc D0 D1 (calc(rRelDc*rRelD*rD) 0 zD)
    arc D1 D2 (0 calc(rRelDc*rRelD*rD) zD)
    arc D2 D3 (-calc(rRelDc*rRelD*rD) 0 zD)
    arc D3 D0 (0 -calc(rRelDc*rRelD*rD) zD)
    arc D4 D5 (rD 0 zD)
    arc D5 D6 (0 rD zD)
    arc D6 D7 (-rD 0 zD)
    arc D7 D4 (0 -rD zD)

    //Plane E:
    arc E0  E1  (calc(rRelEc*rRelE*rE) 0 zE)
    arc E1  E2  (0 calc(rRelEc*rRelE*rE) zE)
    arc E2  E3  (-calc(rRelEc*rRelE*rE) 0 zE)
    arc E3  E0  (0 -calc(rRelEc*rRelE*rE) zE)
    arc E4  E5  (rE 0 zE)
    arc E5  E6  (0 rE zE)
    arc E6  E7  (-rE 0 zE)
    arc E7  E4  (0 -rE zE)
    arc E8  E9  (rHalfDumpE 0 zE)
    arc E9  E10 (0 rHalfDumpE zE)
    arc E10 E11 (-rHalfDumpE 0 zE)
    arc E11 E8  (0 -rHalfDumpE zE)
    arc E12 E13 (rDumpE 0 zE)
    arc E13 E14 (0 rDumpE zE)
    arc E14 E15 (-rDumpE 0 zE)
    arc E15 E12 (0 -rDumpE zE)

    //Plane F:
    arc F0  F1 (calc(rRelFc*rRelF*rF) 0 zF)
    arc F1  F2 (0 calc(rRelFc*rRelF*rF) zF)
    arc F2  F3 (-calc(rRelFc*rRelF*rF) 0 zF)
    arc F3  F0 (0 -calc(rRelFc*rRelF*rF) zF)
    arc F4  F5 (rF 0 zF)
    arc F5  F6 (0 rF zF)
    arc F6  F7 (-rF 0 zF)
    arc F7  F4 (0 -rF zF)
    arc F8  F9  (rHalfDumpF 0 zF)
    arc F9  F10 (0 rHalfDumpF zF)
    arc F10 F11 (-rHalfDumpF 0 zF)
    arc F11 F8  (0 -rHalfDumpF zF)
    arc F12 F13 (rDumpF 0 zF)
    arc F13 F14 (0 rDumpF zF)
    arc F14 F15 (-rDumpF 0 zF)
    arc F15 F12 (0 -rDumpF zF)

    //Plane G:
    arc G0  G1 (calc(rRelGc*rRelG*rG) 0 zG)
    arc G1  G2 (0 calc(rRelGc*rRelG*rG) zG)
    arc G2  G3 (-calc(rRelGc*rRelG*rG) 0 zG)
    arc G3  G0 (0 -calc(rRelGc*rRelG*rG) zG)
    arc G4  G5 (rG 0 zG)
    arc G5  G6 (0 rG zG)
    arc G6  G7 (-rG 0 zG)
    arc G7  G4 (0 -rG zG)
    arc G8  G9  (rHalfDumpG 0 zG)
    arc G9  G10 (0 rHalfDumpG zG)
    arc G10 G11 (-rHalfDumpG 0 zG)
    arc G11 G8  (0 -rHalfDumpG zG)
    arc G12 G13 (rDumpG 0 zG)
    arc G13 G14 (0 rDumpG zG)
    arc G14 G15 (-rDumpG 0 zG)
    arc G15 G12 (0 -rDumpG zG)

    //Plane H:
    arc H0  H1  (calc(rRelHc*rRelH*rH) 0 zH)
    arc H1  H2  (0 calc(rRelHc*rRelH*rH) zH)
    arc H2  H3  (-calc(rRelHc*rRelH*rH) 0 zH)
    arc H3  H0  (0 -calc(rRelHc*rRelH*rH) zH)
    arc H4  H5  (rH 0 zH)
    arc H5  H6  (0 rH zH)
    arc H6  H7  (-rH 0 zH)
    arc H7  H4  (0 -rH zH)
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
       (H0 H4 H5 H1)
       (H1 H5 H6 H2)
       (H2 H6 H7 H3)
       (H3 H7 H4 H0)
       (H0 H1 H2 H3)
    )
    wall wallDump
    (	
	(E5 E9 E8 E4)
	(E6 E10 E9 E5)
	(E7 E11 E10 E6)
	(E4 E8 E11 E7)
	(E9 E13 E12 E8)
	(E10 E14 E13 E9)
	(E11 E15 E14 E10)
	(E8 E12 E15 E11)
	(E12 E13 F13 F12)
	(E13 E14 F14 F13)
	(E14 E15 F15 F14)
	(E15 E12 F12 F15)
	(F12 F13 G13 G12)
	(F13 F14 G14 G13)
	(F14 F15 G15 G14)
	(F15 F12 G12 G15)
	(G4 G8 G9 G5)
	(G5 G9 G10 G6)
	(G6 G10 G11 G7)
	(G7 G11 G8 G4)
	(G8 G12 G13 G9)
	(G9 G13 G14 G10)
	(G10 G14 G15 G11)
	(G11 G15 G12 G8)
    ) 
wall wallOutletPipe
    (
       (G4 G5 H5 H4)
       (G5 G6 H6 H5)
       (G6 G7 H7 H6)
       (G7 G4 H4 H7)
    )
    wall wallDiffuser
    (
       (D4 D5 E5 E4)
       (D5 D6 E6 E5)
       (D6 D7 E7 E6)
       (D7 D4 E4 E7)
    )
    wall statSwirlWall
    (
       (B4 B5 C5 C4)
       (B5 B6 C6 C5)
       (B6 B7 C7 C6)
       (B7 B4 C4 C7)

       (C4 C5 D5 D4)
       (C5 C6 D6 D5)
       (C6 C7 D7 D6)
       (C7 C4 D4 D7)
    )
    wall rotSwirlWall
    (
       (A4 A5 B5 B4)
       (A5 A6 B6 B5)
       (A6 A7 B7 B6)
       (A7 A4 B4 B7)
    )
);

mergePatchPairs 
(
);

// ************************************************************************* //
