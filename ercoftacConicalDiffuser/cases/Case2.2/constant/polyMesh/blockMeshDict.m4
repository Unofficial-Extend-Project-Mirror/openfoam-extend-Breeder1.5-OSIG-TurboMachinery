// Parametrized test case for the ERCOFTAC diffusor.

// Created by Omar Bounous



//Run using:
//m4 -P blockMeshDict.m4 > blockMeshDict

//m4 definitions:
m4_changecom(//)m4_changequote([,])
m4_define(calc, [m4_esyscmd(perl -e 'use Math::Trig; printf ($1)')])
m4_define(VCOUNT, 0)
m4_define(vlabel, [[// ]Vertex $1 = VCOUNT m4_define($1, VCOUNT)m4_define([VCOUNT], m4_incr(VCOUNT))])

//Mathematical constants:
m4_define(pi, 3.1415926536)
m4_define(sqr,0.707106781186548)


//Geometry
m4_define(openingAngle, 10.0)
m4_define(diffuserLength, 0.51)
m4_define(rIn,0.13)
m4_define(rOut, calc(rIn+diffuserLength*tan(openingAngle*pi/180.0)))
m4_define(rDump, calc(3*rOut))
m4_define(dumpLength, calc(10*rOut))
m4_define(outletPipeLength, calc(6*rOut))

//Grid points (integers!):
m4_define(rNumberOfCells1st, 25)  // in the first O-grid
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

//Plane B:
m4_define(zB, -0.10)
m4_define(rB, rIn)

//Plane C:
m4_define(zC, -0.025)
m4_define(rC, rIn)

//Plane D:
m4_define(zD, 0)
m4_define(rD, rIn)

//Plane E:
m4_define(zE, diffuserLength)
m4_define(rE, rOut)
m4_define(rHalfDumpE, calc(2*rOut))
m4_define(rDumpE, rDump)

//Plane F:
m4_define(zF, calc(diffuserLength+dumpLength/2))
m4_define(rF, rOut)
m4_define(rHalfDumpF, calc(2*rOut))
m4_define(rDumpF, rDump)

//Plane G:
m4_define(zG, calc(diffuserLength+dumpLength))
m4_define(rG, rOut)
m4_define(rHalfDumpG, calc(2*rOut))
m4_define(rDumpG, rDump)

//Plane H:
m4_define(zH, calc(diffuserLength+dumpLength+outletPipeLength))
m4_define(rH, rOut)

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
(0 0 zA) vlabel(A0)
(rA 0 zA) vlabel(A1)
(0 rA zA) vlabel(A2)
(-rA 0 zA) vlabel(A3)
(0 -rA zA) vlabel(A4)

//Plane B:
(0 0 zB) vlabel(B0)
(rB 0 zB) vlabel(B1)
(0 rB zB) vlabel(B2)
(-rB 0 zB) vlabel(B3)
(0 -rB zB) vlabel(B4)

//Plane C:
(0 0 zC) vlabel(C0)
(rC 0 zC) vlabel(C1)
(0 rC zC) vlabel(C2)
(-rC 0 zC) vlabel(C3)
(0 -rC zC) vlabel(C4)

//Plane D:
(0 0 zD) vlabel(D0)
(rD 0 zD) vlabel(D1)
(0 rD zD) vlabel(D2)
(-rD 0 zD) vlabel(D3)
(0 -rD zD) vlabel(D4)

//Plane E:
(0 0 zE) vlabel(E0)
(rE 0 zE) vlabel(E1)
(0 rE zE) vlabel(E2)
(-rE 0 zE) vlabel(E3)
(0 -rE zE) vlabel(E4)
(rHalfDumpE 0 zE) vlabel(E5)
(0 rHalfDumpE zE) vlabel(E6)
(-rHalfDumpE 0 zE) vlabel(E7)
(0 -rHalfDumpE zE) vlabel(E8)
(rDumpE 0 zE) vlabel(E9)
(0 rDumpE zE) vlabel(E10)
(-rDumpE 0 zE) vlabel(E11)
(0 -rDumpE zE) vlabel(E12)

//Plane F:
(0 0 zF) vlabel(F0)
(rF 0 zF) vlabel(F1)
(0 rF zF) vlabel(F2)
(-rF 0 zF) vlabel(F3)
(0 -rF zF) vlabel(F4)
(rHalfDumpF 0 zF) vlabel(F5)
(0 rHalfDumpF zF) vlabel(F6)
(-rHalfDumpF 0 zF) vlabel(F7)
(0 -rHalfDumpF zF) vlabel(F8)
(rDumpF 0 zF) vlabel(F9)
(0 rDumpF zF) vlabel(F10)
(-rDumpF 0 zF) vlabel(F11)
(0 -rDumpF zF) vlabel(F12)

//Plane G:
(0 0 zG) vlabel(G0)
(rG 0 zG) vlabel(G1)
(0 rG zG) vlabel(G2)
(-rG 0 zG) vlabel(G3)
(0 -rG zG) vlabel(G4)
(rHalfDumpG 0 zG) vlabel(G5)
(0 rHalfDumpG zG) vlabel(G6)
(-rHalfDumpG 0 zG) vlabel(G7)
(0 -rHalfDumpG zG) vlabel(G8)
(rDumpG 0 zG) vlabel(G9)
(0 rDumpG zG) vlabel(G10)
(-rDumpG 0 zG) vlabel(G11)
(0 -rDumpG zG) vlabel(G12)

//Plane H:
(0 0 zH) vlabel(H0)
(rH 0 zH) vlabel(H1)
(0 rH zH) vlabel(H2)
(-rH 0 zH) vlabel(H3)
(0 -rH zH) vlabel(H4)
);

// Defining blocks:
blocks
(
    //Blocks between plane A and plane B:
    // block0 - positive x and y 
    hex (A0 A1 A2 A0 B0 B1 B2 B0) AB
    (rNumberOfCells1st tNumberOfCells zABnumberOfCells)
    simpleGrading (rGrading2 1 1)
    // block1 - negative x positive y 
    hex (A0 A2 A3 A0 B0 B2 B3 B0) AB
    (rNumberOfCells1st tNumberOfCells zABnumberOfCells)
    simpleGrading (rGrading2 1 1)
    // block2 - negative x and  y 
    hex (A0 A3 A4 A0 B0 B3 B4 B0) AB
    (rNumberOfCells1st tNumberOfCells zABnumberOfCells)
    simpleGrading (rGrading2 1 1)
    // block3 - positive x negative y 
    hex (A0 A4 A1 A0 B0 B4 B1 B0) AB
    (rNumberOfCells1st tNumberOfCells zABnumberOfCells)
    simpleGrading (rGrading2 1 1)

    //Blocks between plane B and plane C:
    // block0 - positive x and y 
    hex (B0 B1 B2 B0 C0 C1 C2 C0) BC
    (rNumberOfCells1st tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading2 1 1)
    // block1 - negative x positive y 
    hex (B0 B2 B3 B0 C0 C2 C3 C0) BC
    (rNumberOfCells1st tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading2 1 1)
    // block2 - negative x and  y 
    hex (B0 B3 B4 B0 C0 C3 C4 C0) BC
    (rNumberOfCells1st tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading2 1 1)
    // block3 - positive x negative y 
    hex (B0 B4 B1 B0 C0 C4 C1 C0) BC
    (rNumberOfCells1st tNumberOfCells zBCnumberOfCells)
    simpleGrading (rGrading2 1 1)

    //Blocks between plane C and plane D:
    // block0 - positive x and y 
    hex (C0 C1 C2 C0 D0 D1 D2 D0) CD
    (rNumberOfCells1st tNumberOfCells zCDnumberOfCells)
    simpleGrading (rGrading2 1 1)
    // block1 - negative x positive y 
    hex (C0 C2 C3 C0 D0 D2 D3 D0) CD
    (rNumberOfCells1st tNumberOfCells zCDnumberOfCells)
    simpleGrading (rGrading2 1 1)
    // block2 - negative x and  y 
    hex (C0 C3 C4 C0 D0 D3 D4 D0) CD
    (rNumberOfCells1st tNumberOfCells zCDnumberOfCells)
    simpleGrading (rGrading2 1 1)
    // block3 - positive x negative y 
    hex (C0 C4 C1 C0 D0 D4 D1 D0) CD
    (rNumberOfCells1st tNumberOfCells zCDnumberOfCells)
    simpleGrading (rGrading2 1 1)

    //Blocks between plane D and plane E:
    // block0 - positive x and y 
    hex (D0 D1 D2 D0 E0 E1 E2 E0) DE
    (rNumberOfCells1st tNumberOfCells zDEnumberOfCells)
    simpleGrading (rGrading2 1 1)
    // block1 - negative x positive y 
    hex (D0 D2 D3 D0 E0 E2 E3 E0) DE
    (rNumberOfCells1st tNumberOfCells zDEnumberOfCells)
    simpleGrading (rGrading2 1 1)
    // block2 - negative x and  y 
    hex (D0 D3 D4 D0 E0 E3 E4 E0) DE
    (rNumberOfCells1st tNumberOfCells zDEnumberOfCells)
    simpleGrading (rGrading2 1 1)
    // block3 - positive x negative y 
    hex (D0 D4 D1 D0 E0 E4 E1 E0) DE
    (rNumberOfCells1st tNumberOfCells zDEnumberOfCells)
    simpleGrading (rGrading2 1 1)

    //Blocks between plane E and plane F:
    // block0 - positive x and y 
    hex (E0 E1 E2 E0 F0 F1 F2 F0) EF
    (rNumberOfCells1st tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading2 1 zGrading1)
    // block1 - negative x positive y 
    hex (E0 E2 E3 E0 F0 F2 F3 F0) EF
    (rNumberOfCells1st tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading2 1 zGrading1)
    // block2 - negative x and  y 
    hex (E0 E3 E4 E0 F0 F3 F4 F0) EF
    (rNumberOfCells1st tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading2 1 zGrading1)
    // block3 - positive x negative y 
    hex (E0 E4 E1 E0 F0 F4 F1 F0) EF
    (rNumberOfCells1st tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading2 1 zGrading1)
    // block4 - positive x and y 
    hex (E1 E5 E6 E2 F1 F5 F6 F2) EF
    (rNumberOfCells2nd tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading1 1 zGrading1)
    // block5 - negative x positive y 
    hex (E2 E6 E7 E3 F2 F6 F7 F3) EF
    (rNumberOfCells2nd tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading1 1 zGrading1)
    // block6 - negative x and  y 
    hex (E3 E7 E8 E4 F3 F7 F8 F4) EF
    (rNumberOfCells2nd tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading1 1 zGrading1)
    // block7 - positive x negative y 
    hex (E4 E8 E5 E1 F4 F8 F5 F1) EF
    (rNumberOfCells2nd tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading1 1 zGrading1)
    // block8 - positive x and y 
    hex (E5 E9 E10 E6 F5 F9 F10 F6) EF
    (rNumberOfCells3rd tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading2 1 zGrading1)
    // block9 - negative x positive y 
    hex (E6 E10 E11 E7 F6 F10 F11 F7) EF
    (rNumberOfCells3rd tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading2 1 zGrading1)
    // block10 - negative x and  y 
    hex (E7 E11 E12 E8 F7 F11 F12 F8) EF
    (rNumberOfCells3rd tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading2 1 zGrading1)
    // block11 - positive x negative y 
    hex (E8 E12 E9 E5 F8 F12 F9 F5) EF
    (rNumberOfCells3rd tNumberOfCells zEFnumberOfCells)
    simpleGrading (rGrading2 1 zGrading1)

//Blocks between plane F and plane G:
    // block0 - positive x and y 
    hex (F0 F1 F2 F0 G0 G1 G2 G0) FG
    (rNumberOfCells1st tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading2 1 zGrading2)
    // block1 - negative x positive y 
    hex (F0 F2 F3 F0 G0 G2 G3 G0) FG
    (rNumberOfCells1st tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading2 1 zGrading2)
    // block2 - negative x and  y 
    hex (F0 F3 F4 F0 G0 G3 G4 G0) FG
    (rNumberOfCells1st tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading2 1 zGrading2)
    // block3 - positive x negative y 
    hex (F0 F4 F1 F0 G0 G4 G1 G0) FG
    (rNumberOfCells1st tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading2 1 zGrading2)
    // block4 - positive x and y 
    hex (F1 F5 F6 F2 G1 G5 G6 G2) FG
    (rNumberOfCells2nd tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading1 1 zGrading2)
    // block5 - negative x positive y 
    hex (F2 F6 F7 F3 G2 G6 G7 G3) FG
    (rNumberOfCells2nd tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading1 1 zGrading2)
   // block6 - negative x and  y 
    hex (F3 F7 F8 F4 G3 G7 G8 G4) FG
    (rNumberOfCells2nd tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading1 1 zGrading2)
    // block7 - positive x negative y 
    hex (F4 F8 F5 F1 G4 G8 G5 G1) FG
    (rNumberOfCells2nd tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading1 1 zGrading2)
    // block8 - positive x and y 
    hex (F5 F9 F10 F6 G5 G9 G10 G6) FG
    (rNumberOfCells3rd tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading2 1 zGrading2)
    // block9 - negative x positive y 
    hex (F6 F10 F11 F7 G6 G10 G11 G7) FG
    (rNumberOfCells3rd tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading2 1 zGrading2)
    // block10 - negative x and  y 
    hex (F7 F11 F12 F8 G7 G11 G12 G8) FG
    (rNumberOfCells3rd tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading2 1 zGrading2)
    // block11 - positive x negative y 
    hex (F8 F12 F9 F5 G8 G12 G9 G5) EF
    (rNumberOfCells3rd tNumberOfCells zFGnumberOfCells)
    simpleGrading (rGrading2 1 zGrading2)

//Blocks between plane G and plane H:
    // block0 - positive x and y 
    hex (G0 G1 G2 G0 H0 H1 H2 H0) GH
    (rNumberOfCells1st tNumberOfCells zDEnumberOfCells)
    simpleGrading (rGrading2 1 zGrading1)
    // block1 - negative x positive y 
    hex (G0 G2 G3 G0 H0 H2 H3 H0) GH
    (rNumberOfCells1st tNumberOfCells zDEnumberOfCells)
    simpleGrading (rGrading2 1 zGrading1)
    // block2 - negative x and  y 
    hex (G0 G3 G4 G0 H0 H3 H4 H0) GH
    (rNumberOfCells1st tNumberOfCells zDEnumberOfCells)
    simpleGrading (rGrading2 1 zGrading1)
    // block3 - positive x negative y 
    hex (G0 G4 G1 G0 H0 H4 H1 H0) GH
    (rNumberOfCells1st tNumberOfCells zDEnumberOfCells)
    simpleGrading (rGrading2 1 zGrading1)
);

edges
(
    //Plane A:
    arc A1 A2 (calc(rA*sqr) calc(rA*sqr) zA)
    arc A2 A3 (-calc(rA*sqr) calc(rA*sqr) zA)
    arc A3 A4 (-calc(rA*sqr) -calc(rA*sqr) zA)
    arc A4 A1 (calc(rA*sqr) -calc(rA*sqr) zA)

    //Plane B:
    arc B1 B2 (calc(rB*sqr) calc(rB*sqr) zB)
    arc B2 B3 (-calc(rB*sqr) calc(rB*sqr) zB)
    arc B3 B4 (-calc(rB*sqr) -calc(rB*sqr) zB)
    arc B4 B1 (calc(rB*sqr) -calc(rB*sqr) zB)

    //Plane C:
    arc C1 C2 (calc(rC*sqr) calc(rC*sqr) zC)
    arc C2 C3 (-calc(rC*sqr) calc(rC*sqr) zC)
    arc C3 C4 (-calc(rC*sqr) -calc(rC*sqr) zC)
    arc C4 C1 (calc(rC*sqr) -calc(rC*sqr) zC)

    //Plane D:
    arc D1 D2 (calc(rD*sqr) calc(rD*sqr) zD)
    arc D2 D3 (-calc(rD*sqr) calc(rD*sqr) zD)
    arc D3 D4 (-calc(rD*sqr) -calc(rD*sqr) zD)
    arc D4 D1 (calc(rD*sqr) -calc(rD*sqr) zD)

    //Plane E:
    arc E1 E2 (calc(rE*sqr) calc(rE*sqr) zE)
    arc E2 E3 (-calc(rE*sqr) calc(rE*sqr) zE)
    arc E3 E4 (-calc(rE*sqr) -calc(rE*sqr) zE)
    arc E4 E1 (calc(rE*sqr) -calc(rE*sqr) zE)
    arc E5 E6 (calc(rHalfDumpE*sqr) calc(rHalfDumpE*sqr) zE)
    arc E6 E7 (-calc(rHalfDumpE*sqr) calc(rHalfDumpE*sqr) zE)
    arc E7 E8 (-calc(rHalfDumpE*sqr) -calc(rHalfDumpE*sqr) zE)
    arc E8 E5 (calc(rHalfDumpE*sqr) -calc(rHalfDumpE*sqr) zE)
    arc E9 E10 (calc(rDumpE*sqr) calc(rDumpE*sqr) zE)
    arc E10 E11 (-calc(rDumpE*sqr) calc(rDumpE*sqr) zE)
    arc E11 E12 (-calc(rDumpE*sqr) -calc(rDumpE*sqr) zE)
    arc E12 E9 (calc(rDumpE*sqr) -calc(rDumpE*sqr) zE)

    //Plane F:
    arc F1 F2 (calc(rF*sqr) calc(rF*sqr) zF)
    arc F2 F3 (-calc(rF*sqr) calc(rF*sqr) zF)
    arc F3 F4 (-calc(rF*sqr) -calc(rF*sqr) zF)
    arc F4 F1 (calc(rF*sqr) -calc(rF*sqr) zF)
    arc F5 F6 (calc(rHalfDumpF*sqr) calc(rHalfDumpF*sqr) zF)
    arc F6 F7 (-calc(rHalfDumpF*sqr) calc(rHalfDumpF*sqr) zF)
    arc F7 F8 (-calc(rHalfDumpF*sqr) -calc(rHalfDumpF*sqr) zF)
    arc F8 F5 (calc(rHalfDumpF*sqr) -calc(rHalfDumpF*sqr) zF)
    arc F9 F10 (calc(rDumpF*sqr) calc(rDumpF*sqr) zF)
    arc F10 F11 (-calc(rDumpF*sqr) calc(rDumpF*sqr) zF)
    arc F11 F12 (-calc(rDumpF*sqr) -calc(rDumpF*sqr) zF)
    arc F12 F9 (calc(rDumpF*sqr) -calc(rDumpF*sqr) zF)

    //Plane G:
    arc G1 G2 (calc(rG*sqr) calc(rG*sqr) zG)
    arc G2 G3 (-calc(rG*sqr) calc(rG*sqr) zG)
    arc G3 G4 (-calc(rG*sqr) -calc(rG*sqr) zG)
    arc G4 G1 (calc(rG*sqr) -calc(rG*sqr) zG)
    arc G5 G6 (calc(rHalfDumpG*sqr) calc(rHalfDumpG*sqr) zG)
    arc G6 G7 (-calc(rHalfDumpG*sqr) calc(rHalfDumpG*sqr) zG)
    arc G7 G8 (-calc(rHalfDumpG*sqr) -calc(rHalfDumpG*sqr) zG)
    arc G8 G5 (calc(rHalfDumpG*sqr) -calc(rHalfDumpG*sqr) zG)
    arc G9 G10 (calc(rDumpG*sqr) calc(rDumpG*sqr) zG)
    arc G10 G11 (-calc(rDumpG*sqr) calc(rDumpG*sqr) zG)
    arc G11 G12 (-calc(rDumpG*sqr) -calc(rDumpG*sqr) zG)
    arc G12 G9 (calc(rDumpG*sqr) -calc(rDumpG*sqr) zG)

    //Plane H:
    arc H1 H2 (calc(rH*sqr) calc(rH*sqr) zH)
    arc H2 H3 (-calc(rH*sqr) calc(rH*sqr) zH)
    arc H3 H4 (-calc(rH*sqr) -calc(rH*sqr) zH)
    arc H4 H1 (calc(rH*sqr) -calc(rH*sqr) zH)
);

// Defining patches:
patches
(
    patch inlet
    (
       (A0 A2 A1 A0)
       (A0 A1 A4 A0)
       (A0 A4 A3 A0)
       (A0 A3 A2 A0)
    )
    patch outlet
    (
       (H0 H1 H2 H0)
       (H0 H2 H3 H0)
       (H0 H3 H4 H0)
       (H0 H4 H1 H0)
    )
    wall wallDump
    (	
       (E1 E2 E6 E5)
       (E2 E3 E7 E6)
       (E3 E4 E8 E7)
       (E4 E1 E5 E8)
       (E5 E6 E10 E9)
       (E6 E7 E11 E10)
       (E7 E8 E12 E11)
       (E8 E5 E9 E12) 
       (E9 E10 F10 F9)
       (E10 E11 F11 F10)
       (E11 E12 F12 F11)
       (E12 E9 F9 F12)
       (F9 F10 G10 G9)
       (F10 F11 G11 G10)
       (F11 F12 G12 G11)
       (F12 F9 G9 G12)
       (G1 G5 G6 G2)
       (G2 G6 G7 G3)
       (G3 G7 G8 G4)
       (G4 G8 G5 G1)
       (G5 G9 G10 G6)
       (G6 G10 G11 G7)
       (G7 G11 G12 G8)
       (G8 G12 G9 G5)

    ) 
wall wallProlongation
    (
       (G1 G2 H2 H1)
       (G2 G3 H3 H2)
       (G3 G4 H4 H3)
       (G4 G1 H1 H4)
    )
    wall wallDiffuser
    (
       (D1 D2 E2 E1)
       (D2 D3 E3 E2)
       (D3 D4 E4 E3)
       (D4 D1 E1 E4)
    )
    wall statSwirlWall
    (
       (B1 B2 C2 C1)
       (B2 B3 C3 C2)
       (B3 B4 C4 C3)
       (B4 B1 C1 C4)
       (C1 C2 D2 D1)
       (C2 C3 D3 D2)
       (C3 C4 D4 D3)
       (C4 C1 D1 D4)
    )
    wall rotSwirlWall
    (
       (A1 A2 B2 B1)
       (A2 A3 B3 B2)
       (A3 A4 B4 B3)
       (A4 A1 B1 B4)
    )
);

mergePatchPairs 
(
);

// ************************************************************************* //
