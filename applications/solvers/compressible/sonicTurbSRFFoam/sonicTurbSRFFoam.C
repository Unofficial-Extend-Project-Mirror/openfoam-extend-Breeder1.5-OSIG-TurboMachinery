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

Application
    sonicTurbSRFFoam

Description
    Transient solver for the RELATIVE, trans-sonic/supersonic, turbulent flow of a
    compressible gas within multiple rotating frames.

Author
    1991-2008 OpenCFD Ltd.
    2009 Oliver Borm <oli.borm@web.de>

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "basicThermo.H"
#include "compressible/RASModel/RASModel.H"
#include "SRFZones.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"

#   include "readTimeControls.H"
#   include "compressibleCourantNo.H"
#   include "setInitialDeltaT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readTimeControls.H"
#       include "readPISOControls.H"
#       include "compressibleCourantNo.H"
#       include "setDeltaT.H"

#       include "rhoEqn.H"

        fvVectorMatrix WxyzEqn
        (
            fvm::ddt(rho, Wxyz)
          + fvm::div(phi, Wxyz)
          + turbulence->divDevRhoReff(Wxyz)
        );

        srfZones.addSu(WxyzEqn);

        solve(WxyzEqn == -fvc::grad(p));

        solve
        (
            fvm::ddt(rho, h)
          + fvm::div(phi, h)
          - fvm::laplacian(turbulence->alphaEff(), h)
         ==
            DpDt
        );

        thermo->correct();

        // --- PISO loop

        for (int corr=0; corr<nCorr; corr++)
        {
            volScalarField rWxyzA = 1.0/WxyzEqn.A();
            Wxyz = rWxyzA*WxyzEqn.H();

            surfaceScalarField phid
            (
                "phid",
                fvc::interpolate(thermo->psi())
               *(
                   (fvc::interpolate(Wxyz) & mesh.Sf())
                 + fvc::ddtPhiCorr(rWxyzA, rho, Wxyz, phi)
               )
            );

            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::ddt(psi, p)
                  + fvm::div(phid, p)
                  - fvm::laplacian(rho*rWxyzA, p)
                );

                pEqn.solve();

                if (nonOrth == nNonOrthCorr)
                {
                    phi = pEqn.flux();
                }
            }

#           include "compressibleContinuityErrs.H"

            Wxyz -= rWxyzA*fvc::grad(p);
            Wxyz.correctBoundaryConditions();
        }

        DpDt = 
            fvc::DDt(surfaceScalarField("phiWxyz", phi/fvc::interpolate(rho)), p);

        turbulence->correct();

        rho = psi*p;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
