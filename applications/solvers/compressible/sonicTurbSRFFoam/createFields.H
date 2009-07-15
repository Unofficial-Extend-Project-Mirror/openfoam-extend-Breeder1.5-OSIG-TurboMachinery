    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<basicThermo> thermo
    (
        basicThermo::New(mesh)
    );

    volScalarField& p = thermo->p();
    volScalarField& h = thermo->h();
    const volScalarField& psi = thermo->psi();

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh
        ),
        thermo->rho()
    );

    Info<< "Reading field Wxyz\n" << endl;
    volVectorField Wxyz
    (
        IOobject
        (
            "Wxyz",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading/calculating face flux field phi\n" << endl;
    surfaceScalarField phi
    (
    IOobject
    (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(rho*Wxyz) & mesh.Sf()
    );

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::RASModel> turbulence
    (
        compressible::RASModel::New
        (
            rho,
            Wxyz,
            phi,
            thermo()
        )
    );

    Info<< "Creating field DpDt\n" << endl;
    volScalarField DpDt = 
        fvc::DDt(surfaceScalarField("phiWxyz", phi/fvc::interpolate(rho)), p);

    Info<< "Creating SRFZones model\n" << endl;
    SRFZones srfZones(mesh);
