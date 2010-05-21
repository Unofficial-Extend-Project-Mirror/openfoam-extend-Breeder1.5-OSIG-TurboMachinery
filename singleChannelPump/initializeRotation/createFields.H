Info<< "Reading Initialization Properties\n" << endl;

IOdictionary initProperties
(
    IOobject
    (
        "initializationProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedVector dimensionedOmega
( 
  initProperties.lookup("omega")
);

dimensionedVector dimensionedOrigin
( 
  initProperties.lookup("origin")
);

dimensionedScalar dimensionedRefRadius
( 
  initProperties.lookup("refRadius")
);



vector omega = dimensionedOmega.value();
vector origin = dimensionedOrigin.value();
scalar refRadius = dimensionedRefRadius.value();
vector rotationAxis = omega / mag(omega);

Info << " Omega = " << omega << endl;
Info << " reference radius = " << refRadius << endl;

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
