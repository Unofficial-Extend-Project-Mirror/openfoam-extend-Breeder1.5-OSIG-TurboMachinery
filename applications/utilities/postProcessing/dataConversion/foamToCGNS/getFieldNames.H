/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright Hydro-Quebec - IREQ, 2008
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
  Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description
  Conversion of CGNS files into Foam's mesh and fields

Authors 
  Robert Magnan,   Hydro-Quebec - IREQ, 2007

\*---------------------------------------------------------------------------*/


{ 
	HashSet<word> volScalarHash;
	HashSet<word> volVectorHash;
	HashSet<word> surfScalarHash;
	HashSet<word> surfVectorHash;
	HashSet<word> sprayScalarHash;
	HashSet<word> sprayVectorHash;

// Add all fields to hashtable
	IOobjectList objects(mesh, runTime.timeName());
	{
		wordList fieldNames(objects.names(volScalarField::typeName));
		forAll(fieldNames, fieldI)
		{
			volScalarHash.insert(fieldNames[fieldI]);
		}
	}
	{
		wordList fieldNames(objects.names(volVectorField::typeName));
		forAll(fieldNames, fieldI)
		{
			volVectorHash.insert(fieldNames[fieldI]);
		}
	}
	{
		wordList fieldNames(objects.names(surfaceScalarField::typeName));
		forAll(fieldNames, fieldI)
		{
			surfScalarHash.insert(fieldNames[fieldI]);
		}
	}
	{
		wordList fieldNames(objects.names(surfaceVectorField::typeName));
		forAll(fieldNames, fieldI)
		{
			surfVectorHash.insert(fieldNames[fieldI]);
		}
	}


// Same for spray
	IOobjectList sprayObjects(mesh, runTime.timeName(), "lagrangian");
	{
		wordList fieldNames(sprayObjects.names(scalarIOField::typeName));
		forAll(fieldNames, fieldI)
		{
			sprayScalarHash.insert(fieldNames[fieldI]);
		}
	}
	{
		wordList fieldNames(sprayObjects.names(vectorIOField::typeName));
		forAll(fieldNames, fieldI)
		{
			sprayVectorHash.insert(fieldNames[fieldI]);
		}
	}

        // Not great, but this will work with either the old or new
	// implementations of HashSet under OpenFOAM-1.5-dev.
	// MB 05/2010
	List<word> volScalarHashToc   = volScalarHash.toc();
	List<word> volVectorHashToc   = volVectorHash.toc();
	List<word> surfScalarHashToc  = surfScalarHash.toc();
	List<word> surfVectorHashToc  = surfVectorHash.toc();
	List<word> sprayScalarHashToc = sprayScalarHash.toc();
	List<word> sprayVectorHashToc = sprayVectorHash.toc();

	volScalarNames   = HashSet<word>(volScalarHashToc);
	volVectorNames   = HashSet<word>(volVectorHashToc);
	surfScalarNames  = HashSet<word>(surfScalarHashToc);
	surfVectorNames  = HashSet<word>(surfVectorHashToc);
	sprayScalarNames = HashSet<word>(sprayScalarHashToc);
	sprayVectorNames = HashSet<word>(sprayVectorHashToc);

#if 0
	// For OpenFOAM 1.6-dev, we should use this instead
	// MB 05/2010
	volScalarNames   = HashSet<word>(volScalarHash.toc());
	volVectorNames   = HashSet<word>(volVectorHash.toc());
	surfScalarNames  = HashSet<word>(surfScalarHash.toc());
	surfVectorNames  = HashSet<word>(surfVectorHash.toc());
	sprayScalarNames = HashSet<word>(sprayScalarHash.toc());
	sprayVectorNames = HashSet<word>(sprayVectorHash.toc());
#endif

}
//Info << "volScalarNames  : " << volScalarNames << endl;
//Info << "volVectorNames  : " << volVectorNames << endl;
//Info << "surfScalarNames : " << surfScalarNames << endl;
//Info << "surfVectorNames : " << surfVectorNames << endl;
//Info << "sprayScalarNames: " << sprayScalarNames << endl;
//Info << "sprayVectorNames: " << sprayVectorNames << endl;
