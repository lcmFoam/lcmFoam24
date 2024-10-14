    IOporosityModelList pZones(mesh);
    Switch pressureImplicitPorosity(false);



volScalarField porosity
(
    IOobject
    (
        "porosity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("NULL", dimless, 1.0),   
    "zeroGradient"
);
porosity.write();
Info << "CO: Reading porosity values and assigning porosity field for LCM filling simulation using Darcy-Forchheimer model" << endl;
IOdictionary porousZones
(
    IOobject
    (
        "porosityProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);
label count = 0;
forAllConstIter(dictionary,porousZones,iter)
{
    if (iter().isDict())
    {
        count++;
    }
}


//Info<< "count: " << count << endl;
label i = 0;
forAllConstIter(dictionary,porousZones,iter)
{
    if (iter().isDict())
    {
        i++;
        const dictionary& modelDict = iter().dict();
        scalar porosity_=readScalar(modelDict.lookup("porosity"));
        word cellZone_ = modelDict.lookup("cellZone");
        Info << "    cellZone(" << i << ")=" << cellZone_ << endl;
        Info << "    porosity(" << i << ")=" << porosity_ << endl;

        const label cellZoneID = mesh.cellZones().findZoneID(cellZone_);    
        const labelList& cells = mesh.cellZones()[cellZoneID];

        Info << "Cells in cellzone " << i << ":" << endl;
        forAll(cells, i)
        {
            const label cell = cells[i];
            Info << cell << endl;
            porosity[cell]=porosity_;
        }
        
        
        //const label cellZoneID = mesh.cellZones().findZoneID(cellZone_);
        //const cellZone& zone = mesh.cellZones()[cellZoneID];
        //const cellZoneMesh& zoneMesh = zone.zoneMesh();
        //const labelList& cellZone = zoneMesh[cellZoneID];  // list of all cell ID's
        //forAll(cellZone,cellI)
        //{
        //    porosity[cellZone[cellI]]=porosity_;
        //}
    }
}


