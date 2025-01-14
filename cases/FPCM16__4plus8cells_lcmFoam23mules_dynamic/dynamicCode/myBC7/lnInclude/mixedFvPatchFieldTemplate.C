/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
    Copyright (C) YEAR AUTHOR, AFFILIATION
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "mixedFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "PatchFunction1.H"

//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = 35c87f70b08012e8d8f55a7284aa2775df550c43
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void myBC7_35c87f70b08012e8d8f55a7284aa2775df550c43(bool load)
{
    if (load)
    {
        // Code that can be explicitly executed after loading
    }
    else
    {
        // Code that can be explicitly executed before unloading
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchScalarField,
    myBC7MixedValueFvPatchScalarField
);

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
myBC7MixedValueFvPatchScalarField::
myBC7MixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(p, iF)
{
    if (false)
    {
        printMessage("Construct myBC7 : patch/DimensionedField");
    }
}


Foam::
myBC7MixedValueFvPatchScalarField::
myBC7MixedValueFvPatchScalarField
(
    const myBC7MixedValueFvPatchScalarField& rhs,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(rhs, p, iF, mapper)
{
    if (false)
    {
        printMessage("Construct myBC7 : patch/DimensionedField/mapper");
    }
}


Foam::
myBC7MixedValueFvPatchScalarField::
myBC7MixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF, dict)
{
    if (false)
    {
        printMessage("Construct myBC7 : patch/dictionary");
    }
}


Foam::
myBC7MixedValueFvPatchScalarField::
myBC7MixedValueFvPatchScalarField
(
    const myBC7MixedValueFvPatchScalarField& rhs
)
:
    parent_bctype(rhs),
    dictionaryContent(rhs)
{
    if (false)
    {
        printMessage("Copy construct myBC7");
    }
}


Foam::
myBC7MixedValueFvPatchScalarField::
myBC7MixedValueFvPatchScalarField
(
    const myBC7MixedValueFvPatchScalarField& rhs,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(rhs, iF)
{
    if (false)
    {
        printMessage("Construct myBC7 : copy/DimensionedField");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::
myBC7MixedValueFvPatchScalarField::
~myBC7MixedValueFvPatchScalarField()
{
    if (false)
    {
        printMessage("Destroy myBC7");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::
myBC7MixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        printMessage("updateCoeffs myBC7");
    }

//{{{ begin code
    #line 52 "/mnt/d/work/LCM_OpenFOAM/OpenFOAM_in_WSL/FPCM16_cases/FPCM16__4plus8cells_lcmFoam23mules_dynamic/0/p_rgh/boundaryField/dynamic"
const tmp<scalarField>& p=patch().lookupPatchField<volScalarField, scalar>("p").patchInternalField(); 
            const tmp<scalarField>& alpha1=patch().lookupPatchField<volScalarField, scalar>("alpha.resin").patchInternalField();
            const fvPatch& boundaryPatch = patch(); 
            const vectorField& Cf = boundaryPatch.Cf(); 
            scalarField& val=this->refValue();
            scalarField& grad=this->refGrad();
            scalarField& frac=this->valueFraction();
            scalar fracVal=1.0;
            forAll(Cf, faceI)
            {
                frac[faceI]=1.0-min(max(alpha1.ref()[faceI],0.0),1.0);  //linear transition
            }
//}}} end code

    this->parent_bctype::updateCoeffs();
}


// ************************************************************************* //

