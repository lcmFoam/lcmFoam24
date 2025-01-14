/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "codedMixedFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = cc40f9179baf5094ae2b9c5b4870ff4fa2b49436
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void myBC7_cc40f9179baf5094ae2b9c5b4870ff4fa2b49436(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchScalarField,
    myBC7MixedValueFvPatchScalarField
);


const char* const myBC7MixedValueFvPatchScalarField::SHA1sum =
    "cc40f9179baf5094ae2b9c5b4870ff4fa2b49436";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

myBC7MixedValueFvPatchScalarField::
myBC7MixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF)
{
    if (false)
    {
        Info<<"construct myBC7 sha1: cc40f9179baf5094ae2b9c5b4870ff4fa2b49436"
            " from patch/DimensionedField\n";
    }
}


myBC7MixedValueFvPatchScalarField::
myBC7MixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct myBC7 sha1: cc40f9179baf5094ae2b9c5b4870ff4fa2b49436"
            " from patch/dictionary\n";
    }
}


myBC7MixedValueFvPatchScalarField::
myBC7MixedValueFvPatchScalarField
(
    const myBC7MixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct myBC7 sha1: cc40f9179baf5094ae2b9c5b4870ff4fa2b49436"
            " from patch/DimensionedField/mapper\n";
    }
}


myBC7MixedValueFvPatchScalarField::
myBC7MixedValueFvPatchScalarField
(
    const myBC7MixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct myBC7 sha1: cc40f9179baf5094ae2b9c5b4870ff4fa2b49436 "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

myBC7MixedValueFvPatchScalarField::
~myBC7MixedValueFvPatchScalarField()
{
    if (false)
    {
        Info<<"destroy myBC7 sha1: cc40f9179baf5094ae2b9c5b4870ff4fa2b49436\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void myBC7MixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs myBC7 sha1: cc40f9179baf5094ae2b9c5b4870ff4fa2b49436\n";
    }

//{{{ begin code
    #line 53 "/home/obertscheider/myOpenFOAM/cases_openfoam2312/FPCM16__4plus8cells_variablethickness/0/p/boundaryField/dynamic"
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

    this->mixedFvPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

