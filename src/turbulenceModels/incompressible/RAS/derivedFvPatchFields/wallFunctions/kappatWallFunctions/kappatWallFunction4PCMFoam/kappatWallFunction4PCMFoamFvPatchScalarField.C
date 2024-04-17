/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "kappatWallFunction4PCMFoamFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void kappatWallFunction4PCMFoamFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    writeEntryIfDifferent<word>(os, "nut", "nut", nutName_);
    os.writeKeyword("Prt") << Prt_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kappatWallFunction4PCMFoamFvPatchScalarField::kappatWallFunction4PCMFoamFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    nutName_("nut"),
    Prt_(0.85)
{}


kappatWallFunction4PCMFoamFvPatchScalarField::kappatWallFunction4PCMFoamFvPatchScalarField
(
    const kappatWallFunction4PCMFoamFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    nutName_(ptf.nutName_),
    Prt_(ptf.Prt_)
{}


kappatWallFunction4PCMFoamFvPatchScalarField::kappatWallFunction4PCMFoamFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    nutName_(dict.lookupOrDefault<word>("nut", "nut")),
    Prt_(dict.lookupOrDefault<scalar>("Prt", 0.85))
{}


kappatWallFunction4PCMFoamFvPatchScalarField::kappatWallFunction4PCMFoamFvPatchScalarField
(
    const kappatWallFunction4PCMFoamFvPatchScalarField& awfpsf
)
:
    fixedValueFvPatchScalarField(awfpsf),
    nutName_(awfpsf.nutName_),
    Prt_(awfpsf.Prt_)
{}


kappatWallFunction4PCMFoamFvPatchScalarField::kappatWallFunction4PCMFoamFvPatchScalarField
(
    const kappatWallFunction4PCMFoamFvPatchScalarField& awfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(awfpsf, iF),
    nutName_(awfpsf.nutName_),
    Prt_(awfpsf.Prt_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void kappatWallFunction4PCMFoamFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField& nutw = lookupPatchField<volScalarField, scalar>(nutName_);
    const scalarField& rhow = lookupPatchField<volScalarField, scalar>("rho");
    const scalarField& Cpw = lookupPatchField<volScalarField, scalar>("Cp");

    operator==(rhow*Cpw*nutw/Prt_);
    fixedValueFvPatchScalarField::updateCoeffs();
}


void kappatWallFunction4PCMFoamFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, kappatWallFunction4PCMFoamFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
