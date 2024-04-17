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

#include "timeVaryingRobin4PCMFoam.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeVaryingRobin4PCMFoam::timeVaryingRobin4PCMFoam
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    TName_("T"),
    kappaName_("undefined-Kappa"),
    timeSeriesqo_(),
    timeSeriesho_(),
    timeSeriesTo_()
{}


Foam::timeVaryingRobin4PCMFoam::timeVaryingRobin4PCMFoam
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    kappaName_(dict.lookup("Kappa")),
    timeSeriesqo_(dict.subDict("qo")),
    timeSeriesho_(dict.subDict("ho")),
    timeSeriesTo_(dict.subDict("To"))
{
    if (dict.found("gradient"))
    {
        gradient() = scalarField("gradient", dict, p.size());
        fixedGradientFvPatchScalarField::updateCoeffs();
        fixedGradientFvPatchScalarField::evaluate();
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


Foam::timeVaryingRobin4PCMFoam::timeVaryingRobin4PCMFoam
(
    const timeVaryingRobin4PCMFoam& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    kappaName_(ptf.kappaName_),
    timeSeriesqo_(ptf.timeSeriesqo_),
    timeSeriesho_(ptf.timeSeriesho_),
    timeSeriesTo_(ptf.timeSeriesTo_)
{}


Foam::timeVaryingRobin4PCMFoam::timeVaryingRobin4PCMFoam
(
    const timeVaryingRobin4PCMFoam& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    kappaName_(ptf.kappaName_),
    timeSeriesqo_(ptf.timeSeriesqo_),
    timeSeriesho_(ptf.timeSeriesho_),
    timeSeriesTo_(ptf.timeSeriesTo_)
{}


Foam::timeVaryingRobin4PCMFoam::timeVaryingRobin4PCMFoam
(
    const timeVaryingRobin4PCMFoam& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    kappaName_(ptf.kappaName_),
    timeSeriesqo_(ptf.timeSeriesqo_),
    timeSeriesho_(ptf.timeSeriesho_),
    timeSeriesTo_(ptf.timeSeriesTo_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeVaryingRobin4PCMFoam::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //- get the currernt time and time step
    scalar t_ = this->db().time().timeOutputValue();


    const fvPatchScalarField& Tp = lookupPatchField<volScalarField, scalar>(TName_);
    const fvPatchScalarField& kappap = lookupPatchField<volScalarField, scalar>(kappaName_);

    gradient() = ( timeSeriesqo_(t_) - timeSeriesho_(t_) * (Tp - timeSeriesTo_(t_)) )/kappap;
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::timeVaryingRobin4PCMFoam::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("Kappa") << kappaName_ << token::END_STATEMENT << nl;
    /*timeSeriesqo_.write(os);
    timeSeriesho_.write(os);
    timeSeriesTo_.write(os);*/

    writeEntry("value", os);
    gradient().writeEntry("gradient", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        timeVaryingRobin4PCMFoam
    );
}

// ************************************************************************* //
