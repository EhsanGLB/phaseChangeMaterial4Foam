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

#include "timeVaryingMixedRobin4PCMFoam.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeVaryingMixedRobin4PCMFoam::timeVaryingMixedRobin4PCMFoam
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    TName_("T"),
    kappaName_("undefined-Kappa"),
    absorptivity_(0.0),
    emissivity_(0.0),
    timeSeriesqs_(),
    timeSeriesho_(),
    timeSeriesTo_(),
    timeSeriesTsky_()
{}


Foam::timeVaryingMixedRobin4PCMFoam::timeVaryingMixedRobin4PCMFoam
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    kappaName_(dict.lookup("Kappa")),
    absorptivity_(readScalar(dict.lookup("absorptivity"))),
    emissivity_(readScalar(dict.lookup("emissivity"))),
    timeSeriesqs_(dict.subDict("qs")),
    timeSeriesho_(dict.subDict("ho")),
    timeSeriesTo_(dict.subDict("To")),
    timeSeriesTsky_(dict.subDict("Tsky"))
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


Foam::timeVaryingMixedRobin4PCMFoam::timeVaryingMixedRobin4PCMFoam
(
    const timeVaryingMixedRobin4PCMFoam& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    kappaName_(ptf.kappaName_),
    absorptivity_(ptf.absorptivity_),
    emissivity_(ptf.emissivity_),
    timeSeriesqs_(ptf.timeSeriesqs_),
    timeSeriesho_(ptf.timeSeriesho_),
    timeSeriesTo_(ptf.timeSeriesTo_),
    timeSeriesTsky_(ptf.timeSeriesTsky_)
{}


Foam::timeVaryingMixedRobin4PCMFoam::timeVaryingMixedRobin4PCMFoam
(
    const timeVaryingMixedRobin4PCMFoam& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    kappaName_(ptf.kappaName_),
    absorptivity_(ptf.absorptivity_),
    emissivity_(ptf.emissivity_),
    timeSeriesqs_(ptf.timeSeriesqs_),
    timeSeriesho_(ptf.timeSeriesho_),
    timeSeriesTo_(ptf.timeSeriesTo_),
    timeSeriesTsky_(ptf.timeSeriesTsky_)
{}


Foam::timeVaryingMixedRobin4PCMFoam::timeVaryingMixedRobin4PCMFoam
(
    const timeVaryingMixedRobin4PCMFoam& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    kappaName_(ptf.kappaName_),
    absorptivity_(ptf.absorptivity_),
    emissivity_(ptf.emissivity_),
    timeSeriesqs_(ptf.timeSeriesqs_),
    timeSeriesho_(ptf.timeSeriesho_),
    timeSeriesTo_(ptf.timeSeriesTo_),
    timeSeriesTsky_(ptf.timeSeriesTsky_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeVaryingMixedRobin4PCMFoam::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //constant parameters
    const scalar sigmaSB_(5.67e-8);


    //- get the currernt time and time step
    scalar t_ = this->db().time().timeOutputValue();


    const fvPatchScalarField& Tp = lookupPatchField<volScalarField, scalar>(TName_);
    const fvPatchScalarField& kappap = lookupPatchField<volScalarField, scalar>(kappaName_);

    gradient() = ( absorptivity_ * timeSeriesqs_(t_) - timeSeriesho_(t_) * (Tp - timeSeriesTo_(t_)) - emissivity_ * sigmaSB_ * (pow(Tp, 4) - pow(timeSeriesTsky_(t_), 4)) )/kappap;
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::timeVaryingMixedRobin4PCMFoam::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("Kappa") << kappaName_ << token::END_STATEMENT << nl;
    os.writeKeyword("absorptivity") << absorptivity_ << token::END_STATEMENT << nl;
    os.writeKeyword("emissivity") << emissivity_ << token::END_STATEMENT << nl;
    /*timeSeriesqs_.write(os);
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
        timeVaryingMixedRobin4PCMFoam
    );
}

// ************************************************************************* //
