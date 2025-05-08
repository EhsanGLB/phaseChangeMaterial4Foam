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

#include "buongiorno4PCMFoam.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::buongiorno4PCMFoam::buongiorno4PCMFoam
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    TName_("T"),
    alphaName_("alpha"),
    jo_(0.0),
    ho_(0.0),
    alphao_(0.0),
    DB_(0.0),
    DT_(0.0)
{}


Foam::buongiorno4PCMFoam::buongiorno4PCMFoam
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    alphaName_(dict.lookupOrDefault<word>("alpha", "alpha")),
    jo_(),
    ho_(),
    alphao_(),
    DB_(readScalar(dict.lookup("DB"))),
    DT_(readScalar(dict.lookup("DT")))
{
    Istream& isjo_ = dict.lookup("jo");
    isjo_.format(IOstream::ASCII);
    isjo_ >> jo_;

    Istream& isho_ = dict.lookup("ho");
    isho_.format(IOstream::ASCII);
    isho_ >> ho_;

    Istream& isalphao_ = dict.lookup("alphao");
    isalphao_.format(IOstream::ASCII);
    isalphao_ >> alphao_;

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


Foam::buongiorno4PCMFoam::buongiorno4PCMFoam
(
    const buongiorno4PCMFoam& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    alphaName_(ptf.alphaName_),
    jo_(ptf.jo_),
    ho_(ptf.ho_),
    alphao_(ptf.alphao_),
    DB_(ptf.DB_),
    DT_(ptf.DT_)
{}


Foam::buongiorno4PCMFoam::buongiorno4PCMFoam
(
    const buongiorno4PCMFoam& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    alphaName_(ptf.alphaName_),
    jo_(ptf.jo_),
    ho_(ptf.ho_),
    alphao_(ptf.alphao_),
    DB_(ptf.DB_),
    DT_(ptf.DT_)
{}


Foam::buongiorno4PCMFoam::buongiorno4PCMFoam
(
    const buongiorno4PCMFoam& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    alphaName_(ptf.alphaName_),
    jo_(ptf.jo_),
    ho_(ptf.ho_),
    alphao_(ptf.alphao_),
    DB_(ptf.DB_),
    DT_(ptf.DT_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::buongiorno4PCMFoam::joFunction(scalar t)
{
    int n_ = 0;
    scalar f_ = 0.0;

    while( n_ < jo_.size() )
    {
        f_ += jo_[n_]*pow(t, n_);
        n_ += 1;
    }

    return f_;
}


Foam::scalar Foam::buongiorno4PCMFoam::hoFunction(scalar t)
{
    int n_ = 0;
    scalar f_ = 0.0;

    while( n_ < ho_.size() )
    {
        f_ += ho_[n_]*pow(t, n_);
        n_ += 1;
    }

    return f_;
}


Foam::scalar Foam::buongiorno4PCMFoam::alphaoFunction(scalar t)
{
    int n_ = 0;
    scalar f_ = 0.0;

    while( n_ < alphao_.size() )
    {
        f_ += alphao_[n_]*pow(t, n_);
        n_ += 1;
    }

    return f_;
}


void Foam::buongiorno4PCMFoam::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //- get the currernt time and time step
    scalar t_ = this->db().time().value();

    const fvPatchScalarField& Tp = lookupPatchField<volScalarField, scalar>(TName_);
    const fvPatchScalarField& alphap = lookupPatchField<volScalarField, scalar>(alphaName_);

    gradient() = joFunction(t_) -  hoFunction(t_) * (alphap - alphaoFunction(t_)) - ((DT_/Tp)/DB_)*Tp.snGrad();
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::buongiorno4PCMFoam::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("alpha") << alphaName_ << token::END_STATEMENT << nl;
    os.writeKeyword("jo") << jo_ << token::END_STATEMENT << nl;
    os.writeKeyword("ho") << ho_ << token::END_STATEMENT << nl;
    os.writeKeyword("alphao") << alphao_ << token::END_STATEMENT << nl;
    os.writeKeyword("DB") << DB_ << token::END_STATEMENT << nl;
    os.writeKeyword("DT") << DT_ << token::END_STATEMENT << nl;

    writeEntry("value", os);
    gradient().writeEntry("gradient", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        buongiorno4PCMFoam
    );
}

// ************************************************************************* //
