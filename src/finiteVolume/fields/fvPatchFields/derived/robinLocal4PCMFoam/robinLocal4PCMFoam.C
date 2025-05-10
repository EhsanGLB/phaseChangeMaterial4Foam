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

#include "robinLocal4PCMFoam.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::robinLocal4PCMFoam::robinLocal4PCMFoam
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    TName_("T"),
    kappaName_("undefined-Kappa"),
    qo_(0.0),
    ho_(0.0),
    To_(0.0),
    dir_(0, 0, 0),
    p1_(0.0),
    p2_(0.0)
{}


Foam::robinLocal4PCMFoam::robinLocal4PCMFoam
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    kappaName_(dict.lookup("Kappa")),
    qo_(),
    ho_(),
    To_(),
    dir_(dict.lookup("dir")),
    p1_(readScalar(dict.lookup("p1"))),
    p2_(readScalar(dict.lookup("p2")))
{
    Istream& isqo_ = dict.lookup("qo");
    isqo_.format(IOstream::ASCII);
    isqo_ >> qo_;

    Istream& isho_ = dict.lookup("ho");
    isho_.format(IOstream::ASCII);
    isho_ >> ho_;

    Istream& isTo_ = dict.lookup("To");
    isTo_.format(IOstream::ASCII);
    isTo_ >> To_;

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


Foam::robinLocal4PCMFoam::robinLocal4PCMFoam
(
    const robinLocal4PCMFoam& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    kappaName_(ptf.kappaName_),
    qo_(ptf.qo_),
    ho_(ptf.ho_),
    To_(ptf.To_),
    dir_(ptf.dir_),
    p1_(ptf.p1_),
    p2_(ptf.p2_)
{}


Foam::robinLocal4PCMFoam::robinLocal4PCMFoam
(
    const robinLocal4PCMFoam& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    kappaName_(ptf.kappaName_),
    qo_(ptf.qo_),
    ho_(ptf.ho_),
    To_(ptf.To_),
    dir_(ptf.dir_),
    p1_(ptf.p1_),
    p2_(ptf.p2_)
{}


Foam::robinLocal4PCMFoam::robinLocal4PCMFoam
(
    const robinLocal4PCMFoam& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    kappaName_(ptf.kappaName_),
    qo_(ptf.qo_),
    ho_(ptf.ho_),
    To_(ptf.To_),
    dir_(ptf.dir_),
    p1_(ptf.p1_),
    p2_(ptf.p2_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::robinLocal4PCMFoam::qoFunction(scalar t)
{
    int n_ = 0;
    scalar f_ = 0.0;

    while( n_ < qo_.size() )
    {
        f_ += qo_[n_]*pow(t, n_);
        n_ += 1;
    }

    return f_;
}


Foam::scalar Foam::robinLocal4PCMFoam::hoFunction(scalar t)
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


Foam::scalar Foam::robinLocal4PCMFoam::ToFunction(scalar t)
{
    int n_ = 0;
    scalar f_ = 0.0;

    while( n_ < To_.size() )
    {
        f_ += To_[n_]*pow(t, n_);
        n_ += 1;
    }

    return f_;
}


void Foam::robinLocal4PCMFoam::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //- get the currernt time and time step
    scalar t_ = this->db().time().value();


    const fvPatchScalarField& Tp = lookupPatchField<volScalarField, scalar>(TName_);
    const fvPatchScalarField& kappap = lookupPatchField<volScalarField, scalar>(kappaName_);
    scalarField Cf_ = patch().Cf() & dir_;

    forAll(Cf_, i)
    {
        if ( p1_ <= Cf_[i] && Cf_[i] < p2_)
        {
            gradient()[i] = ( qoFunction(t_) - hoFunction(t_) * (Tp[i] - ToFunction(t_)) )/kappap[i];
        }

        else {gradient()[i] = 0.0;}
    }

    //gradient() = ( qoFunction(t_) - hoFunction(t_) * (Tp - ToFunction(t_)) )/kappap;
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::robinLocal4PCMFoam::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("Kappa") << kappaName_ << token::END_STATEMENT << nl;
    os.writeKeyword("qo") << qo_ << token::END_STATEMENT << nl;
    os.writeKeyword("ho") << ho_ << token::END_STATEMENT << nl;
    os.writeKeyword("To") << To_ << token::END_STATEMENT << nl;
    os.writeKeyword("dir") << dir_ << token::END_STATEMENT << nl;
    os.writeKeyword("p1") << p1_ << token::END_STATEMENT << nl;
    os.writeKeyword("p2") << p2_ << token::END_STATEMENT << nl;

    writeEntry("value", os);
    gradient().writeEntry("gradient", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        robinLocal4PCMFoam
    );
}

// ************************************************************************* //
