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

#include "powerLawRadialVelocity4PCMFoam.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

powerLawRadialVelocity4PCMFoam::powerLawRadialVelocity4PCMFoam
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    Vm_(0),
    R_(0.0),
    m_(0.0),
    longDir_(0, 0, 0)
{}


powerLawRadialVelocity4PCMFoam::powerLawRadialVelocity4PCMFoam
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    Vm_(),
    R_(readScalar(dict.lookup("R"))),
    m_(readScalar(dict.lookup("m"))),
    longDir_(dict.lookup("longDir"))
{
    Istream& isVm_ = dict.lookup("Vm");
    isVm_.format(IOstream::ASCII);
    isVm_ >> Vm_;

    longDir_ /= mag(longDir_);

    evaluate();
}


powerLawRadialVelocity4PCMFoam::powerLawRadialVelocity4PCMFoam
(
    const powerLawRadialVelocity4PCMFoam& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    Vm_(ptf.Vm_),
    R_(ptf.R_),
    m_(ptf.m_),
    longDir_(ptf.longDir_)
{}


powerLawRadialVelocity4PCMFoam::powerLawRadialVelocity4PCMFoam
(
    const powerLawRadialVelocity4PCMFoam& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
    Vm_(fcvpvf.Vm_),
    R_(fcvpvf.R_),
    m_(fcvpvf.m_),
    longDir_(fcvpvf.longDir_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::powerLawRadialVelocity4PCMFoam::VmFunction(scalar t)
{
    int n_ = 0;
    scalar f_ = 0.0;

    while( n_ < Vm_.size() )
    {
        f_ += Vm_[n_]*pow(t, n_);
        n_ += 1;
    }

    return f_;
}


void powerLawRadialVelocity4PCMFoam::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    //- get the currernt time and time step
    scalar t_ = this->db().time().value();

    boundBox bb_(patch().patch().localPoints(), true);
    vector ctr_ = 0.5*(bb_.max() + bb_.min());
    const vectorField& c_ = patch().Cf();
    scalarField rp_ = mag(c_ - ctr_);
    scalar Cm_ = (1 + 1/m_)*(2 + 1/m_) / 2;

    vectorField::operator=(longDir_*Cm_*VmFunction(t_)*pow( (1.0 - (rp_/R_)), 1/m_));
    fixedValueFvPatchVectorField::updateCoeffs();
}


// Write
void powerLawRadialVelocity4PCMFoam::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("Vm") << Vm_ << token::END_STATEMENT << nl;
    os.writeKeyword("longDir") << longDir_ << token::END_STATEMENT << nl;
    os.writeKeyword("R") << R_ << token::END_STATEMENT << nl;
    os.writeKeyword("m") << m_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, powerLawRadialVelocity4PCMFoam);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
