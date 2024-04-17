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

#include "parabolicRadialVelocity4PCMFoam.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

parabolicRadialVelocity4PCMFoam::parabolicRadialVelocity4PCMFoam
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    Vm_(0),
    R_(0.0),
    longDir_(0, 0, 0)
{}


parabolicRadialVelocity4PCMFoam::parabolicRadialVelocity4PCMFoam
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    Vm_(),
    R_(readScalar(dict.lookup("R"))),
    longDir_(dict.lookup("longDir"))
{
    Istream& isVm_ = dict.lookup("Vm");
    isVm_.format(IOstream::ASCII);
    isVm_ >> Vm_;

    longDir_ /= mag(longDir_);

    evaluate();
}


parabolicRadialVelocity4PCMFoam::parabolicRadialVelocity4PCMFoam
(
    const parabolicRadialVelocity4PCMFoam& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    Vm_(ptf.Vm_),
    R_(ptf.R_),
    longDir_(ptf.longDir_)
{}


parabolicRadialVelocity4PCMFoam::parabolicRadialVelocity4PCMFoam
(
    const parabolicRadialVelocity4PCMFoam& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
    Vm_(fcvpvf.Vm_),
    R_(fcvpvf.R_),
    longDir_(fcvpvf.longDir_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::parabolicRadialVelocity4PCMFoam::VmFunction(scalar t)
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


void parabolicRadialVelocity4PCMFoam::updateCoeffs()
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

    vectorField::operator=(longDir_*2.0*VmFunction(t_)*(1.0 - pow(rp_/R_, 2.0)));
    fixedValueFvPatchVectorField::updateCoeffs();
}


// Write
void parabolicRadialVelocity4PCMFoam::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("Vm") << Vm_ << token::END_STATEMENT << nl;
    os.writeKeyword("longDir") << longDir_ << token::END_STATEMENT << nl;
    os.writeKeyword("R") << R_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, parabolicRadialVelocity4PCMFoam);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
