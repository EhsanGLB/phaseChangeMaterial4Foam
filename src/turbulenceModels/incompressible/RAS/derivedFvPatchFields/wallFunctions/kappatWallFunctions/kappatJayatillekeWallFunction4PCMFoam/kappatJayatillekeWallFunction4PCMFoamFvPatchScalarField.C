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

#include "kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField::tolerance_ = 0.01;
label kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField::maxIters_ = 10;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField::checkType()
{
    if (!patch().isWall())
    {
        FatalErrorInFunction
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


void kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    os.writeKeyword("Prt") << Prt_ << token::END_STATEMENT << nl;
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
}


scalar kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField::Psmooth
(
    const scalar Prat
) const
{
    return 9.24*(pow(Prat, 0.75) - 1.0)*(1.0 + 0.28*exp(-0.007*Prat));
}


scalar kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField::yPlusTherm
(
    const scalar P,
    const scalar Prat
) const
{
    scalar ypt = 11.0;

    for (int i=0; i<maxIters_; i++)
    {
        scalar f = ypt - (log(E_*ypt)/kappa_ + P)/Prat;
        scalar df = 1.0 - 1.0/(ypt*kappa_*Prat);
        scalar yptNew = ypt - f/df;

        if (yptNew < VSMALL)
        {
            return 0;
        }
        else if (mag(yptNew - ypt) < tolerance_)
        {
            return yptNew;
        }
        else
        {
            ypt = yptNew;
        }
     }

    return ypt;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField::kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Prt_(0.85),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8)
{
    checkType();
}


kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField::kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField
(
    const kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    Prt_(ptf.Prt_),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_)
{
    checkType();
}


kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField::kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    Prt_(dict.lookupOrDefault<scalar>("Prt", 0.85)),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8))
{
    checkType();
}


kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField::kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField
(
    const kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField& awfpsf
)
:
    fixedValueFvPatchScalarField(awfpsf),
    Prt_(awfpsf.Prt_),
    Cmu_(awfpsf.Cmu_),
    kappa_(awfpsf.kappa_),
    E_(awfpsf.E_)
{
    checkType();
}


kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField::kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField
(
    const kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField& awfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(awfpsf, iF),
    Prt_(awfpsf.Prt_),
    Cmu_(awfpsf.Cmu_),
    kappa_(awfpsf.kappa_),
    E_(awfpsf.E_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();
    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");
    const scalar Cmu25 = pow(Cmu_, 0.25);
    const scalarField& y = rasModel.y()[patchi];
    const scalarField& nuw = lookupPatchField<volScalarField, scalar>("nu");
    const scalarField& rhow = lookupPatchField<volScalarField, scalar>("rho");
    const scalarField& Cpw = lookupPatchField<volScalarField, scalar>("Cp");
    const volScalarField& k = db().lookupObject<volScalarField>("k");
    const IOdictionary& transportProperties = db().lookupObject<IOdictionary>("transportProperties");
    dimensionedScalar Prd(transportProperties.lookup("Pr"));
    const scalar Pr = Prd.value();

    // Populate boundary values
    scalarField& kappatw = *this;
    forAll(kappatw, facei)
    {
        label faceCellI = patch().faceCells()[facei];

        // y+
        scalar yPlus = Cmu25*sqrt(k[faceCellI])*y[facei]/nuw[facei];

        // Molecular-to-turbulent Prandtl number ratio
        scalar Prat = Pr/Prt_;

        // Thermal sublayer thickness
        scalar P = Psmooth(Prat);
        scalar yPlusTherm = this->yPlusTherm(P, Prat);

        // Update turbulent thermal conductivity
        if (yPlus > yPlusTherm)
        {
            scalar nu = nuw[facei];
            scalar kt = nu*(yPlus/(Prt_*(log(E_*yPlus)/kappa_ + P)) - 1/Pr)*(rhow[facei]*Cpw[facei]);
            kappatw[facei] = max(0.0, kt);
        }
        else
        {
            kappatw[facei] = 0.0;
        }
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
