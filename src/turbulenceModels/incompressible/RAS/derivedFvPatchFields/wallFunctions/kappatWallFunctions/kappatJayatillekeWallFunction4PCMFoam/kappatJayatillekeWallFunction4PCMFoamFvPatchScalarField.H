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

Class
    Foam::incompressible::RASModels::kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField

Description
    Boundary condition for turbulent thermal diffusivity when using wall
    functions
    - replicates OpenFOAM v1.5 (and earlier) behaviour
    - Turbulent Prandtl number defaults to 0.85 if unspecified

SourceFiles
    kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField_H
#define kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
            Class kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Turbulent Prandtl number
        scalar Prt_;

        //- Cmu coefficient
        scalar Cmu_;

        //- Von Karman constant
        scalar kappa_;

        //- E coefficient
        scalar E_;


        // Solution parameters

            static scalar tolerance_;
            static label maxIters_;

    // Private member functions

        //- Check the type of the patch
        void checkType();



protected:

    // Protected Member Functions

        //- Write local wall function variables
        virtual void writeLocalEntries(Ostream&) const;

        //- `P' function
        scalar Psmooth(const scalar Prat) const;

        //- Calculate y+ at the edge of the thermal laminar sublayer
        scalar yPlusTherm
        (
            const scalar P,
            const scalar Prat
        ) const;


public:

    //- Runtime type information
    TypeName("kappatJayatillekeWallFunction4PCMFoam");


    // Constructors

        //- Construct from patch and internal field
        kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField
        //  onto a new patch
        kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField
        (
            const kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField
        (
            const kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField
        (
            const kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new kappatJayatillekeWallFunction4PCMFoamFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
