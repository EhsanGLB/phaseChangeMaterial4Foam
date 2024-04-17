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
    Foam::buongiorno4PCMFoam

Description

    @verbatim
        wall
        {
            type		buongiorno4PCMFoam;
            alpha		alphaEPCM;
            jo			(0.0);
            ho			(0.0);
            alphao		(0.0);
            DB			1e-13;
            DT			1e-13;
        }
    @endverbatim

SourceFiles
    buongiorno4PCMFoam.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef buongiorno4PCMFoams_H
#define buongiorno4PCMFoams_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class nonUniformBuongiornoFvPatch Declaration
\*---------------------------------------------------------------------------*/

class buongiorno4PCMFoam
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Name of the temperature field
        word TName_;

        //- Name of the alpha field
        word alphaName_;

        //- jo
        scalarList jo_;

        //- convective heat transfer
        scalarList ho_;

        //- Out alpha
        scalarList alphao_;

        //- Brownian diffusivity 
        scalar DB_;

        //- Thermophoresis diffusivity
        scalar DT_;


public:

    //- Runtime type information
    TypeName("buongiorno4PCMFoam");


    // Construcalphaors

        //- Construct from patch and internal field
        buongiorno4PCMFoam
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        buongiorno4PCMFoam
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given buongiorno4PCMFoam onalphao
        //  a new patch
        buongiorno4PCMFoam
        (
            const buongiorno4PCMFoam&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        buongiorno4PCMFoam
        (
            const buongiorno4PCMFoam&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new buongiorno4PCMFoam(*this)
            );
        }

        //- Construct as copy setting internal field reference
        buongiorno4PCMFoam
        (
            const buongiorno4PCMFoam&,
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
                new buongiorno4PCMFoam(*this, iF)
            );
        }


    // Member functions

        //- Calculate the jo(time)
        virtual scalar joFunction(scalar time);

        //- Calculate the ho(time)
        virtual scalar hoFunction(scalar time);

        //- Calculate the alphao(time)
        virtual scalar alphaoFunction(scalar time);

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
