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
    Foam::robinLocal4PCMFoam

Description

    @verbatim
        wall
        {
            type		robinLocal4PCMFoam;
            Kappa		kappa;
            qo			(q0 q1 q2);
            ho			(h0 h1 h2);
            To			(T0 T1 T2);
            dir			(1 0 0);
            p1			0;
            p2			1;
        }
    @endverbatim

SourceFiles
    robinLocal4PCMFoam.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef robinLocal4PCMFoam_H
#define robinLocal4PCMFoam_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class nonUniformRobinFvPatch Declaration
\*---------------------------------------------------------------------------*/

class robinLocal4PCMFoam
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Name of the temperature field
        word TName_;

        //- Name of the thermal conductivity field
        word kappaName_;

        //- Flux depends on time
        scalarList qo_;

        //- convective heat transfer
        scalarList ho_;

        //- Out temperature
        scalarList To_;

        //- Direction
        vector dir_;

        //- p1
        scalar p1_;

        //- p2
        scalar p2_;


public:

    //- Runtime type information
    TypeName("robinLocal4PCMFoam");


    // Constructors

        //- Construct from patch and internal field
        robinLocal4PCMFoam
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        robinLocal4PCMFoam
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given robinLocal4PCMFoam onto
        //  a new patch
        robinLocal4PCMFoam
        (
            const robinLocal4PCMFoam&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        robinLocal4PCMFoam
        (
            const robinLocal4PCMFoam&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new robinLocal4PCMFoam(*this)
            );
        }

        //- Construct as copy setting internal field reference
        robinLocal4PCMFoam
        (
            const robinLocal4PCMFoam&,
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
                new robinLocal4PCMFoam(*this, iF)
            );
        }


    // Member functions

        //- Calculate the qo(time)
        virtual scalar qoFunction(scalar time);

        //- Calculate the ho(time)
        virtual scalar hoFunction(scalar time);

        //- Calculate the To(time)
        virtual scalar ToFunction(scalar time);

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
