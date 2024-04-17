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
    Foam::timeVaryingRobin4PCMFoam

Description

    @verbatim
        wall
        {
            type                    timeVaryingRobin4PCMFoam;
            Kappa                   kappa;
            qo
            {
                fileName        "$FOAM_CASE/constant/qovsTime";
                outOfBounds     clamp;           // (error|warn|clamp|repeat)
            }
            ho
            {
                fileName        "$FOAM_CASE/constant/hovsTime";
                outOfBounds     clamp;
            }
            To
            {
                fileName        "$FOAM_CASE/constant/TovsTime";
                outOfBounds     clamp;
            }
        }
    @endverbatim

SourceFiles
    timeVaryingRobin4PCMFoam.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingRobin4PCMFoam_H
#define timeVaryingRobin4PCMFoam_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"
#include "Switch.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class nonUniformRobinFvPatch Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingRobin4PCMFoam
:
    public fixedGradientFvPatchScalarField
{
    // Private data

        //- Name of the temperature field
        word TName_;

        //- Name of the thermal conductivity field
        word kappaName_;

        //- Table of time vs heat flux, including the bounding treatment
        interpolationTable<scalar> timeSeriesqo_;

        //- Table of time vs convective heat transfer coefficient, including the bounding treatment
        interpolationTable<scalar> timeSeriesho_;

        //- Table of time vs temperature, including the bounding treatment
        interpolationTable<scalar> timeSeriesTo_;


public:

    //- Runtime type information
    TypeName("timeVaryingRobin4PCMFoam");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingRobin4PCMFoam
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingRobin4PCMFoam
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given timeVaryingRobin4PCMFoam onto
        //  a new patch
        timeVaryingRobin4PCMFoam
        (
            const timeVaryingRobin4PCMFoam&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingRobin4PCMFoam
        (
            const timeVaryingRobin4PCMFoam&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new timeVaryingRobin4PCMFoam(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingRobin4PCMFoam
        (
            const timeVaryingRobin4PCMFoam&,
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
                new timeVaryingRobin4PCMFoam(*this, iF)
            );
        }


    // Member functions

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
