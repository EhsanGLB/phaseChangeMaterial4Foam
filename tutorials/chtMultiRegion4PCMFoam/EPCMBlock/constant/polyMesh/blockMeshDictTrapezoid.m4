/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | foam-extend: Open Source CFD                    |
|  \\    /   O peration     | Version:     4.1                                |
|   \\  /    A nd           | Web:         http://www.foam-extend.org         |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//------------------------------- nanoFluid4Foam project -------------------------------//
//Author
    //Ehsan Golab, SUT. All rights reserved.
    //Ehsan1996Golab@gmail.com

//--------------------------------------------------------------------------------------//

//Run using:
//m4 -P blockMeshDict.m4 > blockMeshDict

//m4 definitions:
m4_changecom(//)m4_changequote([,])
m4_define(calc, [m4_esyscmd(perl -e 'use Math::Trig; printf ($1)')])

m4_define(VCOUNT, 0)
m4_define(vlabel, [[// ]Vertex $1 = VCOUNT m4_define($1, VCOUNT)m4_define([VCOUNT], m4_incr(VCOUNT))])

m4_define(pi, 3.14159265358979323844)
m4_define(rad, [calc($1*pi/180.0)])

//variables
m4_define(widthT, 200)
m4_define(widthB, 200)
m4_define(height, 500)
m4_define(length, 150)
m4_define(numberOfCell_x, 80)
m4_define(numberOfCell_y, 20)
m4_define(numberOfCell_z, 60)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
convertToMeters 0.001;

vertices
(
    //-Back
    (    0.0    0.0    calc(-1*length/2)    ) vlabel(B0)
    (    widthB    0.0    calc(-1*length/2)    ) vlabel(B1)
    (    widthT    height    calc(-1*length/2)    ) vlabel(B2)
    (    0.0    height    calc(-1*length/2)    ) vlabel(B3)

    //-Front
    (    0.0    0.0    calc(length/2)    ) vlabel(F0)
    (    widthB    0.0    calc(length/2)    ) vlabel(F1)
    (    widthT    height    calc(length/2)    ) vlabel(F2)
    (    0.0    height    calc(length/2)    ) vlabel(F3)
);

blocks
(
    hex (B0 B1 B2 B3 F0 F1 F2 F3) (numberOfCell_x numberOfCell_y numberOfCell_z) simpleGrading (1 1 1)
);

edges
(
);

boundary
(
    top
    {
        type patch;
        faces
        (
            (B3 F3 F2 B2)
        );
    }

    bottom
    {
        type patch;
        faces
        (
            (B0 B1 F1 F0)
        );
    }

    right
    {
        type patch;
        faces
        (
            (B1 B2 F2 F1)
        );
    }

    left
    {
        type patch;
        faces
        (
            (B0 F0 F3 B3)
        );
    }

    front
    {
        type patch;
        faces
        (
            (F0 F1 F2 F3)
        );
    }

    back
    {
        type patch;
        faces
        (
            (B0 B3 B2 B1)
        );
    }
);

mergePatchPairs
(
);

// ************************************************************************* //
