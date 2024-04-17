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

//------------------------------- PCM4Foam project -------------------------------//
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
m4_define(startAngleOffset, rad(45.0))
m4_define(rOut, calc(tan(60.0*pi/180.0)))

//variables
m4_define(R1, 0.5)
m4_define(R2, 1.0)
m4_define(length, 0.1)
m4_define(numberOfCell_x, 40)
m4_define(numberOfCell_y, 40)
m4_define(numberOfCell_z, 1)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
convertToMeters 1.0;

vertices
(
    //-Back
    (    R1    0    calc(-1*length/2)    ) vlabel(B0)
    (    calc(R1*cos(rad(45.0)))    calc(R1*sin(rad(45.0)))    calc(-1*length/2)    ) vlabel(B1)
    (    0    R1    calc(-1*length/2)    ) vlabel(B2)
    (    calc(R1*cos(rad(135.0)))    calc(R1*sin(rad(135.0)))    calc(-1*length/2) ) vlabel(B3)
    (    calc(-1*R1)    0 calc(-1*length/2)    ) vlabel(B4)
    (    R2    0    calc(-1*length/2)    ) vlabel(B5)
    (    calc(R2*cos(rad(45.0)))    calc(R2*sin(rad(45.0)))    calc(-1*length/2)    ) vlabel(B6)
    (    0    R2    calc(-1*length/2)    ) vlabel(B7)
    (    calc(R2*cos(rad(135.0)))    calc(R2*sin(rad(135.0)))    calc(-1*length/2)    ) vlabel(B8)
    (    calc(-1*R2)    0    calc(-1*length/2)    ) vlabel(B9)

    //-Front
    (    R1    0    calc(length/2)    ) vlabel(F0)
    (    calc(R1*cos(rad(45.0)))    calc(R1*sin(rad(45.0)))    calc(length/2)    ) vlabel(F1)
    (    0    R1    calc(length/2)    ) vlabel(F2)
    (    calc(R1*cos(rad(135.0)))    calc(R1*sin(rad(135.0)))    calc(length/2)    ) vlabel(F3)
    (    calc(-1*R1)    0    calc(length/2)    ) vlabel(F4)
    (    R2    0    calc(length/2)    ) vlabel(F5)
    (    calc(R2*cos(rad(45.0)))    calc(R2*sin(rad(45.0)))    calc(length/2)    ) vlabel(F6)
    (    0    R2    calc(length/2)    ) vlabel(F7)
    (    calc(R2*cos(rad(135.0)))    calc(R2*sin(rad(135.0)))    calc(length/2)    ) vlabel(F8)
    (    calc(-1*R2)    0    calc(length/2)    ) vlabel(F9)
);

blocks
(
    hex (B0 B5 B6 B1 F0 F5 F6 F1) (numberOfCell_x numberOfCell_y numberOfCell_z) simpleGrading (1 1 1)
    hex (B1 B6 B7 B2 F1 F6 F7 F2) (numberOfCell_x numberOfCell_y numberOfCell_z) simpleGrading (1 1 1)
    hex (B2 B7 B8 B3 F2 F7 F8 F3) (numberOfCell_x numberOfCell_y numberOfCell_z) simpleGrading (1 1 1)
    hex (B3 B8 B9 B4 F3 F8 F9 F4) (numberOfCell_x numberOfCell_y numberOfCell_z) simpleGrading (1 1 1)
);

edges
(
    //-Back
    arc B0 B1 (    calc(R1*cos(rad(22.5)))    calc(R1*sin(rad(22.5)))    calc(-1*length/2)    )
    arc B1 B2 (    calc(R1*cos(rad(67.5)))    calc(R1*sin(rad(67.5)))    calc(-1*length/2)    )
    arc B2 B3 (    calc(R1*cos(rad(112.5)))    calc(R1*sin(rad(112.5)))    calc(-1*length/2)    )
    arc B3 B4 (    calc(R1*cos(rad(157.5)))    calc(R1*sin(rad(157.5)))    calc(-1*length/2)    )
    arc B5 B6 (    calc(R2*cos(rad(22.5)))    calc(R2*sin(rad(22.5)))    calc(-1*length/2)    )
    arc B6 B7 (    calc(R2*cos(rad(67.5)))    calc(R2*sin(rad(67.5)))    calc(-1*length/2)    )
    arc B7 B8 (    calc(R2*cos(rad(112.5)))    calc(R2*sin(rad(112.5)))    calc(-1*length/2)    )
    arc B8 B9 (    calc(R2*cos(rad(157.5)))    calc(R2*sin(rad(157.5)))    calc(-1*length/2)    )
    //-Front
    arc F0 F1 (    calc(R1*cos(rad(22.5)))    calc(R1*sin(rad(22.5)))    calc(length/2)    )
    arc F1 F2 (    calc(R1*cos(rad(67.5)))    calc(R1*sin(rad(67.5)))    calc(length/2)    )
    arc F2 F3 (    calc(R1*cos(rad(112.5)))    calc(R1*sin(rad(112.5)))    calc(length/2)    )
    arc F3 F4 (    calc(R1*cos(rad(157.5)))    calc(R1*sin(rad(157.5)))    calc(length/2)    )
    arc F5 F6 (    calc(R2*cos(rad(22.5)))    calc(R2*sin(rad(22.5)))    calc(length/2)    )
    arc F6 F7 (    calc(R2*cos(rad(67.5)))    calc(R2*sin(rad(67.5)))    calc(length/2)    )
    arc F7 F8 (    calc(R2*cos(rad(112.5)))    calc(R2*sin(rad(112.5)))    calc(length/2)    )
    arc F8 F9 (    calc(R2*cos(rad(157.5)))    calc(R2*sin(rad(157.5)))    calc(length/2)    )
);

boundary
(
    innerWall
    {
        type wall;
        faces
        (
            (B1 B0 F0 F1)
            (B2 B1 F1 F2)
            (B3 B2 F2 F3)
            (B4 B3 F3 F4)
        );
    }

    outerWall
    {
        type wall;
        faces
        (
            (B5 B6 F6 F5)
            (B6 B7 F7 F6)
            (B7 B8 F8 F7)
            (B8 B9 F9 F8)
        );
    }

    horizontalWalls
    {
        type wall;
        faces
        (
            (B0 B5 F5 F0)
            (B9 B4 F4 F9)
        );
    }

    frontAndBack
    {
        type empty;
        faces
        (
            (B0 B1 B6 B5)
            (B1 B2 B7 B6)
            (B2 B3 B8 B7)
            (B3 B4 B9 B8)
            (F0 F1 F6 F5)
            (F1 F2 F7 F6)
            (F2 F3 F8 F7)
            (F3 F4 F9 F8)
        );
    }
);

mergePatchPairs
(
);

// ************************************************************************* //
