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

//dimensions of internal square [x, y, z]->[W, H, L]
//width of bottom
m4_define(IntSqWB, 4.0)
//width of top
m4_define(IntSqWT, 4.0)
//height
m4_define(IntSqH, 100.0)

// number of cell
m4_define(IntSqWN, 12)
m4_define(IntSqHN, 50)



//dimensions of cylinder
m4_define(cyRB, 5.0)
m4_define(cyRT, 5.0)
m4_define(cyH, 100.0)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
convertToMeters 0.001;

vertices
(
  //-square
    //-Bottom of square
    (    calc(IntSqWB/2)    calc(-1*IntSqH/2)    calc(-1*IntSqWB/2)    ) vlabel(IntSqB0)
    (    calc(IntSqWB/2)    calc(-1*IntSqH/2)    calc(IntSqWB/2)    ) vlabel(IntSqB1)
    (    calc(-1*IntSqWB/2)    calc(-1*IntSqH/2)    calc(IntSqWB/2)    ) vlabel(IntSqB2)
    (    calc(-1*IntSqWB/2)    calc(-1*IntSqH/2)    calc(-1*IntSqWB/2)    ) vlabel(IntSqB3)

    //-Top of square
    (    calc(IntSqWT/2)    calc(IntSqH/2)    calc(-1*IntSqWT/2)    ) vlabel(IntSqT0)
    (    calc(IntSqWT/2)    calc(IntSqH/2)    calc(IntSqWT/2)    ) vlabel(IntSqT1)
    (    calc(-1*IntSqWT/2)    calc(IntSqH/2)    calc(IntSqWT/2)    ) vlabel(IntSqT2)
    (    calc(-1*IntSqWT/2)    calc(IntSqH/2)    calc(-1*IntSqWT/2)    ) vlabel(IntSqT3)

  //-cylinder
    //-North
    //-Bottom of cylinder
    (    calc(cyRB*cos(rad(45.0)))    calc(-1*cyH/2)    calc(-1*cyRB*cos(rad(45.0)))    ) vlabel(cyB0)
    (    calc(cyRB*cos(rad(45.0)))    calc(-1*cyH/2)    calc(cyRB*cos(rad(45.0)))    ) vlabel(cyB1)
    (    calc(-1*cyRB*cos(rad(45.0)))    calc(-1*cyH/2)    calc(cyRB*cos(rad(45.0)))    ) vlabel(cyB2)
    (    calc(-1*cyRB*cos(rad(45.0)))    calc(-1*cyH/2)    calc(-1*cyRB*cos(rad(45.0)))    ) vlabel(cyB3)
    //-Top of cylinder
    (    calc(cyRT*cos(rad(45.0)))    calc(cyH/2)    calc(-1*cyRT*cos(rad(45.0)))    ) vlabel(cyT0)
    (    calc(cyRT*cos(rad(45.0)))    calc(cyH/2)    calc(cyRT*cos(rad(45.0)))    ) vlabel(cyT1)
    (    calc(-1*cyRT*cos(rad(45.0)))    calc(cyH/2)    calc(cyRT*cos(rad(45.0)))    ) vlabel(cyT2)
    (    calc(-1*cyRT*cos(rad(45.0)))    calc(cyH/2)    calc(-1*cyRT*cos(rad(45.0)))    ) vlabel(cyT3)


);

blocks
(
    hex (IntSqB1 IntSqB2 IntSqT2 IntSqT1 IntSqB0 IntSqB3 IntSqT3 IntSqT0) (IntSqWN IntSqHN IntSqWN) simpleGrading (1 1 1)
    hex (IntSqB0 IntSqB3 IntSqT3 IntSqT0 cyB0 cyB3 cyT3 cyT0) (IntSqWN IntSqHN IntSqWN) simpleGrading (1 1 1)
    hex (IntSqB3 IntSqB2 IntSqT2 IntSqT3 cyB3 cyB2 cyT2 cyT3) (IntSqWN IntSqHN IntSqWN) simpleGrading (1 1 1)
    hex (IntSqB2 IntSqB1 IntSqT1 IntSqT2 cyB2 cyB1 cyT1 cyT2) (IntSqWN IntSqHN IntSqWN) simpleGrading (1 1 1)
    hex (IntSqB1 IntSqB0 IntSqT0 IntSqT1 cyB1 cyB0 cyT0 cyT1) (IntSqWN IntSqHN IntSqWN) simpleGrading (1 1 1)
);

edges
(
  //-cylinder
    //-Bottom of cylinder
    arc cyB0 cyB1 (    cyRB    calc(-1*cyH/2)    0    )
    arc cyB1 cyB2 (    0    calc(-1*cyH/2)    cyRB    )
    arc cyB2 cyB3 (    calc(-1*cyRB)    calc(-1*cyH/2)    0    )
    arc cyB3 cyB0 (    0    calc(-1*cyH/2)    calc(-1*cyRB)    )

    //-Top of cylinder
    arc cyT0 cyT1 (    cyRT    calc(cyH/2)    0    )
    arc cyT1 cyT2 (    0    calc(cyH/2)    cyRT    )
    arc cyT2 cyT3 (    calc(-1*cyRT)    calc(cyH/2)    0    )
    arc cyT3 cyT0 (    0    calc(cyH/2)    calc(-1*cyRT)    )
);

boundary
(
    inlet
    {
        type wall;
        faces
        (
            (IntSqB0 IntSqB1 IntSqB2 IntSqB3)
            (IntSqB0 cyB0 cyB1 IntSqB1)
            (IntSqB3 cyB3 cyB0 IntSqB0)
            (IntSqB2 cyB2 cyB3 IntSqB3)
            (IntSqB1 cyB1 cyB2 IntSqB2)
        );
    }

    outlet
    {
        type wall;
        faces
        (
            (IntSqT0 IntSqT1 IntSqT2 IntSqT3)
            (IntSqT0 cyT0 cyT1 IntSqT1)
            (IntSqT3 cyT3 cyT0 IntSqT0)
            (IntSqT2 cyT2 cyT3 IntSqT3)
            (IntSqT1 cyT1 cyT2 IntSqT2)
        );
    }

    cylinderWall
    {
        type wall;
        faces
        (
            (cyB1 cyB0 cyT0 cyT1)
            (cyB0 cyB3 cyT3 cyT0)
            (cyB3 cyB2 cyT2 cyT3)
            (cyB2 cyB1 cyT1 cyT2)
        );
    }
);

mergePatchPairs
(
);

// ************************************************************************* //
