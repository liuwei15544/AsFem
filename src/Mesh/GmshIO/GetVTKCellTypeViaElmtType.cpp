//*****************************************************
//*** AsFem: a simple finite element method program ***
//*** Copyright (C) 2018 walkandthinker             ***
//*** Contact: walkandthinker@gmail.com             ***
//*****************************************************
//* This file is part of the ASFEM framework
//* All rights reserved, see COPYRIGHT for full restrictions
//* Licensed under GPL 3.0, please see LICENSE for details
//******************************************************
//
// Created by walkandthinker on 21.04.19.
// get the vtk cell type via gmsh's elmt type

#include "Mesh/GmshIO.h"

int GmshIO::GetVTKCellTypeViaElmtType(int elmttype) const
{
    // All the types are list in:
    // https://vtk.org/doc/release/5.6/html/a02550.html#b1d6fd1f3177b8a2a32bb018807151f873c553260f450a4c939475b66d05daa4

    switch (elmttype)
    {
        case 1:
            // 2-node line.
            // 1--2
            return 3;
        case 2:
            // 3-node triangle.
            return 5;
        case 3:
            // 4-node quadrangle.
            return 9;
        case 4:
            // 4-node tetrahedron.
            return 10;
        case 5:
            // 8-node hexahedron.
            return 12;
        case 7:
            // 5-node pyramid
            // https://vtk.org/doc/nightly/html/classvtkPyramid.html#details
            return 14;
        case 8:
            // 3-node second order line (2 nodes associated with the vertices and
            // 1 with the edge).
            return 21;
        case 9:
            // 6-node second order triangle (3 nodes associated with the vertices
            // and 3 with the edges).
            // https://vtk.org/doc/nightly/html/classvtkQuadraticTriangle.html
            return 22;
        case 10:
            // 9-node second order quadrangle (4 nodes associated with the ver-
            // tices, 4 with the edges and 1 with the face).
            return 28;
        case 11:
            // 10-node second order tetrahedron (4 nodes associated with the ver-
            // tices and 6 with the edges).
            // https://vtk.org/doc/nightly/html/classvtkQuadraticTetra.html
            return 24;
        case 12:
            // 27-node second order hexahedron (8 nodes associated with the ver-
            // tices, 12 with the edges, 6 with the faces and 1 with the volume).
            // https://vtk.org/doc/nightly/html/classvtkTriQuadraticHexahedron.html
            return 29;
        case 15:
            // 1-node point.
            return 1;
        case 16:
            // 8-node second order quadrangle (4 nodes associated with the ver-
            // tices and 4 with the edges).
            return 23;
        case 17:
            // 20-node second order hexahedron (8 nodes associated with the ver-
            // tices and 12 with the edges).
            return 25;
        case 19:
            // 13-node second order pyramid (5 nodes associated with the vertices
            // and 8 with the edges).
            // https://vtk.org/doc/nightly/html/classvtkQuadraticPyramid.html#details
            return 27;
        case 26:
            // 4-node third order edge (2 nodes associated with the vertices, 2
            // internal to the edge)
            return 4;
        case 27:
            // 5-node fourth order edge (2 nodes associated with the vertices, 3
            // internal to the edge)
            return 4;
        case 28:
            // 6-node fifth order edge (2 nodes associated with the vertices, 4
            // internal to the edge)
            return 4;
        default:
            Msg_Gmsh_UnsupportElmtType(elmttype);
            Msg_ExitProgram();
    }
}

