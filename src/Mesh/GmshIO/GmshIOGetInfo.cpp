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
// Created by walkandthinker on 08.04.19.
// Define the gmsh io class's get functions(read .msh file from gmsh)

#include "Mesh/GmshIO.h"

int GmshIO::GetNodesNumViaElmtType(int elmttype) const
{
    switch (elmttype)
    {
        case 1:
            // 2-node line.
            return 2;
        case 2:
            // 3-node triangle.
            return 3;
        case 3:
            // 4-node quadrangle.
            return 4;
        case 4:
            // 4-node tetrahedron.
            return 4;
        case 5:
            // 8-node hexahedron.
            return 8;
        case 6:
            // 6-node prism
            return 6;
        case 7:
            // 5-node pyramid
            return 5;
        case 8:
            // 3-node second order line (2 nodes associated with the vertices and
            // 1 with the edge).
            return 3;
        case 9:
            // 6-node second order triangle (3 nodes associated with the vertices
            // and 3 with the edges).
            return 6;
        case 10:
            // 9-node second order quadrangle (4 nodes associated with the ver-
            // tices, 4 with the edges and 1 with the face).
            return 9;
        case 11:
            // 10-node second order tetrahedron (4 nodes associated with the ver-
            // tices and 6 with the edges).
            return 10;
        case 12:
            // 27-node second order hexahedron (8 nodes associated with the ver-
            // tices, 12 with the edges, 6 with the faces and 1 with the volume).
            return 27;
        case 13:
            // 18-node second order prism (6 nodes associated with the vertices,
            // 9 with the edges and 3 with the quadrangular faces).
            return 18;
        case 14:
            // 14-node second order pyramid (5 nodes associated with the vertices,
            // 8 with the edges and 1 with the quadrangular face).
            return 14;
        case 15:
            // 1-node point.
            return 1;
        case 16:
            // 8-node second order quadrangle (4 nodes associated with the ver-
            // tices and 4 with the edges).
            return 8;
        case 17:
            // 20-node second order hexahedron (8 nodes associated with the ver-
            // tices and 12 with the edges).
            return 20;
        case 18:
            // 15-node second order prism (6 nodes associated with the vertices
            // and 9 with the edges).
            return 15;
        case 19:
            // 13-node second order pyramid (5 nodes associated with the vertices
            // and 8 with the edges).
            return 13;
        case 20:
            // 9-node third order incomplete triangle (3 nodes associated with the
            // vertices, 6 with the edges)
            return 9;
        case 21:
            // 10-node third order triangle (3 nodes associated with the vertices,
            // 6 with the edges, 1 with the face)
            return 10;
        case 22:
            // 12-node fourth order incomplete triangle (3 nodes associated with
            // the vertices, 9 with the edges)
            return 12;
        case 23:
            // 15-node fourth order triangle (3 nodes associated with the vertices,
            // 9 with the edges, 3 with the face)
            return 15;
        case 24:
            // 15-node fifth order incomplete triangle (3 nodes associated with the
            // vertices, 12 with the edges)
            return 15;
        case 25:
            // 21-node fifth order complete triangle (3 nodes associated with the
            // vertices, 12 with the edges, 6 with the face)
            return 21;
        case 26:
            // 4-node third order edge (2 nodes associated with the vertices, 2
            // internal to the edge)
            return 4;
        case 27:
            // 5-node fourth order edge (2 nodes associated with the vertices, 3
            // internal to the edge)
            return 5;
        case 28:
            // 6-node fifth order edge (2 nodes associated with the vertices, 4
            // internal to the edge)
            return 6;
        case 29:
            // 20-node third order tetrahedron (4 nodes associated with the ver-
            // tices, 12 with the edges, 4 with the faces)
            return 20;
        case 30:
            // 35-node fourth order tetrahedron (4 nodes associated with the ver-
            // tices, 18 with the edges, 12 with the faces, 1 in the volume)
            return 35;
        case 31:
            // 56-node fifth order tetrahedron (4 nodes associated with the ver-
            // tices, 24 with the edges, 24 with the faces, 4 in the volume)
            return 56;
        case 92:
            // 64-node third order hexahedron (8 nodes associated with the ver-
            // tices, 24 with the edges, 24 with the faces, 8 in the volume)
            return 64;
        case 93:
            // 125-node fourth order hexahedron (8 nodes associated with the
            // vertices, 36 with the edges, 54 with the faces, 27 in the volume)
            return 125;
        default:
            Msg_Gmsh_UnsupportElmtType(elmttype);
            Msg_ExitProgram();
    }
}

//********************************************************
int GmshIO::GetElmtDimViaElmtType(int elmttype) const
{
    // here the dim is not the real dimension, is just the physical dim
    // i.e. point   is 0D
    //      line    is 1D
    //      surface is 2D
    //      volume  is 3D
    // but line can be a 3D line(1D line element in 3D space)
    // so here we just return the physical dim, not the real dim!!!
    switch (elmttype)
    {
        case 1:
            // 2-node line.
            return 1;
        case 2:
            // 3-node triangle.
            return 2;
        case 3:
            // 4-node quadrangle.
            return 2;
        case 4:
            // 4-node tetrahedron.
            return 3;
        case 5:
            // 8-node hexahedron.
            return 3;
        case 6:
            // 6-node prism
            return 3;
        case 7:
            // 5-node pyramid
            return 3;
        case 8:
            // 3-node second order line (2 nodes associated with the vertices and
            // 1 with the edge).
            return 1;
        case 9:
            // 6-node second order triangle (3 nodes associated with the vertices
            // and 3 with the edges).
            return 2;
        case 10:
            // 9-node second order quadrangle (4 nodes associated with the ver-
            // tices, 4 with the edges and 1 with the face).
            return 2;
        case 11:
            // 10-node second order tetrahedron (4 nodes associated with the ver-
            // tices and 6 with the edges).
            return 3;
        case 12:
            // 27-node second order hexahedron (8 nodes associated with the ver-
            // tices, 12 with the edges, 6 with the faces and 1 with the volume).
            return 3;
        case 13:
            // 18-node second order prism (6 nodes associated with the vertices,
            // 9 with the edges and 3 with the quadrangular faces).
            return 2;
        case 14:
            // 14-node second order pyramid (5 nodes associated with the vertices,
            // 8 with the edges and 1 with the quadrangular face).
            return 2;
        case 15:
            // 1-node point.
            return 0;
        case 16:
            // 8-node second order quadrangle (4 nodes associated with the ver-
            // tices and 4 with the edges).
            return 2;
        case 17:
            // 20-node second order hexahedron (8 nodes associated with the ver-
            // tices and 12 with the edges).
            return 3;
        case 18:
            // 15-node second order prism (6 nodes associated with the vertices
            // and 9 with the edges).
            return 3;
        case 19:
            // 13-node second order pyramid (5 nodes associated with the vertices
            // and 8 with the edges).
            return 3;
        case 20:
            // 9-node third order incomplete triangle (3 nodes associated with the
            // vertices, 6 with the edges)
            return 2;
        case 21:
            // 10-node third order triangle (3 nodes associated with the vertices,
            // 6 with the edges, 1 with the face)
            return 2;
        case 22:
            // 12-node fourth order incomplete triangle (3 nodes associated with
            // the vertices, 9 with the edges)
            return 2;
        case 23:
            // 15-node fourth order triangle (3 nodes associated with the vertices,
            // 9 with the edges, 3 with the face)
            return 2;
        case 24:
            // 15-node fifth order incomplete triangle (3 nodes associated with the
            // vertices, 12 with the edges)
            return 2;
        case 25:
            // 21-node fifth order complete triangle (3 nodes associated with the
            // vertices, 12 with the edges, 6 with the face)
            return 2;
        case 26:
            // 4-node third order edge (2 nodes associated with the vertices, 2
            // internal to the edge)
            return 1;
        case 27:
            // 5-node fourth order edge (2 nodes associated with the vertices, 3
            // internal to the edge)
            return 1;
        case 28:
            // 6-node fifth order edge (2 nodes associated with the vertices, 4
            // internal to the edge)
            return 1;
        case 29:
            // 20-node third order tetrahedron (4 nodes associated with the ver-
            // tices, 12 with the edges, 4 with the faces)
            return 3;
        case 30:
            // 35-node fourth order tetrahedron (4 nodes associated with the ver-
            // tices, 18 with the edges, 12 with the faces, 1 in the volume)
            return 3;
        case 31:
            // 56-node fifth order tetrahedron (4 nodes associated with the ver-
            // tices, 24 with the edges, 24 with the faces, 4 in the volume)
            return 3;
        case 92:
            // 64-node third order hexahedron (8 nodes associated with the ver-
            // tices, 24 with the edges, 24 with the faces, 8 in the volume)
            return 3;
        case 93:
            // 125-node fourth order hexahedron (8 nodes associated with the
            // vertices, 36 with the edges, 54 with the faces, 27 in the volume)
            return 3;
        default:
            Msg_Gmsh_UnsupportElmtType(elmttype);
            Msg_ExitProgram();
    }
}


//*********************************************************
vector<int> GmshIO::GetNodeOrderViaElmtType(int elmttype) const
{
    // in gmsh, the index starts from zero, but we set it to start from 1 !!!
    vector<int> temp;
    switch (elmttype)
    {
        case 1:
            // 2-node line
            temp=vector<int>{1+0,1+1};
            return temp;
        case 2:
            // 3-node triangle.
            temp=vector<int>{1+0,1+1,2+1};
            return temp;
        case 3:
            // 4-node quadrangle.
            temp=vector<int>{1+0,1+1,2+1};
            return temp;
        case 4:
            // 4-node tetrahedron.
            temp=vector<int>{1+0,1+1,1+2,1+3};
            return temp;
//        case 5:
//            // 8-node hexahedron.
//            return 8;
//        case 6:
//            // 6-node prism
//            return 6;
//        case 7:
//            // 5-node pyramid
//            return 5;
//        case 8:
//            // 3-node second order line (2 nodes associated with the vertices and
//            // 1 with the edge).
//            return 3;
//        case 9:
//            // 6-node second order triangle (3 nodes associated with the vertices
//            // and 3 with the edges).
//            return 6;
//        case 10:
//            // 9-node second order quadrangle (4 nodes associated with the ver-
//            // tices, 4 with the edges and 1 with the face).
//            return 9;
//        case 11:
//            // 10-node second order tetrahedron (4 nodes associated with the ver-
//            // tices and 6 with the edges).
//            return 10;
//        case 12:
//            // 27-node second order hexahedron (8 nodes associated with the ver-
//            // tices, 12 with the edges, 6 with the faces and 1 with the volume).
//            return 27;
//        case 13:
//            // 18-node second order prism (6 nodes associated with the vertices,
//            // 9 with the edges and 3 with the quadrangular faces).
//            return 18;
//        case 14:
//            // 14-node second order pyramid (5 nodes associated with the vertices,
//            // 8 with the edges and 1 with the quadrangular face).
//            return 14;
//        case 15:
//            // 1-node point.
//            return 1;
//        case 16:
//            // 8-node second order quadrangle (4 nodes associated with the ver-
//            // tices and 4 with the edges).
//            return 8;
//        case 17:
//            // 20-node second order hexahedron (8 nodes associated with the ver-
//            // tices and 12 with the edges).
//            return 20;
//        case 18:
//            // 15-node second order prism (6 nodes associated with the vertices
//            // and 9 with the edges).
//            return 15;
//        case 19:
//            // 13-node second order pyramid (5 nodes associated with the vertices
//            // and 8 with the edges).
//            return 13;
//        case 20:
//            // 9-node third order incomplete triangle (3 nodes associated with the
//            // vertices, 6 with the edges)
//            return 9;
//        case 21:
//            // 10-node third order triangle (3 nodes associated with the vertices,
//            // 6 with the edges, 1 with the face)
//            return 10;
//        case 22:
//            // 12-node fourth order incomplete triangle (3 nodes associated with
//            // the vertices, 9 with the edges)
//            return 12;
//        case 23:
//            // 15-node fourth order triangle (3 nodes associated with the vertices,
//            // 9 with the edges, 3 with the face)
//            return 15;
//        case 24:
//            // 15-node fifth order incomplete triangle (3 nodes associated with the
//            // vertices, 12 with the edges)
//            return 15;
//        case 25:
//            // 21-node fifth order complete triangle (3 nodes associated with the
//            // vertices, 12 with the edges, 6 with the face)
//            return 21;
//        case 26:
//            // 4-node third order edge (2 nodes associated with the vertices, 2
//            // internal to the edge)
//            return 4;
//        case 27:
//            // 5-node fourth order edge (2 nodes associated with the vertices, 3
//            // internal to the edge)
//            return 5;
//        case 28:
//            // 6-node fifth order edge (2 nodes associated with the vertices, 4
//            // internal to the edge)
//            return 6;
//        case 29:
//            // 20-node third order tetrahedron (4 nodes associated with the ver-
//            // tices, 12 with the edges, 4 with the faces)
//            return 20;
//        case 30:
//            // 35-node fourth order tetrahedron (4 nodes associated with the ver-
//            // tices, 18 with the edges, 12 with the faces, 1 in the volume)
//            return 35;
//        case 31:
//            // 56-node fifth order tetrahedron (4 nodes associated with the ver-
//            // tices, 24 with the edges, 24 with the faces, 4 in the volume)
//            return 56;
//        case 92:
//            // 64-node third order hexahedron (8 nodes associated with the ver-
//            // tices, 24 with the edges, 24 with the faces, 8 in the volume)
//            return 64;
//        case 93:
//            // 125-node fourth order hexahedron (8 nodes associated with the
//            // vertices, 36 with the edges, 54 with the faces, 27 in the volume)
//            return 125;
        default:
            Msg_Gmsh_UnsupportElmtType(elmttype);
            Msg_ExitProgram();
    }
}