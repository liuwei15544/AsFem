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
// Created by walkandthinker on 15.08.18.
// Generate 3D mesh for AsFem

#include "Mesh/Mesh.h"

void Mesh::Create3DMesh()
{
    MeshCreated=false;
    BCMeshCreated=false;

    double dx,dy;
    int e,i,j,k;
    int i1,i2,i3,i4,i5,i6,i7,i8,i9;


    if(MeshType=="quad4")
    {
        dx=(Xmax-Xmin)/Nx;
        dy=(Ymax-Ymin)/Ny;

        nElmts=Nx*Ny;

        nNodes=(Nx+1)*(Ny+1);
        nNodesPerElmt=4;
        nNodesPerBCElmt=2;

        NodeCoords.resize(nNodes*4,0.0);
        Conn.resize(nElmts*nNodesPerElmt,0);

        for(j=1;j<=Ny+1;++j)
        {
            for(i=1;i<=Nx+1;++i)
            {
                k=(j-1)*(Nx+1)+i;
                NodeCoords[(k-1)*4  ]=1.0;
                NodeCoords[(k-1)*4+1]=Xmin+(i-1)*dx;
                NodeCoords[(k-1)*4+2]=Ymin+(j-1)*dy;
                NodeCoords[(k-1)*4+3]=0.0;
            }
        }

        // Create Connectivity matrix
        for(j=1;j<=Ny;j++)
        {
            for(i=1;i<=Nx;i++)
            {
                e=(j-1)*Nx+i;
                i1=(j-1)*(Nx+1)+i;
                i2=i1+1;
                i3=i2+Nx+1;
                i4=i3-1;

                Conn[(e-1)*nNodesPerElmt+1-1]=i1;
                Conn[(e-1)*nNodesPerElmt+2-1]=i2;
                Conn[(e-1)*nNodesPerElmt+3-1]=i3;
                Conn[(e-1)*nNodesPerElmt+4-1]=i4;

            }
        }
        VTKCellType=9;
    }
    else if(MeshType=="quad8")
    {
        // for 2D-8 nodes mesh
        dx=(Xmax-Xmin)/(2.0*Nx);
        dy=(Ymax-Ymin)/(2.0*Ny);

        nElmts=Nx*Ny;
        nNodes=(2*Nx+1)*(2*Ny+1)-nElmts;
        nNodesPerElmt=8;
        nNodesPerBCElmt=3;

        NodeCoords.resize(nNodes*4,0.0);
        Conn.resize(nElmts*nNodesPerElmt,0);

        for(j=1;j<=Ny;j++)
        {
            // for bottom line of each element
            for(i=1;i<=2*Nx+1;i++)
            {
                k=(j-1)*(2*Nx+1+Nx+1)+i;

                NodeCoords[(k-1)*4  ]=1.0;
                NodeCoords[(k-1)*4+1]=Xmin+(i-1)*dx;
                NodeCoords[(k-1)*4+2]=Ymin+(j-1)*2*dy;
                NodeCoords[(k-1)*4+3]=0.0;
            }
            // for middle line of each element
            for(i=1;i<=Nx+1;i++)
            {
                k=(j-1)*(2*Nx+1+Nx+1)+2*Nx+1+i;

                NodeCoords[(k-1)*4  ]=1.0;
                NodeCoords[(k-1)*4+1]=Xmin+(i-1)*2*dx;
                NodeCoords[(k-1)*4+2]=Ymin+(j-1)*2*dy+dy;
                NodeCoords[(k-1)*4+3]=0.0;

            }
        }
        // for the last top line
        j=Ny+1;
        for(i=1;i<=2*Nx+1;i++)
        {
            k=(j-1)*(2*Nx+1+Nx+1)+i;

            NodeCoords[(k-1)*4  ]=1.0;
            NodeCoords[(k-1)*4+1]=Xmin+(i-1)*dx;
            NodeCoords[(k-1)*4+2]=Ymin+(j-1)*2*dy;
            NodeCoords[(k-1)*4+3]=0.0;
        }

        // Create Connectivity matrix
        for(j=1;j<=Ny;j++)
        {
            for(i=1;i<=Nx;i++)
            {
                e=(j-1)*Nx+i;
                i1=(j-1)*(2*Nx+1+Nx+1)+2*i-1;
                i2=i1+2;
                i3=i2+(2*Nx+1+Nx+1);
                i4=i3-2;

                i5=i1+1;
                i6=i2+(2*Nx+1)-i;
                i7=i3-1;
                i8=i1+(2*Nx+1)-(i-1);

                Conn[(e-1)*nNodesPerElmt+1-1]=i1;
                Conn[(e-1)*nNodesPerElmt+2-1]=i2;
                Conn[(e-1)*nNodesPerElmt+3-1]=i3;
                Conn[(e-1)*nNodesPerElmt+4-1]=i4;
                Conn[(e-1)*nNodesPerElmt+5-1]=i5;
                Conn[(e-1)*nNodesPerElmt+6-1]=i6;
                Conn[(e-1)*nNodesPerElmt+7-1]=i7;
                Conn[(e-1)*nNodesPerElmt+8-1]=i8;
            }
        }

        VTKCellType=23;
    }
    else if(MeshType=="quad9")
    {
        dx=(Xmax-Xmin)/(2.0*Nx);
        dy=(Ymax-Ymin)/(2.0*Ny);

        nElmts=Nx*Ny;
        nNodes=(2*Nx+1)*(2*Ny+1);
        nNodesPerElmt=9;
        nNodesPerBCElmt=3;

        NodeCoords.resize(nNodes*4,0.0);
        Conn.resize(nElmts*nNodesPerElmt,0);

        for(j=1;j<=2*Ny+1;j++)
        {
            for(i=1;i<=2*Nx+1;i++)
            {
                k=(j-1)*(2*Nx+1)+i;

                NodeCoords[(k-1)*4  ]=1.0;
                NodeCoords[(k-1)*4+1]=Xmin+(i-1)*dx;
                NodeCoords[(k-1)*4+2]=Ymin+(j-1)*dy;
                NodeCoords[(k-1)*4+3]=0.0;

            }
        }
        // Create Connectivity matrix
        for(j=1;j<=Ny;j++)
        {
            for(i=1;i<=Nx;i++)
            {
                e=(j-1)*Nx+i;
                i1=(j-1)*2*(2*Nx+1)+2*i-1;
                i2=i1+2;
                i3=i2+2*(2*Nx+1);
                i4=i3-2;

                i5=i1+1;
                i6=i2+(2*Nx+1);
                i7=i3-1;
                i8=i1+(2*Nx+1);
                i9=i8+1;

                Conn[(e-1)*nNodesPerElmt+1-1]=i1;
                Conn[(e-1)*nNodesPerElmt+2-1]=i2;
                Conn[(e-1)*nNodesPerElmt+3-1]=i3;
                Conn[(e-1)*nNodesPerElmt+4-1]=i4;
                Conn[(e-1)*nNodesPerElmt+5-1]=i5;
                Conn[(e-1)*nNodesPerElmt+6-1]=i6;
                Conn[(e-1)*nNodesPerElmt+7-1]=i7;
                Conn[(e-1)*nNodesPerElmt+8-1]=i8;
                Conn[(e-1)*nNodesPerElmt+9-1]=i9;

            }
        }

        VTKCellType=28;
    }
}

