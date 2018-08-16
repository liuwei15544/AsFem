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

    double dx,dy,dz;
    int e,i,j,k,kk;
    int i1,i2,i3,i4,i5,i6,i7,i8,i9;
    int i10,i11,i12,i13,i14,i15,i16,i17,i18,i19;
    int i20,i21,i22,i23,i24,i25,i26,i27;


    if(MeshType=="hex8")
    {
        dx=(Xmax-Xmin)/Nx;
        dy=(Ymax-Ymin)/Ny;
        dz=(Zmax-Zmin)/Nz;

        nElmts=Nx*Ny*Nz;
        nNodes=(Nx+1)*(Ny+1)*(Nz+1);
        nNodesPerElmt=8;
        nNodesPerBCElmt=4;

        NodeCoords.resize(nNodes*4,0.0);
        Conn.resize(nElmts*nNodesPerElmt,0);

        for(k=1;k<=Nz+1;k++)
        {
            for(j=1;j<=Ny+1;++j)
            {
                for(i=1;i<=Nx+1;++i)
                {
                    kk=(j-1)*(Nx+1)+i+(k-1)*(Nx+1)*(Ny+1);
                    NodeCoords[(kk-1)*4  ]=1.0;
                    NodeCoords[(kk-1)*4+1]=Xmin+(i-1)*dx;
                    NodeCoords[(kk-1)*4+2]=Ymin+(j-1)*dy;
                    NodeCoords[(kk-1)*4+3]=Zmin+(k-1)*dz;
                }
            }
        }

        // Create Connectivity matrix
        kk=0;
        for(k=1;k<=Nz;k++)
        {
            for(j=1;j<=Ny;j++)
            {
                for(i=1;i<=Nx;i++)
                {
                    e=(j-1)*Nx+i+(k-1)*Nx*Ny;
                    i1=(j-1)*(Nx+1)+i+(k-1)*(Nx+1)*(Ny+1);
                    i2=i1+1;
                    i3=i2+Nx+1;
                    i4=i3-1;
                    i5=i1+(Nx+1)*(Ny+1);
                    i6=i2+(Nx+1)*(Ny+1);
                    i7=i3+(Nx+1)*(Ny+1);
                    i8=i4+(Nx+1)*(Ny+1);

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
        }
        VTKCellType=12;
    }
    else if(MeshType=="hex20")
    {
        // for 2D-8 nodes mesh
        dx=(Xmax-Xmin)/(2.0*Nx);
        dy=(Ymax-Ymin)/(2.0*Ny);
        dz=(Zmax-Zmin)/(2.0*Nz);

        nElmts=Nx*Ny*Nz;
        int nLayer1Nodes=(2*Nx+1)*(2*Ny+1)-Nx*Ny;// for norm layer
        int nLayer2Nodes=(Nx+1)*(Ny+1);          // for middle layer

        nNodes=nLayer1Nodes*(Nz+1)+nLayer2Nodes*Nz;
        nNodesPerElmt=20;
        nNodesPerBCElmt=8;

        NodeCoords.resize(nNodes*4,0.0);
        Conn.resize(nElmts*nNodesPerElmt,0);

        for(k=1;k<=Nz;++k)
        {
            // First for normal layer
            for(j=1;j<=Ny;++j)
            {
                // for bottom line of each element
                for(i=1;i<=2*Nx+1;i++)
                {
                    kk=(j-1)*(2*Nx+1+Nx+1)+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                    NodeCoords[4*(kk-1)+0]=1.0;
                    NodeCoords[4*(kk-1)+1]=Xmin+(i-1)*dx;
                    NodeCoords[4*(kk-1)+2]=Ymin+(j-1)*2*dy;
                    NodeCoords[4*(kk-1)+3]=Zmin+(k-1)*2*dz;
                }
                // for middle line of each element
                for(i=1;i<=Nx+1;++i)
                {
                    kk=(j-1)*(2*Nx+1+Nx+1)+2*Nx+1+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                    NodeCoords[4*(kk-1)+0]=1.0;
                    NodeCoords[4*(kk-1)+1]=Xmin+(i-1)*2*dx;
                    NodeCoords[4*(kk-1)+2]=Ymin+(j-1)*2*dy+dy;
                    NodeCoords[4*(kk-1)+3]=Zmin+(k-1)*2*dz;
                }
            }
            // for top line
            j=Ny+1;
            for(i=1;i<=2*Nx+1;i++)
            {
                kk=(j-1)*(2*Nx+1+Nx+1)+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                NodeCoords[4*(kk-1)+0]=1.0;
                NodeCoords[4*(kk-1)+1]=domain.BottomBottomLeftPoint().X()+(i-1)*dx;
                NodeCoords[4*(kk-1)+2]=domain.BottomBottomLeftPoint().Y()+(j-1)*2*dy;
                NodeCoords(3,kk-1)=domain.BottomBottomLeftPoint().Z()+(k-1)*2*dz;
            }
            // Then for middle type layer
            for(j=1;j<=ny+1;++j)
            {
                for(i=1;i<=nx+1;++i)
                {
                    kk=(j-1)*(nx+1)+i+(k-1)*(nLayer1Nodes+nLayer2Nodes)+nLayer1Nodes;
                    NodeCoords(0,kk-1)=1.0;
                    NodeCoords(1,kk-1)=domain.BottomBottomLeftPoint().X()+(i-1)*2*dx;
                    NodeCoords(2,kk-1)=domain.BottomBottomLeftPoint().Y()+(j-1)*2*dy;
                    NodeCoords(3,kk-1)=domain.BottomBottomLeftPoint().Z()+(k-1)*2*dz+dz;
                }
            }
        }
        // for the last top layer
        k=nz+1;
        for(j=1;j<=ny;++j)
        {
            // for bottom line of each element
            for(i=1;i<=2*nx+1;i++)
            {
                kk=(j-1)*(2*nx+1+nx+1)+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                NodeCoords(0,kk-1)=1.0;
                NodeCoords(1,kk-1)=domain.BottomBottomLeftPoint().X()+(i-1)*dx;
                NodeCoords(2,kk-1)=domain.BottomBottomLeftPoint().Y()+(j-1)*2*dy;
                NodeCoords(3,kk-1)=domain.BottomBottomLeftPoint().Z()+(k-1)*2*dz;
            }
            // for middle line of each element
            for(i=1;i<=nx+1;++i)
            {
                kk=(j-1)*(2*nx+1+nx+1)+2*nx+1+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                NodeCoords(0,kk-1)=1.0;
                NodeCoords(1,kk-1)=domain.BottomBottomLeftPoint().X()+(i-1)*2*dx;
                NodeCoords(2,kk-1)=domain.BottomBottomLeftPoint().Y()+(j-1)*2*dy+dy;
                NodeCoords(3,kk-1)=domain.BottomBottomLeftPoint().Z()+(k-1)*2*dz;
            }
        }
        // for top line
        j=ny+1;
        for(i=1;i<=2*nx+1;i++)
        {
            kk=(j-1)*(2*nx+1+nx+1)+i+(k-1)*(nLayer1Nodes+nLayer2Nodes);
            NodeCoords(0,kk-1)=1.0;
            NodeCoords(1,kk-1)=domain.BottomBottomLeftPoint().X()+(i-1)*dx;
            NodeCoords(2,kk-1)=domain.BottomBottomLeftPoint().Y()+(j-1)*2*dy;
            NodeCoords(3,kk-1)=domain.BottomBottomLeftPoint().Z()+(k-1)*2*dz;
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

