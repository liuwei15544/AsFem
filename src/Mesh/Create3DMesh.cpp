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
                NodeCoords[4*(kk-1)+1]=Xmin+(i-1)*dx;
                NodeCoords[4*(kk-1)+2]=Ymin+(j-1)*2*dy;
                NodeCoords[4*(kk-1)+3]=Zmin+(k-1)*2*dz;
            }
            // Then for middle type layer
            for(j=1;j<=Ny+1;++j)
            {
                for(i=1;i<=Nx+1;++i)
                {
                    kk=(j-1)*(Nx+1)+i+(k-1)*(nLayer1Nodes+nLayer2Nodes)+nLayer1Nodes;
                    NodeCoords[4*(kk-1)+0]=1.0;
                    NodeCoords[4*(kk-1)+1]=Xmin+(i-1)*2*dx;
                    NodeCoords[4*(kk-1)+2]=Ymin+(j-1)*2*dy;
                    NodeCoords[4*(kk-1)+3]=Zmin+(k-1)*2*dz+dz;
                }
            }
        }
        // for the last top layer
        k=Nz+1;
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
            NodeCoords[4*(kk-1)+1]=Xmin+(i-1)*dx;
            NodeCoords[4*(kk-1)+2]=Ymin+(j-1)*2*dy;
            NodeCoords[4*(kk-1)+3]=Zmin+(k-1)*2*dz;
        }

        // Create Connectivity matrix
        for(k=1;k<=Nz;++k)
        {
            for(j=1;j<=Ny;++j)
            {
                for(i=1;i<=Nx;++i)
                {
                    e=(j-1)*Nx+i+(k-1)*Nx*Ny;
                    i1=(j-1)*(2*Nx+1+Nx+1)+2*i-1+(k-1)*(nLayer1Nodes+nLayer2Nodes);
                    i2=i1+2;
                    i3=i2+(2*Nx+1+Nx+1);
                    i4=i3-2;

                    i5=i1+nLayer1Nodes+nLayer2Nodes;
                    i6=i2+nLayer1Nodes+nLayer2Nodes;
                    i7=i3+nLayer1Nodes+nLayer2Nodes;
                    i8=i4+nLayer1Nodes+nLayer2Nodes;

                    i9 =i1+1;
                    i10=i2+(2*Nx+1-i);
                    i11=i3-1;
                    i12=i10-1;

                    i13= i9+nLayer1Nodes+nLayer2Nodes;
                    i14=i10+nLayer1Nodes+nLayer2Nodes;
                    i15=i11+nLayer1Nodes+nLayer2Nodes;
                    i16=i12+nLayer1Nodes+nLayer2Nodes;

                    i17=i1+nLayer1Nodes-(i-1+(j-1)*(Nx+Nx+1));
                    i18=i17+1;
                    i19=i18+Nx+1;
                    i20=i19-1;


                    Conn[(e-1)*nNodesPerElmt+1-1]=i1;
                    Conn[(e-1)*nNodesPerElmt+2-1]=i2;
                    Conn[(e-1)*nNodesPerElmt+3-1]=i3;
                    Conn[(e-1)*nNodesPerElmt+4-1]=i4;
                    Conn[(e-1)*nNodesPerElmt+5-1]=i5;
                    Conn[(e-1)*nNodesPerElmt+6-1]=i6;
                    Conn[(e-1)*nNodesPerElmt+7-1]=i7;
                    Conn[(e-1)*nNodesPerElmt+8-1]=i8;
                    Conn[(e-1)*nNodesPerElmt+9-1]=i9;
                    Conn[(e-1)*nNodesPerElmt+10-1]=i10;

                    Conn[(e-1)*nNodesPerElmt+11-1]=i11;
                    Conn[(e-1)*nNodesPerElmt+12-1]=i12;
                    Conn[(e-1)*nNodesPerElmt+13-1]=i13;
                    Conn[(e-1)*nNodesPerElmt+14-1]=i14;
                    Conn[(e-1)*nNodesPerElmt+15-1]=i15;
                    Conn[(e-1)*nNodesPerElmt+16-1]=i16;
                    Conn[(e-1)*nNodesPerElmt+17-1]=i17;
                    Conn[(e-1)*nNodesPerElmt+18-1]=i18;
                    Conn[(e-1)*nNodesPerElmt+19-1]=i19;
                    Conn[(e-1)*nNodesPerElmt+20-1]=i20;

                }
            }
        }

        VTKCellType=25;
    }
    else if(MeshType=="hex27")
    {
        // for 3D-27 nodes mesh
        dx=(Xmax-Xmin)/(2.0*Nx);
        dy=(Ymax-Ymin)/(2.0*Ny);
        dz=(Zmax-Zmin)/(2.0*Nz);


        nElmts=Nx*Ny*Nz;
        int nLayerNodes=(2*Nx+1)*(2*Ny+1);// for norm layer



        nNodes=(2*Nz+1)*nLayerNodes;
        nNodesPerElmt=27;
        nNodesPerBCElmt=9;


        NodeCoords.resize(4*nNodes,0.0);
        Conn.resize(nElmts*nNodesPerElmt,0);



        for(k=1;k<=Nz;++k)
        {
            // For first layer
            for(j=1;j<=2*Ny+1;++j)
            {
                // for bottom line of each element
                for(i=1;i<=2*Nx+1;++i)
                {
                    kk=(j-1)*(2*Nx+1)+i+(k-1)*2*nLayerNodes;
                    NodeCoords[(kk-1)*4+0]=1.0;
                    NodeCoords[(kk-1)*4+1]=Xmin+(i-1)*dx;
                    NodeCoords[(kk-1)*4+2]=Ymin+(j-1)*dy;
                    NodeCoords[(kk-1)*4+3]=Zmin+(k-1)*2*dz;
                }
            }
            // Then for second layer
            for(j=1;j<=2*Ny+1;++j)
            {
                for(i=1;i<=2*Nx+1;++i)
                {
                    kk=(j-1)*(2*Nx+1)+i+(k-1)*2*nLayerNodes+nLayerNodes;
                    NodeCoords[(kk-1)*4+0]=1.0;
                    NodeCoords[(kk-1)*4+1]=Xmin+(i-1)*dx;
                    NodeCoords[(kk-1)*4+2]=Ymin+(j-1)*dy;
                    NodeCoords[(kk-1)*4+3]=Zmin+(k-1)*2*dz+dz;
                }
            }
        }

        // for the last top layer
        k=Nz+1;
        for(j=1;j<=2*Ny+1;++j)
        {
            // for bottom line of each element
            for(i=1;i<=2*Nx+1;++i)
            {
                kk=(j-1)*(2*Nx+1)+i+(k-1)*2*nLayerNodes;
                NodeCoords[(kk-1)*4+0]=1.0;
                NodeCoords[(kk-1)*4+1]=Xmin+(i-1)*dx;
                NodeCoords[(kk-1)*4+2]=Ymin+(j-1)*dy;
                NodeCoords[(kk-1)*4+3]=Zmin+(k-1)*2*dz;
            }
        }



        // Create Connectivity matrix
        for(k=1;k<=Nz;++k)
        {
            for(j=1;j<=Ny;++j)
            {
                for(i=1;i<=Nx;++i)
                {
                    e=(j-1)*Nx+i+(k-1)*Nx*Ny;
                    i1=(j-1)*2*(2*Nx+1)+2*i-1+(k-1)*2*nLayerNodes;
                    i2=i1+2;
                    i3=i2+(2*Nx+1)*2;
                    i4=i3-2;

                    i5=i1+2*nLayerNodes;
                    i6=i2+2*nLayerNodes;
                    i7=i3+2*nLayerNodes;
                    i8=i4+2*nLayerNodes;

                    i9 =i1+1;
                    i10=i2+(2*Nx+1);
                    i11=i3-1;
                    i12=i1+(2*Nx+1);

                    i13=i5+1;
                    i14=i6+(2*Nx+1);
                    i15=i7-1;
                    i16=i5+(2*Nx+1);

                    i17=i1+nLayerNodes;
                    i18=i2+nLayerNodes;
                    i19=i3+nLayerNodes;
                    i20=i4+nLayerNodes;

                    i21=i17+(2*Nx+1);
                    i22=i21+2;

                    //i23=i20+1;
                    //i24=i17+1;

                    i23=i17+1;
                    i24=i20+1;

                    i25=i12+1;
                    i26=i16+1;

                    i27=i21+1;


                    Conn[(e-1)*nNodesPerElmt+1-1]=i1;
                    Conn[(e-1)*nNodesPerElmt+2-1]=i2;
                    Conn[(e-1)*nNodesPerElmt+3-1]=i3;
                    Conn[(e-1)*nNodesPerElmt+4-1]=i4;
                    Conn[(e-1)*nNodesPerElmt+5-1]=i5;
                    Conn[(e-1)*nNodesPerElmt+6-1]=i6;
                    Conn[(e-1)*nNodesPerElmt+7-1]=i7;
                    Conn[(e-1)*nNodesPerElmt+8-1]=i8;
                    Conn[(e-1)*nNodesPerElmt+9-1]=i9;
                    Conn[(e-1)*nNodesPerElmt+10-1]=i10;

                    Conn[(e-1)*nNodesPerElmt+11-1]=i11;
                    Conn[(e-1)*nNodesPerElmt+12-1]=i12;
                    Conn[(e-1)*nNodesPerElmt+13-1]=i13;
                    Conn[(e-1)*nNodesPerElmt+14-1]=i14;
                    Conn[(e-1)*nNodesPerElmt+15-1]=i15;
                    Conn[(e-1)*nNodesPerElmt+16-1]=i16;
                    Conn[(e-1)*nNodesPerElmt+17-1]=i17;
                    Conn[(e-1)*nNodesPerElmt+18-1]=i18;
                    Conn[(e-1)*nNodesPerElmt+19-1]=i19;
                    Conn[(e-1)*nNodesPerElmt+20-1]=i20;

                    Conn[(e-1)*nNodesPerElmt+21-1]=i21;
                    Conn[(e-1)*nNodesPerElmt+22-1]=i22;
                    Conn[(e-1)*nNodesPerElmt+23-1]=i23;
                    Conn[(e-1)*nNodesPerElmt+24-1]=i24;
                    Conn[(e-1)*nNodesPerElmt+25-1]=i25;
                    Conn[(e-1)*nNodesPerElmt+26-1]=i26;
                    Conn[(e-1)*nNodesPerElmt+27-1]=i27;
                }
            }
        }
        VTKCellType=29;
    }

    // TODO: split mesh for 3D mesh

    MeshCreated=true;
    BCMeshCreated=false;
}

