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
// Created by walkandthinker on 25.08.18.
// Define the output system in AsFem

#include "OutputSystem/OutputSystem.h"

void OutputSystem::WriteUToVTUFile(Mesh &mesh,EquationSystem &equationSystem, Vec &U)
{
    string filename;

    // for ghost node
    VecScatterCreateToAll(U,&scatter,&Useq);
    VecScatterBegin(scatter,U,Useq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(scatter,U,Useq,INSERT_VALUES,SCATTER_FORWARD);

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    if(rank==0)
    {
        if(InputFileName.size()<1)
        {
            filename="result.vtu";
        }
        else if(InputFileName.size()<5)
        {
            filename=InputFileName+".vtu";
        }
        else
        {
            filename=InputFileName.substr(0,5)+".vtu";
        }

        int i,j,iInd,e;

        PetscFOpen(PETSC_COMM_SELF,filename.c_str(),"w",&fd);

        PetscFPrintf(PETSC_COMM_SELF,fd,"<?xml version=\"1.0\"?>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"<UnstructuredGrid>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"<Piece NumberOfPoints=\"%8d\" NumberOfCells=\"%8d\">\n",mesh.GetNodesNum(),mesh.GetElmtsNum());
        PetscFPrintf(PETSC_COMM_SELF,fd,"<Points>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"<DataArray type=\"Float64\" Name=\"nodes\"  NumberOfComponents=\"3\"  format=\"ascii\">\n");

        for(i=1;i<=mesh.GetNodesNum();++i)
        {
            PetscFPrintf(PETSC_COMM_SELF,fd,"%14.6e %14.6e %14.6e\n",
                         mesh.GetIthNodeJthCoord(i,1),
                         mesh.GetIthNodeJthCoord(i,2),
                         mesh.GetIthNodeJthCoord(i,3));
        }
        PetscFPrintf(PETSC_COMM_SELF,fd,"</DataArray>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"</Points>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"<Cells>\n");

        //***********************************************************
        PetscFPrintf(PETSC_COMM_SELF,fd,"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        for(e=1;e<=mesh.GetElmtsNum();++e)
        {
            for(j=1;j<=mesh.GetNodesNumPerElmt();++j)
            {
                PetscFPrintf(PETSC_COMM_SELF,fd,"%8d ",mesh.GetIthConnJthIndex(e,j)-1);
            }
            PetscFPrintf(PETSC_COMM_SELF,fd,"\n");
        }
        PetscFPrintf(PETSC_COMM_SELF,fd,"</DataArray>\n");


        PetscFPrintf(PETSC_COMM_SELF,fd,"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        PetscInt offset=0;
        for(e=1;e<=mesh.GetElmtsNum();++e)
        {
            offset+=mesh.GetNodesNumPerElmt();
            PetscFPrintf(PETSC_COMM_SELF,fd,"%8d\n",offset);
        }
        PetscFPrintf(PETSC_COMM_SELF,fd,"</DataArray>\n");

        PetscFPrintf(PETSC_COMM_SELF,fd,"<DataArray type=\"Int32\" Name=\"types\"  NumberOfComponents=\"1\"  format=\"ascii\">\n");
        PetscInt VTKCellType=mesh.GetVTKCellType();
        for(e=1;e<=mesh.GetElmtsNum();++e)
        {
            PetscFPrintf(PETSC_COMM_SELF,fd,"%3d\n",VTKCellType);
        }
        PetscFPrintf(PETSC_COMM_SELF,fd,"</DataArray>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"</Cells>\n");

        PetscFPrintf(PETSC_COMM_SELF,fd,"<PointData Scalar=\"sol\">\n");
        PetscScalar value;
        for(j=1;j<=equationSystem.GetDofsNumPerNode();j++)
        {
            PetscFPrintf(PETSC_COMM_SELF,fd,"<DataArray type=\"Float64\"  Name=\" %s \"  NumberOfComponents=\"1\" format=\"ascii\">\n",equationSystem.GetDofNameByIndex(j).c_str());
            for(i=0;i<mesh.GetNodesNum();++i)
            {
                iInd=i*equationSystem.GetDofsNumPerNode()+j-1;
                VecGetValues(Useq,1,&iInd,&value);
                PetscFPrintf(PETSC_COMM_SELF,fd,"%14.6e\n",value);
            }
            PetscFPrintf(PETSC_COMM_SELF,fd,"</DataArray>\n\n");
        }
        PetscFPrintf(PETSC_COMM_SELF,fd,"</PointData>\n");

        PetscFPrintf(PETSC_COMM_SELF,fd,"</Piece>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"</UnstructuredGrid>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"</VTKFile>\n");


        PetscFClose(PETSC_COMM_SELF,fd);
    }

    VecScatterDestroy(&scatter);
    VecDestroy(&Useq);
}

//********************************************************
void OutputSystem::WriteUToVTUFile(int step,Mesh &mesh,EquationSystem &equationSystem, Vec &U)
{
    string filename;

    // for ghost node
    VecScatterCreateToAll(U,&scatter,&Useq);
    VecScatterBegin(scatter,U,Useq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(scatter,U,Useq,INSERT_VALUES,SCATTER_FORWARD);

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    if(rank==0)
    {
        if(InputFileName.size()<1)
        {
            filename="result";
        }
        else if(InputFileName.size()<5)
        {
            filename=InputFileName;
        }
        else
        {
            filename=InputFileName.substr(0,5);
        }

        ostringstream ss;
        ss<<setfill('0')<<setw(8)<<step;
        filename=filename+ss.str()+".vtu";

        int i,j,iInd,e;

        PetscFOpen(PETSC_COMM_SELF,filename.c_str(),"w",&fd);

        PetscFPrintf(PETSC_COMM_SELF,fd,"<?xml version=\"1.0\"?>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"<UnstructuredGrid>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"<Piece NumberOfPoints=\"%8d\" NumberOfCells=\"%8d\">\n",mesh.GetNodesNum(),mesh.GetElmtsNum());
        PetscFPrintf(PETSC_COMM_SELF,fd,"<Points>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"<DataArray type=\"Float64\" Name=\"nodes\"  NumberOfComponents=\"3\"  format=\"ascii\">\n");

        for(i=1;i<=mesh.GetNodesNum();++i)
        {
            PetscFPrintf(PETSC_COMM_SELF,fd,"%14.6e %14.6e %14.6e\n",
                         mesh.GetIthNodeJthCoord(i,1),
                         mesh.GetIthNodeJthCoord(i,2),
                         mesh.GetIthNodeJthCoord(i,3));
        }
        PetscFPrintf(PETSC_COMM_SELF,fd,"</DataArray>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"</Points>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"<Cells>\n");

        //***********************************************************
        PetscFPrintf(PETSC_COMM_SELF,fd,"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        for(e=1;e<=mesh.GetElmtsNum();++e)
        {
            for(j=1;j<=mesh.GetNodesNumPerElmt();++j)
            {
                PetscFPrintf(PETSC_COMM_SELF,fd,"%8d ",mesh.GetIthConnJthIndex(e,j)-1);
            }
            PetscFPrintf(PETSC_COMM_SELF,fd,"\n");
        }
        PetscFPrintf(PETSC_COMM_SELF,fd,"</DataArray>\n");


        PetscFPrintf(PETSC_COMM_SELF,fd,"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        PetscInt offset=0;
        for(e=1;e<=mesh.GetElmtsNum();++e)
        {
            offset+=mesh.GetNodesNumPerElmt();
            PetscFPrintf(PETSC_COMM_SELF,fd,"%8d\n",offset);
        }
        PetscFPrintf(PETSC_COMM_SELF,fd,"</DataArray>\n");

        PetscFPrintf(PETSC_COMM_SELF,fd,"<DataArray type=\"Int32\" Name=\"types\"  NumberOfComponents=\"1\"  format=\"ascii\">\n");
        PetscInt VTKCellType=mesh.GetVTKCellType();
        for(e=1;e<=mesh.GetElmtsNum();++e)
        {
            PetscFPrintf(PETSC_COMM_SELF,fd,"%3d\n",VTKCellType);
        }
        PetscFPrintf(PETSC_COMM_SELF,fd,"</DataArray>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"</Cells>\n");

        PetscFPrintf(PETSC_COMM_SELF,fd,"<PointData Scalar=\"sol\">\n");
        PetscScalar value;
        for(j=1;j<=equationSystem.GetDofsNumPerNode();j++)
        {
            PetscFPrintf(PETSC_COMM_SELF,fd,"<DataArray type=\"Float64\"  Name=\" %s \"  NumberOfComponents=\"1\" format=\"ascii\">\n",equationSystem.GetDofNameByIndex(j).c_str());
            for(i=0;i<mesh.GetNodesNum();++i)
            {
                iInd=i*equationSystem.GetDofsNumPerNode()+j-1;
                VecGetValues(Useq,1,&iInd,&value);
                PetscFPrintf(PETSC_COMM_SELF,fd,"%14.6e\n",value);
            }
            PetscFPrintf(PETSC_COMM_SELF,fd,"</DataArray>\n\n");
        }
        PetscFPrintf(PETSC_COMM_SELF,fd,"</PointData>\n");

        PetscFPrintf(PETSC_COMM_SELF,fd,"</Piece>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"</UnstructuredGrid>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"</VTKFile>\n");


        PetscFClose(PETSC_COMM_SELF,fd);
    }

    VecScatterDestroy(&scatter);
    VecDestroy(&Useq);
}

//****************************************
void OutputSystem::WriteUAndProjToVTUFile(Mesh &mesh,
                                          EquationSystem &equationSystem,
                                          Vec &U, Vec &Proj)
{
    string filename;

    // for ghost node
    VecScatterCreateToAll(U,&scatter,&Useq);
    VecScatterBegin(scatter,U,Useq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(scatter,U,Useq,INSERT_VALUES,SCATTER_FORWARD);

    VecScatterCreateToAll(Proj,&scatterproj,&PROJseq);
    VecScatterBegin(scatterproj,Proj,PROJseq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(scatterproj,Proj,PROJseq,INSERT_VALUES,SCATTER_FORWARD);

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    if(rank==0)
    {
        if(InputFileName.size()<1)
        {
            filename="result.vtu";
        }
        else if(InputFileName.size()<5)
        {
            filename=InputFileName+".vtu";
        }
        else
        {
            filename=InputFileName.substr(0,5)+".vtu";
        }

        int i,j,iInd,e;

        PetscFOpen(PETSC_COMM_SELF,filename.c_str(),"w",&fd);

        PetscFPrintf(PETSC_COMM_SELF,fd,"<?xml version=\"1.0\"?>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"<UnstructuredGrid>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"<Piece NumberOfPoints=\"%8d\" NumberOfCells=\"%8d\">\n",mesh.GetNodesNum(),mesh.GetElmtsNum());
        PetscFPrintf(PETSC_COMM_SELF,fd,"<Points>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"<DataArray type=\"Float64\" Name=\"nodes\"  NumberOfComponents=\"3\"  format=\"ascii\">\n");

        for(i=1;i<=mesh.GetNodesNum();++i)
        {
            PetscFPrintf(PETSC_COMM_SELF,fd,"%14.6e %14.6e %14.6e\n",
                         mesh.GetIthNodeJthCoord(i,1),
                         mesh.GetIthNodeJthCoord(i,2),
                         mesh.GetIthNodeJthCoord(i,3));
        }
        PetscFPrintf(PETSC_COMM_SELF,fd,"</DataArray>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"</Points>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"<Cells>\n");

        //***********************************************************
        PetscFPrintf(PETSC_COMM_SELF,fd,"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        for(e=1;e<=mesh.GetElmtsNum();++e)
        {
            for(j=1;j<=mesh.GetNodesNumPerElmt();++j)
            {
                PetscFPrintf(PETSC_COMM_SELF,fd,"%8d ",mesh.GetIthConnJthIndex(e,j)-1);
            }
            PetscFPrintf(PETSC_COMM_SELF,fd,"\n");
        }
        PetscFPrintf(PETSC_COMM_SELF,fd,"</DataArray>\n");


        PetscFPrintf(PETSC_COMM_SELF,fd,"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        PetscInt offset=0;
        for(e=1;e<=mesh.GetElmtsNum();++e)
        {
            offset+=mesh.GetNodesNumPerElmt();
            PetscFPrintf(PETSC_COMM_SELF,fd,"%8d\n",offset);
        }
        PetscFPrintf(PETSC_COMM_SELF,fd,"</DataArray>\n");

        PetscFPrintf(PETSC_COMM_SELF,fd,"<DataArray type=\"Int32\" Name=\"types\"  NumberOfComponents=\"1\"  format=\"ascii\">\n");
        PetscInt VTKCellType=mesh.GetVTKCellType();
        for(e=1;e<=mesh.GetElmtsNum();++e)
        {
            PetscFPrintf(PETSC_COMM_SELF,fd,"%3d\n",VTKCellType);
        }
        PetscFPrintf(PETSC_COMM_SELF,fd,"</DataArray>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"</Cells>\n");

        PetscFPrintf(PETSC_COMM_SELF,fd,"<PointData Scalar=\"sol\">\n");
        PetscScalar value,weight;
        //*******************************
        //*** for solution
        //*******************************
        for(j=1;j<=equationSystem.GetDofsNumPerNode();j++)
        {
            PetscFPrintf(PETSC_COMM_SELF,fd,"<DataArray type=\"Float64\"  Name=\" %s \"  NumberOfComponents=\"1\" format=\"ascii\">\n",equationSystem.GetDofNameByIndex(j).c_str());
            for(i=0;i<mesh.GetNodesNum();++i)
            {
                iInd=i*equationSystem.GetDofsNumPerNode()+j-1;
                VecGetValues(Useq,1,&iInd,&value);
                PetscFPrintf(PETSC_COMM_SELF,fd,"%14.6e\n",value);
            }
            PetscFPrintf(PETSC_COMM_SELF,fd,"</DataArray>\n\n");
        }

        //*******************************
        //*** for projection
        //*******************************
        for(j=1;j<=12;j++)
        {
            PetscFPrintf(PETSC_COMM_SELF,fd,"<DataArray type=\"Float64\"  Name=\"Proj-%2d \"  NumberOfComponents=\"1\" format=\"ascii\">\n",j);
            for(i=0;i<mesh.GetNodesNum();++i)
            {
                iInd=i*13+j;
                VecGetValues(PROJseq,1,&iInd,&value);
                iInd=i*13+0;
                VecGetValues(PROJseq,1,&iInd,&weight);
                if(fabs(value/weight)>1.0e-12)
                {
                    PetscFPrintf(PETSC_COMM_SELF,fd,"%14.6e\n",value/weight);
                }
                else
                {
                    PetscFPrintf(PETSC_COMM_SELF,fd,"%14.6e\n",1.0e-12);
                }
            }
            PetscFPrintf(PETSC_COMM_SELF,fd,"</DataArray>\n\n");
        }
        PetscFPrintf(PETSC_COMM_SELF,fd,"</PointData>\n");

        PetscFPrintf(PETSC_COMM_SELF,fd,"</Piece>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"</UnstructuredGrid>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"</VTKFile>\n");


        PetscFClose(PETSC_COMM_SELF,fd);
    }

    VecScatterDestroy(&scatter);
    VecScatterDestroy(&scatterproj);
    VecDestroy(&Useq);
    VecDestroy(&PROJseq);
}

//****************************************
void OutputSystem::WriteUAndProjToVTUFile(int step,Mesh &mesh,
                                          EquationSystem &equationSystem,
                                          Vec &U, Vec &Proj)
{
    string filename;

    // for ghost node
    VecScatterCreateToAll(U,&scatter,&Useq);
    VecScatterBegin(scatter,U,Useq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(scatter,U,Useq,INSERT_VALUES,SCATTER_FORWARD);

    VecScatterCreateToAll(Proj,&scatterproj,&PROJseq);
    VecScatterBegin(scatterproj,Proj,PROJseq,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(scatterproj,Proj,PROJseq,INSERT_VALUES,SCATTER_FORWARD);

    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    if(rank==0)
    {
        if(InputFileName.size()<1)
        {
            filename="result";
        }
        else if(InputFileName.size()<5)
        {
            filename=InputFileName;
        }
        else
        {
            filename=InputFileName.substr(0,5);
        }
        ostringstream ss;
        ss<<setfill('0')<<setw(8)<<step;
        filename=filename+ss.str()+".vtu";


        int i,j,iInd,e;

        PetscFOpen(PETSC_COMM_SELF,filename.c_str(),"w",&fd);

        PetscFPrintf(PETSC_COMM_SELF,fd,"<?xml version=\"1.0\"?>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"<UnstructuredGrid>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"<Piece NumberOfPoints=\"%8d\" NumberOfCells=\"%8d\">\n",mesh.GetNodesNum(),mesh.GetElmtsNum());
        PetscFPrintf(PETSC_COMM_SELF,fd,"<Points>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"<DataArray type=\"Float64\" Name=\"nodes\"  NumberOfComponents=\"3\"  format=\"ascii\">\n");

        for(i=1;i<=mesh.GetNodesNum();++i)
        {
            PetscFPrintf(PETSC_COMM_SELF,fd,"%14.6e %14.6e %14.6e\n",
                         mesh.GetIthNodeJthCoord(i,1),
                         mesh.GetIthNodeJthCoord(i,2),
                         mesh.GetIthNodeJthCoord(i,3));
        }
        PetscFPrintf(PETSC_COMM_SELF,fd,"</DataArray>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"</Points>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"<Cells>\n");

        //***********************************************************
        PetscFPrintf(PETSC_COMM_SELF,fd,"<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        for(e=1;e<=mesh.GetElmtsNum();++e)
        {
            for(j=1;j<=mesh.GetNodesNumPerElmt();++j)
            {
                PetscFPrintf(PETSC_COMM_SELF,fd,"%8d ",mesh.GetIthConnJthIndex(e,j)-1);
            }
            PetscFPrintf(PETSC_COMM_SELF,fd,"\n");
        }
        PetscFPrintf(PETSC_COMM_SELF,fd,"</DataArray>\n");


        PetscFPrintf(PETSC_COMM_SELF,fd,"<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n");
        PetscInt offset=0;
        for(e=1;e<=mesh.GetElmtsNum();++e)
        {
            offset+=mesh.GetNodesNumPerElmt();
            PetscFPrintf(PETSC_COMM_SELF,fd,"%8d\n",offset);
        }
        PetscFPrintf(PETSC_COMM_SELF,fd,"</DataArray>\n");

        PetscFPrintf(PETSC_COMM_SELF,fd,"<DataArray type=\"Int32\" Name=\"types\"  NumberOfComponents=\"1\"  format=\"ascii\">\n");
        PetscInt VTKCellType=mesh.GetVTKCellType();
        for(e=1;e<=mesh.GetElmtsNum();++e)
        {
            PetscFPrintf(PETSC_COMM_SELF,fd,"%3d\n",VTKCellType);
        }
        PetscFPrintf(PETSC_COMM_SELF,fd,"</DataArray>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"</Cells>\n");

        PetscFPrintf(PETSC_COMM_SELF,fd,"<PointData Scalar=\"sol\">\n");
        PetscScalar value,weight;
        //*******************************
        //*** for solution
        //*******************************
        for(j=1;j<=equationSystem.GetDofsNumPerNode();j++)
        {
            PetscFPrintf(PETSC_COMM_SELF,fd,"<DataArray type=\"Float64\"  Name=\" %s \"  NumberOfComponents=\"1\" format=\"ascii\">\n",equationSystem.GetDofNameByIndex(j).c_str());
            for(i=0;i<mesh.GetNodesNum();++i)
            {
                iInd=i*equationSystem.GetDofsNumPerNode()+j-1;
                VecGetValues(Useq,1,&iInd,&value);
                PetscFPrintf(PETSC_COMM_SELF,fd,"%14.6e\n",value);
            }
            PetscFPrintf(PETSC_COMM_SELF,fd,"</DataArray>\n\n");
        }

        //*******************************
        //*** for projection
        //*******************************
        for(j=1;j<=12;j++)
        {
            PetscFPrintf(PETSC_COMM_SELF,fd,"<DataArray type=\"Float64\"  Name=\"Proj-%2d \"  NumberOfComponents=\"1\" format=\"ascii\">\n",j);
            for(i=0;i<mesh.GetNodesNum();++i)
            {
                iInd=i*13+j;
                VecGetValues(PROJseq,1,&iInd,&value);
                iInd=i*13+0;
                VecGetValues(PROJseq,1,&iInd,&weight);
                if(step==0)
                {
                    PetscFPrintf(PETSC_COMM_SELF,fd,"%14.6e\n",value);
                }
                else
                {
                    if(fabs(value/weight)>1.0e-12)
                    {
                        PetscFPrintf(PETSC_COMM_SELF,fd,"%14.6e\n",value/weight);
                    }
                    else
                    {
                        PetscFPrintf(PETSC_COMM_SELF,fd,"%14.6e\n",1.0e-12);
                    }
                }
            }
            PetscFPrintf(PETSC_COMM_SELF,fd,"</DataArray>\n\n");
        }
        PetscFPrintf(PETSC_COMM_SELF,fd,"</PointData>\n");

        PetscFPrintf(PETSC_COMM_SELF,fd,"</Piece>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"</UnstructuredGrid>\n");
        PetscFPrintf(PETSC_COMM_SELF,fd,"</VTKFile>\n");


        PetscFClose(PETSC_COMM_SELF,fd);
    }

    VecScatterDestroy(&scatter);
    VecScatterDestroy(&scatterproj);
    VecDestroy(&Useq);
    VecDestroy(&PROJseq);
}


