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
// Created by walkandthinker on 18.08.18.
// read mesh block from input file

#include "InputSystem/InputSystem.h"

bool InputSystem::ReadMeshBlock(Mesh &mesh)
{
    int linenum=0;
    string str,line,line0="mesh";
    int blockstartlinenum;
    vector<double> numbers;
    int dim;
    double xmin,xmax,ymin,ymax,zmin,zmax;
    string meshtype;
    bool ReadMeshSuccess=false;
    int nx,ny,nz;

    // Read the first comment line
    getline(in,line);linenum+=1;
    if(IsBracketMatch(in,line0,blockstartlinenum))
    {

        GotoLine(in,blockstartlinenum);linenum=blockstartlinenum-1;
        getline(in,line);linenum+=1;// read [mesh]

        str=RemoveSpace(line);
        while(str.compare(0,2,"[]")!=0)
        {
            // read type
            getline(in,line);linenum+=1;
            str=RemoveSpace(line);
            while(str.size()<1)
            {
                getline(in,line);linenum+=1;
                str=RemoveSpace(line);
            }
            if(str.compare(0,10,"type=asfem")==0)
            {
                // read built-in mesh
                getline(in,line);linenum+=1;// read dim
                str=RemoveSpace(line);
                while(str.size()<1)
                {
                    getline(in,line);linenum+=1;
                    str=RemoveSpace(line);
                }
                if(str.compare(0,4,"dim=")!=0)
                {
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find 'dim=' line-%-3d      ***\n",linenum);
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given dim=1,2,3 !!!  ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    return false;
                }
                numbers=SplitNum(line);
                if(numbers.size()<1)
                {
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given dim=1,2,3 !!!  ***\n");
                    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                    return false;
                }
                else
                {
                    dim=int(numbers[0]);
                    if(dim<1||dim>3)
                    {
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: dim=%2d is invalid in line-%-3d!!!***\n",dim,linenum);
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given dim=1,2,3 !!!  ***\n");
                        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                        return false;
                    }
                    else if(dim==1)
                    {
                        mesh.SetDim(1);

                        getline(in,line);linenum+=1;// read xmin
                        str=RemoveSpace(line);
                        while(str.size()<1)
                        {
                            getline(in,line);linenum+=1;
                            str=RemoveSpace(line);
                        }
                        if(str.compare(0,5,"xmin=")==0)
                        {
                            numbers=SplitNum(line);
                            if(numbers.size()<1)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given dim=1,2,3 !!!  ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                            xmin=numbers[0];
                        }
                        else
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find xmin=  in line-%3d***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given xmin=value !!! ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }

                        //***for xmax
                        getline(in,line);linenum+=1;// read xmax
                        str=RemoveSpace(line);
                        while(str.size()<1)
                        {
                            getline(in,line);linenum+=1;
                            str=RemoveSpace(line);
                        }
                        if(str.compare(0,5,"xmax=")==0)
                        {
                            numbers=SplitNum(line);
                            if(numbers.size()<1)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given dim=1,2,3 !!!  ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                            xmax=numbers[0];
                        }
                        else
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find xmax=  in line-%3d***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given xmax=value !!! ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }

                        if(xmin>xmax)
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: invalid mesh information        ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        xmin=%8.3f>xmax=%8.3f !!! ***\n",xmin,xmax);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        please check your input file    ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }

                        mesh.SetXmin(xmin);
                        mesh.SetXmax(xmax);

                        //***for nx
                        getline(in,line);linenum+=1;// read nx
                        str=RemoveSpace(line);
                        while(str.size()<1)
                        {
                            getline(in,line);linenum+=1;
                            str=RemoveSpace(line);
                        }
                        if(str.compare(0,3,"nx=")==0)
                        {
                            numbers=SplitNum(line);
                            if(numbers.size()<1)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given nx=int value   ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                            nx=numbers[0];
                            if(nx<1)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: nx=%5d(line-%3d) is invalid !!***\n",nx,linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given nx=value(>0)   ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                        }
                        else
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: 'nx=' can't found in line=%-3d   ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given nx=value(>0)   ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                        mesh.SetNx(nx);

                        //**** Read mesh type
                        getline(in,line);linenum+=1;//
                        str=RemoveSpace(line);
                        while(str.size()<1)
                        {
                            getline(in,line);linenum+=1;
                            str=RemoveSpace(line);
                        }
                        if(str.compare(0,9,"meshtype=")==0)
                        {
                            if(str.find("edge2")!=string::npos)
                            {
                                meshtype="edge2";
                            }
                            else if(str.find("edge3")!=string::npos)
                            {
                                meshtype="edge3";
                            }
                            else if(str.find("edge4")!=string::npos)
                            {
                                meshtype="edge4";
                            }
                            else
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported mesh type(line-%-3d) ***\n",linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        meshtype=edge2,3,4 is expected!!***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                        }

                        mesh.SetMeshType(meshtype);

                        mesh.CreateMesh();
                        ReadMeshSuccess=true;
                        return true;
                    }
                    else if(dim==2)
                    {
                        mesh.SetDim(2);

                        getline(in,line);linenum+=1;// read xmin
                        str=RemoveSpace(line);
                        while(str.size()<1)
                        {
                            getline(in,line);linenum+=1;
                            str=RemoveSpace(line);
                        }
                        if(str.compare(0,5,"xmin=")==0)
                        {
                            numbers=SplitNum(line);
                            if(numbers.size()<1)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given dim=1,2,3 !!!  ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                            xmin=numbers[0];
                        }
                        else
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find 'xmin='  in line-%-3d ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given xmin=value !!! ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }

                        //***for xmax
                        getline(in,line);linenum+=1;// read xmax
                        str=RemoveSpace(line);
                        while(str.size()<1)
                        {
                            getline(in,line);linenum+=1;
                            str=RemoveSpace(line);
                        }
                        if(str.compare(0,5,"xmax=")==0)
                        {
                            numbers=SplitNum(line);
                            if(numbers.size()<1)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given dim=1,2,3 !!!  ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                            xmax=numbers[0];
                        }
                        else
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find 'xmax='  in line-%-3d ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given xmax=value !!! ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }

                        if(xmin>xmax)
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: invalid mesh information        ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        xmin=%8.3f>xmax=%8.3f !!! ***\n",xmin,xmax);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        please check your input file    ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                        mesh.SetXmin(xmin);
                        mesh.SetXmax(xmax);

                        //****************************
                        getline(in,line);linenum+=1;// read ymin
                        str=RemoveSpace(line);
                        while(str.size()<1)
                        {
                            getline(in,line);linenum+=1;
                            str=RemoveSpace(line);
                        }
                        if(str.compare(0,5,"ymin=")==0)
                        {
                            numbers=SplitNum(line);
                            if(numbers.size()<1)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given dim=1,2,3 !!!  ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                            ymin=numbers[0];
                        }
                        else
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find 'ymin=' in line-%3d  ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given ymin=value !!! ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }

                        //***for ymax
                        getline(in,line);linenum+=1;// read ymax
                        str=RemoveSpace(line);
                        while(str.size()<1)
                        {
                            getline(in,line);linenum+=1;
                            str=RemoveSpace(line);
                        }
                        if(str.compare(0,5,"ymax=")==0)
                        {
                            numbers=SplitNum(line);
                            if(numbers.size()<1)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given dim=1,2,3 !!!  ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                            ymax=numbers[0];
                        }
                        else
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find 'ymax=' in line-%3d  ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given ymax=value !!! ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }

                        if(ymin>ymax)
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: invalid mesh information        ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        ymin=%8.3f>ymax=%8.3f !!! ***\n",ymin,ymax);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        please check your input file    ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                        mesh.SetYmin(ymin);
                        mesh.SetYmax(ymax);

                        //***for nx
                        getline(in,line);linenum+=1;// read nx
                        str=RemoveSpace(line);
                        while(str.size()<1)
                        {
                            getline(in,line);linenum+=1;
                            str=RemoveSpace(line);
                        }
                        if(str.compare(0,3,"nx=")==0)
                        {
                            numbers=SplitNum(line);
                            if(numbers.size()<1)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given nx=int value   ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                            nx=numbers[0];
                            if(nx<1)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: nx=%5d(line-%3d) is invalid !!***\n",nx,linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given nx=value(>0)   ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                        }
                        else
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: 'nx=' can't found in line=%-3d   ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given nx=value(>0)   ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                        mesh.SetNx(nx);

                        //***for ny
                        getline(in,line);linenum+=1;// read ny
                        str=RemoveSpace(line);
                        while(str.size()<1)
                        {
                            getline(in,line);linenum+=1;
                            str=RemoveSpace(line);
                        }
                        if(str.compare(0,3,"ny=")==0)
                        {
                            numbers=SplitNum(line);
                            if(numbers.size()<1)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given nx=int value   ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                            ny=numbers[0];
                            if(ny<1)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: ny=%5d(line-%3d) is invalid !!***\n",nx,linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given ny=value(>0)   ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                        }
                        else
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: 'ny=' can't found in line=%-3d   ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given ny=value(>0)   ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                        mesh.SetNy(ny);

                        //**** Read mesh type
                        getline(in,line);linenum+=1;//
                        str=RemoveSpace(line);
                        while(str.size()<1)
                        {
                            getline(in,line);linenum+=1;
                            str=RemoveSpace(line);
                        }
                        if(str.compare(0,9,"meshtype=")==0)
                        {
                            if(str.find("quad4")!=string::npos)
                            {
                                meshtype="quad4";
                            }
                            else if(str.find("quad8")!=string::npos)
                            {
                                meshtype="quad8";
                            }
                            else if(str.find("quad9")!=string::npos)
                            {
                                meshtype="quad9";
                            }
                            else
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported mesh type(line-%-3d) ***\n",linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***       'meshtype=quad4,8,9'is expected!!***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                        }
                        else
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error:can't find 'meshtype=' in line%-3d***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***       'meshtype=quad4,8,9'is expected!!***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                        mesh.SetMeshType(meshtype);

                        mesh.CreateMesh();
                        ReadMeshSuccess=true;
                        return ReadMeshSuccess;
                    }
                    else if(dim==3)
                    {
                        mesh.SetDim(3);

                        getline(in,line);linenum+=1;// read xmin
                        str=RemoveSpace(line);
                        while(str.size()<1)
                        {
                            getline(in,line);linenum+=1;
                            str=RemoveSpace(line);
                        }
                        if(str.compare(0,5,"xmin=")==0)
                        {
                            numbers=SplitNum(line);
                            if(numbers.size()<1)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given dim=1,2,3 !!!  ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                            xmin=numbers[0];
                        }
                        else
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find 'xmin='  in line-%-3d ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given xmin=value !!! ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                        //***for xmax
                        getline(in,line);linenum+=1;// read xmax
                        str=RemoveSpace(line);
                        while(str.size()<1)
                        {
                            getline(in,line);linenum+=1;
                            str=RemoveSpace(line);
                        }
                        if(str.compare(0,5,"xmax=")==0)
                        {
                            numbers=SplitNum(line);
                            if(numbers.size()<1)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given dim=1,2,3 !!!  ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                            xmax=numbers[0];
                        }
                        else
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find 'xmax='  in line-%-3d ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given xmax=value !!! ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }

                        if(xmin>xmax)
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: invalid mesh information        ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        xmin=%8.3f>xmax=%8.3f !!! ***\n",xmin,xmax);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        please check your input file    ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                        mesh.SetXmin(xmin);
                        mesh.SetXmax(xmax);

                        //****************************
                        getline(in,line);linenum+=1;// read ymin
                        str=RemoveSpace(line);
                        while(str.size()<1)
                        {
                            getline(in,line);linenum+=1;
                            str=RemoveSpace(line);
                        }
                        if(str.compare(0,5,"ymin=")==0)
                        {
                            numbers=SplitNum(line);
                            if(numbers.size()<1)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given ymin=value !!! ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                            ymin=numbers[0];
                        }
                        else
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find 'ymin=' in line-%3d  ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given ymin=value !!! ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }

                        //***for ymax
                        getline(in,line);linenum+=1;// read ymax
                        str=RemoveSpace(line);
                        while(str.size()<1)
                        {
                            getline(in,line);linenum+=1;
                            str=RemoveSpace(line);
                        }
                        if(str.compare(0,5,"ymax=")==0)
                        {
                            numbers=SplitNum(line);
                            if(numbers.size()<1)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given ymax=value !!! ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                            ymax=numbers[0];
                        }
                        else
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find 'ymax=' in line-%3d  ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given ymax=value !!! ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }

                        if(ymin>ymax)
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: invalid mesh information        ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        ymin=%8.3f>ymax=%8.3f !!! ***\n",ymin,ymax);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        please check your input file    ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                        mesh.SetYmin(ymin);
                        mesh.SetYmax(ymax);

                        //****************************
                        getline(in,line);linenum+=1;// read zmin
                        str=RemoveSpace(line);
                        while(str.size()<1)
                        {
                            getline(in,line);linenum+=1;
                            str=RemoveSpace(line);
                        }
                        if(str.compare(0,5,"zmin=")==0)
                        {
                            numbers=SplitNum(line);
                            if(numbers.size()<1)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given zmin=value !!! ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                            zmin=numbers[0];
                        }
                        else
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find 'zmin=' in line-%3d  ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given zmin=value !!! ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }

                        //***for zmax
                        getline(in,line);linenum+=1;// read zmax
                        str=RemoveSpace(line);
                        while(str.size()<1)
                        {
                            getline(in,line);linenum+=1;
                            str=RemoveSpace(line);
                        }
                        if(str.compare(0,5,"zmax=")==0)
                        {
                            numbers=SplitNum(line);
                            if(numbers.size()<1)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given zmax=value !!! ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                            zmax=numbers[0];
                        }
                        else
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find 'zmax=' in line-%3d  ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given zmax=value !!! ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }

                        if(zmin>zmax)
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: invalid mesh information        ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        zmin=%8.3f>zmax=%8.3f !!! ***\n",zmin,zmax);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        please check your input file    ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                        mesh.SetZmin(zmin);
                        mesh.SetZmax(zmax);

                        //***for nx
                        getline(in,line);linenum+=1;// read nx
                        str=RemoveSpace(line);
                        while(str.size()<1)
                        {
                            getline(in,line);linenum+=1;
                            str=RemoveSpace(line);
                        }
                        if(str.compare(0,3,"nx=")==0)
                        {
                            numbers=SplitNum(line);
                            if(numbers.size()<1)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given nx=int value   ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                            nx=numbers[0];
                            if(nx<1)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: nx=%5d(line-%3d) is invalid !!***\n",nx,linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given nx=value(>0)   ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                        }
                        else
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: 'nx=' can't found in line=%-3d   ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given nx=value(>0)   ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                        mesh.SetNx(nx);

                        //***for ny
                        getline(in,line);linenum+=1;// read ny
                        str=RemoveSpace(line);
                        while(str.size()<1)
                        {
                            getline(in,line);linenum+=1;
                            str=RemoveSpace(line);
                        }
                        if(str.compare(0,3,"ny=")==0)
                        {
                            numbers=SplitNum(line);
                            if(numbers.size()<1)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given ny=int value   ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                            ny=numbers[0];
                            if(ny<1)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: ny=%5d(line-%3d) is invalid !!***\n",nx,linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given ny=value(>0)   ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                        }
                        else
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: 'ny=' can't found in line=%-3d   ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given ny=value(>0)   ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                        mesh.SetNy(ny);

                        //***for nz
                        getline(in,line);linenum+=1;// read nz
                        str=RemoveSpace(line);
                        while(str.size()<1)
                        {
                            getline(in,line);linenum+=1;
                            str=RemoveSpace(line);
                        }
                        if(str.compare(0,3,"nz=")==0)
                        {
                            numbers=SplitNum(line);
                            if(numbers.size()<1)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: can't find number in line-%-3d  ***\n",linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given nz=int value   ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                            nz=numbers[0];
                            if(nz<1)
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: nz=%5d(line-%3d) is invalid !!***\n",nz,linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given nz=value(>0)   ***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                        }
                        else
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: 'nz=' can't found in line=%-3d   ***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        you should given nz=value(>0)   ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                        mesh.SetNz(nz);

                        //**** Read mesh type
                        getline(in,line);linenum+=1;//
                        str=RemoveSpace(line);
                        while(str.size()<1)
                        {
                            getline(in,line);linenum+=1;
                            str=RemoveSpace(line);
                        }
                        if(str.compare(0,9,"meshtype=")==0)
                        {
                            if(str.find("hex8")!=string::npos)
                            {
                                meshtype="hex8";
                            }
                            else if(str.find("hex20")!=string::npos)
                            {
                                meshtype="hex20";
                            }
                            else if(str.find("hex27")!=string::npos)
                            {
                                meshtype="hex27";
                            }
                            else
                            {
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: unsupported mesh type(line-%-3d) ***\n",linenum);
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        meshtype=hex8,20,27 is expected!***\n");
                                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                                return false;
                            }
                        }
                        else
                        {
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error:can't find 'meshtype=' in line%-3d***\n",linenum);
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***       meshtype=hex8,20,27 is expected! ***\n");
                            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                            return false;
                        }
                        mesh.SetMeshType(meshtype);

                        mesh.CreateMesh();

                        ReadMeshSuccess=true;
                        return ReadMeshSuccess;
                    }
                }

            }
        }
    }
    else
    {
        ReadMeshSuccess=false;
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***----------------------------------------***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***--- can't find matched [mesh]/[] block--***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***--- please check your input file!!!   --***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!!!                          ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    return false;
}