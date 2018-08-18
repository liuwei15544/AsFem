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
// get mesh information of AsFem

#include "Mesh/Mesh.h"

int Mesh::GetIthConnJthIndex(int e, int j) const
{
    if(!MeshCreated)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: mesh is not generated !!!            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        can't get local connect info!!!      ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                 ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscFinalize();
        abort();
    }

    if(e<1||e>nElmts)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%8d is out of nElmts=%8d !         ***\n",e,nElmts);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                 ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscFinalize();
        abort();
    }

    if(j<1||j>nNodesPerElmt)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: j=%3d is out of nNodesPerElmt=%3d !!!***\n",j,nNodesPerElmt);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        j should be 1~%3d !!!                ***\n",nNodesPerElmt);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                 ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscFinalize();
        abort();
    }

    return Conn[(e-1)*nNodesPerElmt+j-1];
}
//***********************************
double Mesh::GetIthNodeJthCoord(int i, int j) const
{
    if(!MeshCreated)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: mesh is not generated !!!            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        can't get local connect info!!!      ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                 ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscFinalize();
        abort();
    }

    if(i<1||i>nNodes)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: i=%6d is out of nNodes=%6d!          ***\n",i,nNodes);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                                 ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***************************************************\n");
        PetscFinalize();
        abort();
    }
    if(j<0||j>3)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: j=%2d is out 0~3\n",i,nNodes);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        j should be 0->w,1->x,2->y,3->z!\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    return NodeCoords[4*(i-1)+j];
}
//*************************************
void Mesh::GetLocalCoords(int e, double (&coords)[27][4]) const
{
    if(!MeshCreated)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: mesh is not generated !!!       ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"***        can't get local connect info!!! ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }

    if(e<1||e>nElmts)
    {
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%8d is out of nElmts=%8d\n",e,nElmts);
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
        PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
        PetscFinalize();
        abort();
    }
    int i,j;
    for(i=1;i<=nNodesPerElmt;i++)
    {
        j=GetIthConnJthIndex(e,i);
        coords[i-1][0]=GetIthNodeJthCoord(j,0);
        coords[i-1][1]=GetIthNodeJthCoord(j,1);
        coords[i-1][2]=GetIthNodeJthCoord(j,2);
        coords[i-1][3]=GetIthNodeJthCoord(j,3);
    }
}

//******************************************
//*** For boundary mesh information
//******************************************
int Mesh::GetSideBCElmtNum(string sidename) const
{
    if(sidename=="left")
    {
        return int(LeftBCConn.size()/nNodesPerBCElmt);
    }
    else if(sidename=="right")
    {
        return int(RightBCConn.size()/nNodesPerBCElmt);
    }
    else if(sidename=="bottom")
    {
        return int(BottomBCConn.size()/nNodesPerBCElmt);
    }
    else if(sidename=="top")
    {
        return int(TopBCConn.size()/nNodesPerBCElmt);
    }
    else
    {
        bool IsSideNameInList=false;
        for(unsigned int i=0;i<BCMeshSet.size();i++)
        {
            if(BCMeshSet[i].first==sidename)
            {
                IsSideNameInList=true;
                return int(BCMeshSet[i].second.size()/nNodesPerBCElmt);
            }
        }

        if(!IsSideNameInList)
        {
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: sidename=%12s is not in the list!!!",sidename.c_str());
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
            PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
            PetscFinalize();
            abort();
        }
    }

    return 0;
}

//**********************************************
int Mesh::GetSideBCIthConnJthIndex(string sidename, int i, int j) const
{
    bool IsSideNameInList=false;
    if(nDim==1)
    {
        if(sidename=="left")
        {
            return LeftBCConn[(i-1)*nNodesPerBCElmt+j-1];
        }
        else if(sidename=="right")
        {
            return RightBCConn[(i-1)*nNodesPerBCElmt+j-1];
        }
        else
        {
            IsSideNameInList=false;
            for(unsigned int e=0;e<BCMeshSet.size();e++)
            {
                if(BCMeshSet[e].first==sidename)
                {
                    IsSideNameInList=true;
                    return BCMeshSet[e].second[(i-1)*nNodesPerBCElmt+j-1];
                }
            }

            if(!IsSideNameInList)
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: sidename=%12s is not in the list!!!",sidename);
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
        }
    }
    else if(nDim==2)
    {
        if(sidename=="left")
        {
            return LeftBCConn[(i-1)*nNodesPerBCElmt+j-1];
        }
        else if(sidename=="right")
        {
            return RightBCConn[(i-1)*nNodesPerBCElmt+j-1];
        }
        else if(sidename=="bottom")
        {
            return BottomBCConn[(i-1)*nNodesPerBCElmt+j-1];
        }
        else if(sidename=="top")
        {
            return TopBCConn[(i-1)*nNodesPerBCElmt+j-1];
        }
        else
        {
            IsSideNameInList=false;
            for(unsigned int e=0;e<BCMeshSet.size();e++)
            {
                if(BCMeshSet[e].first==sidename)
                {
                    IsSideNameInList=true;
                    return BCMeshSet[e].second[(i-1)*nNodesPerBCElmt+j-1];
                }
            }

            if(!IsSideNameInList)
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: sidename=%12s is not in the list!!!",sidename.c_str());
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
        }
    }
    else if(nDim==3)
    {
        if(sidename=="left")
        {
            return LeftBCConn[(i-1)*nNodesPerBCElmt+j-1];
        }
        else if(sidename=="right")
        {
            return RightBCConn[(i-1)*nNodesPerBCElmt+j-1];
        }
        else if(sidename=="bottom")
        {
            return BottomBCConn[(i-1)*nNodesPerBCElmt+j-1];
        }
        else if(sidename=="top")
        {
            return TopBCConn[(i-1)*nNodesPerBCElmt+j-1];
        }
        else if(sidename=="back")
        {
            return BackBCConn[(i-1)*nNodesPerBCElmt+j-1];
        }
        else if(sidename=="front")
        {
            return BackBCConn[(i-1)*nNodesPerBCElmt+j-1];
        }
        else
        {
            IsSideNameInList=false;
            for(unsigned int e=0;e<BCMeshSet.size();e++)
            {
                if(BCMeshSet[e].first==sidename)
                {
                    IsSideNameInList=true;
                    return BCMeshSet[e].second[(i-1)*nNodesPerBCElmt+j-1];
                }
            }

            if(!IsSideNameInList)
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: sidename=%12s is not in the list!!!",sidename.c_str());
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
        }
    }
    return -1;
}

//******************************
double Mesh::GetSideBCIthNodeJthCoord(string sidename, int e, int i, int j) const
{
    bool IsSideNameInList;
    int ii,jj;
    if(nDim==1)
    {
        if(sidename=="left")
        {
            ii=GetSideBCIthConnJthIndex("left",e,i);
            return GetIthNodeJthCoord(ii,j);
        }
        else if(sidename=="right")
        {
            ii=GetSideBCIthConnJthIndex("right",e,i);
            return GetIthNodeJthCoord(ii,j);
        }
        else
        {
            IsSideNameInList=false;
            for(unsigned int ii=0;ii<BCMeshSet.size();ii++)
            {
                if(BCMeshSet[ii].first==sidename)
                {
                    IsSideNameInList=true;
                    jj=BCMeshSet[ii].second[(e-1)*nNodesPerBCElmt+i-1];
                    return GetIthNodeJthCoord(jj,j);
                }
            }

            if(!IsSideNameInList)
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: sidename=%12s is not in the list!!!",sidename.c_str());
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
        }
    }
    else if(nDim==2)
    {
        if(sidename=="left")
        {
            ii=GetSideBCIthConnJthIndex("left",e,i);
            return GetIthNodeJthCoord(ii,j);
        }
        else if(sidename=="right")
        {
            ii=GetSideBCIthConnJthIndex("right",e,i);
            return GetIthNodeJthCoord(ii,j);
        }
        else if(sidename=="bottom")
        {
            ii=GetSideBCIthConnJthIndex("bottom",e,i);
            return GetIthNodeJthCoord(ii,j);
        }
        else if(sidename=="top")
        {
            ii=GetSideBCIthConnJthIndex("top",e,i);
            return GetIthNodeJthCoord(ii,j);
        }
        else
        {
            IsSideNameInList=false;
            for(unsigned int ii=0;ii<BCMeshSet.size();ii++)
            {
                if(BCMeshSet[ii].first==sidename)
                {
                    IsSideNameInList=true;
                    jj=BCMeshSet[ii].second[(e-1)*nNodesPerBCElmt+i-1];
                    return GetIthNodeJthCoord(jj,j);
                }
            }

            if(!IsSideNameInList)
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: sidename=%12s is not in the list!!!",sidename.c_str());
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
        }
    }
    else if(nDim==3)
    {
        if(sidename=="left")
        {
            ii=GetSideBCIthConnJthIndex("left",e,i);
            return GetIthNodeJthCoord(ii,j);
        }
        else if(sidename=="right")
        {
            ii=GetSideBCIthConnJthIndex("right",e,i);
            return GetIthNodeJthCoord(ii,j);
        }
        else if(sidename=="bottom")
        {
            ii=GetSideBCIthConnJthIndex("bottom",e,i);
            return GetIthNodeJthCoord(ii,j);
        }
        else if(sidename=="top")
        {
            ii=GetSideBCIthConnJthIndex("top",e,i);
            return GetIthNodeJthCoord(ii,j);
        }
        else if(sidename=="back")
        {
            ii=GetSideBCIthConnJthIndex("back",e,i);
            return GetIthNodeJthCoord(ii,j);
        }
        else if(sidename=="front")
        {
            ii=GetSideBCIthConnJthIndex("front",e,i);
            return GetIthNodeJthCoord(ii,j);
        }
        else
        {
            IsSideNameInList=false;
            for(unsigned int ii=0;ii<BCMeshSet.size();ii++)
            {
                if(BCMeshSet[ii].first==sidename)
                {
                    IsSideNameInList=true;
                    jj=BCMeshSet[ii].second[(e-1)*nNodesPerBCElmt+i-1];
                    return GetIthNodeJthCoord(jj,j);
                }
            }

            if(!IsSideNameInList)
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: sidename=%12s is not in the list!!!",sidename.c_str());
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
        }
    }

    return 0.0;
}

//*****************************************
void Mesh::GetLocalBCCoords(string sidename, int e, double (&coords)[27][4]) const
{
    bool IsSideNameInList;
    int i,j;
    if(nDim==1)
    {
        if(sidename=="left")
        {
            if(e<1||e>GetSideBCElmtNum("left"))
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%8d is out of nElmts=%8d\n",e,GetSideBCElmtNum("left"));
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
            else
            {
                for(i=1;i<=nNodesPerBCElmt;i++)
                {
                    j=GetSideBCIthConnJthIndex("left",e,i);
                    coords[i-1][0]=GetIthNodeJthCoord(j,0);
                    coords[i-1][1]=GetIthNodeJthCoord(j,1);
                    coords[i-1][2]=GetIthNodeJthCoord(j,2);
                    coords[i-1][3]=GetIthNodeJthCoord(j,3);
                }
            }
        }
        else if(sidename=="right")
        {
            if(e<1||e>GetSideBCElmtNum("right"))
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%8d is out of nElmts=%8d\n",e,GetSideBCElmtNum("right"));
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
            else
            {
                for(i=1;i<=nNodesPerBCElmt;i++)
                {
                    j=GetSideBCIthConnJthIndex("right",e,i);
                    coords[i-1][0]=GetIthNodeJthCoord(j,0);
                    coords[i-1][1]=GetIthNodeJthCoord(j,1);
                    coords[i-1][2]=GetIthNodeJthCoord(j,2);
                    coords[i-1][3]=GetIthNodeJthCoord(j,3);
                }
            }
        }
        else
        {
            IsSideNameInList=false;
            for(unsigned int ii=0;ii<BCMeshSet.size();ii++)
            {
                if(BCMeshSet[ii].first==sidename)
                {
                    IsSideNameInList=true;
                    for(i=1;i<=nNodesPerBCElmt;i++)
                    {
                        j=GetSideBCIthConnJthIndex(sidename,e,i);
                        coords[i-1][0]=GetIthNodeJthCoord(j,0);
                        coords[i-1][1]=GetIthNodeJthCoord(j,1);
                        coords[i-1][2]=GetIthNodeJthCoord(j,2);
                        coords[i-1][3]=GetIthNodeJthCoord(j,3);
                    }
                }
            }

            if(!IsSideNameInList)
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: sidename=%12s is not in the list!!!",sidename.c_str());
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
        }
    }
    else if(nDim==2)
    {
        if(sidename=="left")
        {
            if(e<1||e>GetSideBCElmtNum("left"))
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%8d is out of nElmts=%8d\n",e,GetSideBCElmtNum("left"));
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
            else
            {
                for(i=1;i<=nNodesPerBCElmt;i++)
                {
                    j=GetSideBCIthConnJthIndex("left",e,i);
                    coords[i-1][0]=GetIthNodeJthCoord(j,0);
                    coords[i-1][1]=GetIthNodeJthCoord(j,1);
                    coords[i-1][2]=GetIthNodeJthCoord(j,2);
                    coords[i-1][3]=GetIthNodeJthCoord(j,3);
                }
            }
        }
        else if(sidename=="right")
        {
            if(e<1||e>GetSideBCElmtNum("right"))
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%8d is out of nElmts=%8d\n",e,GetSideBCElmtNum("right"));
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
            else
            {
                for(i=1;i<=nNodesPerBCElmt;i++)
                {
                    j=GetSideBCIthConnJthIndex("right",e,i);
                    coords[i-1][0]=GetIthNodeJthCoord(j,0);
                    coords[i-1][1]=GetIthNodeJthCoord(j,1);
                    coords[i-1][2]=GetIthNodeJthCoord(j,2);
                    coords[i-1][3]=GetIthNodeJthCoord(j,3);
                }
            }
        }
        else if(sidename=="bottom")
        {
            if(e<1||e>GetSideBCElmtNum("bottom"))
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%8d is out of nElmts=%8d\n",e,GetSideBCElmtNum("bottom"));
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
            else
            {
                for(i=1;i<=nNodesPerBCElmt;i++)
                {
                    j=GetSideBCIthConnJthIndex("bottom",e,i);
                    coords[i-1][0]=GetIthNodeJthCoord(j,0);
                    coords[i-1][1]=GetIthNodeJthCoord(j,1);
                    coords[i-1][2]=GetIthNodeJthCoord(j,2);
                    coords[i-1][3]=GetIthNodeJthCoord(j,3);
                }
            }
        }
        else if(sidename=="top")
        {
            if(e<1||e>GetSideBCElmtNum("top"))
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%8d is out of nElmts=%8d\n",e,GetSideBCElmtNum("top"));
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
            else
            {
                for(i=1;i<=nNodesPerBCElmt;i++)
                {
                    j=GetSideBCIthConnJthIndex("top",e,i);
                    coords[i-1][0]=GetIthNodeJthCoord(j,0);
                    coords[i-1][1]=GetIthNodeJthCoord(j,1);
                    coords[i-1][2]=GetIthNodeJthCoord(j,2);
                    coords[i-1][3]=GetIthNodeJthCoord(j,3);
                }
            }
        }
        else
        {
            IsSideNameInList=false;
            for(unsigned int ii=0;ii<BCMeshSet.size();ii++)
            {
                if(BCMeshSet[ii].first==sidename)
                {
                    IsSideNameInList=true;
                    for(i=1;i<=nNodesPerBCElmt;i++)
                    {
                        j=GetSideBCIthConnJthIndex(sidename,e,i);
                        coords[i-1][0]=GetIthNodeJthCoord(j,0);
                        coords[i-1][1]=GetIthNodeJthCoord(j,1);
                        coords[i-1][2]=GetIthNodeJthCoord(j,2);
                        coords[i-1][3]=GetIthNodeJthCoord(j,3);
                    }
                }
            }

            if(!IsSideNameInList)
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: sidename=%12s is not in the list!!!",sidename.c_str());
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
        }
    }
    else if(nDim==3)
    {
        if(sidename=="left")
        {
            if(e<1||e>GetSideBCElmtNum("left"))
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%8d is out of nElmts=%8d\n",e,GetSideBCElmtNum("left"));
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
            else
            {
                for(i=1;i<=nNodesPerBCElmt;i++)
                {
                    j=GetSideBCIthConnJthIndex("left",e,i);
                    coords[i-1][0]=GetIthNodeJthCoord(j,0);
                    coords[i-1][1]=GetIthNodeJthCoord(j,1);
                    coords[i-1][2]=GetIthNodeJthCoord(j,2);
                    coords[i-1][3]=GetIthNodeJthCoord(j,3);
                }
            }
        }
        else if(sidename=="right")
        {
            if(e<1||e>GetSideBCElmtNum("right"))
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%8d is out of nElmts=%8d\n",e,GetSideBCElmtNum("right"));
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
            else
            {
                for(i=1;i<=nNodesPerBCElmt;i++)
                {
                    j=GetSideBCIthConnJthIndex("right",e,i);
                    coords[i-1][0]=GetIthNodeJthCoord(j,0);
                    coords[i-1][1]=GetIthNodeJthCoord(j,1);
                    coords[i-1][2]=GetIthNodeJthCoord(j,2);
                    coords[i-1][3]=GetIthNodeJthCoord(j,3);
                }
            }
        }
        else if(sidename=="bottom")
        {
            if(e<1||e>GetSideBCElmtNum("bottom"))
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%8d is out of nElmts=%8d\n",e,GetSideBCElmtNum("bottom"));
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
            else
            {
                for(i=1;i<=nNodesPerBCElmt;i++)
                {
                    j=GetSideBCIthConnJthIndex("bottom",e,i);
                    coords[i-1][0]=GetIthNodeJthCoord(j,0);
                    coords[i-1][1]=GetIthNodeJthCoord(j,1);
                    coords[i-1][2]=GetIthNodeJthCoord(j,2);
                    coords[i-1][3]=GetIthNodeJthCoord(j,3);
                }
            }
        }
        else if(sidename=="top")
        {
            if(e<1||e>GetSideBCElmtNum("top"))
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%8d is out of nElmts=%8d\n",e,GetSideBCElmtNum("top"));
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
            else
            {
                for(i=1;i<=nNodesPerBCElmt;i++)
                {
                    j=GetSideBCIthConnJthIndex("top",e,i);
                    coords[i-1][0]=GetIthNodeJthCoord(j,0);
                    coords[i-1][1]=GetIthNodeJthCoord(j,1);
                    coords[i-1][2]=GetIthNodeJthCoord(j,2);
                    coords[i-1][3]=GetIthNodeJthCoord(j,3);
                }
            }
        }
        else if(sidename=="back")
        {
            if(e<1||e>GetSideBCElmtNum("back"))
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%8d is out of nElmts=%8d\n",e,GetSideBCElmtNum("back"));
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
            else
            {
                for(i=1;i<=nNodesPerBCElmt;i++)
                {
                    j=GetSideBCIthConnJthIndex("back",e,i);
                    coords[i-1][0]=GetIthNodeJthCoord(j,0);
                    coords[i-1][1]=GetIthNodeJthCoord(j,1);
                    coords[i-1][2]=GetIthNodeJthCoord(j,2);
                    coords[i-1][3]=GetIthNodeJthCoord(j,3);
                }
            }
        }
        else if(sidename=="front")
        {
            if(e<1||e>GetSideBCElmtNum("front"))
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: e=%8d is out of nElmts=%8d\n",e,GetSideBCElmtNum("front"));
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
            else
            {
                for(i=1;i<=nNodesPerBCElmt;i++)
                {
                    j=GetSideBCIthConnJthIndex("front",e,i);
                    coords[i-1][0]=GetIthNodeJthCoord(j,0);
                    coords[i-1][1]=GetIthNodeJthCoord(j,1);
                    coords[i-1][2]=GetIthNodeJthCoord(j,2);
                    coords[i-1][3]=GetIthNodeJthCoord(j,3);
                }
            }
        }
        else
        {
            IsSideNameInList=false;
            for(unsigned int ii=0;ii<BCMeshSet.size();ii++)
            {
                if(BCMeshSet[ii].first==sidename)
                {
                    IsSideNameInList=true;
                    for(i=1;i<=nNodesPerBCElmt;i++)
                    {
                        j=GetSideBCIthConnJthIndex(sidename,e,i);
                        coords[i-1][0]=GetIthNodeJthCoord(j,0);
                        coords[i-1][1]=GetIthNodeJthCoord(j,1);
                        coords[i-1][2]=GetIthNodeJthCoord(j,2);
                        coords[i-1][3]=GetIthNodeJthCoord(j,3);
                    }
                }
            }

            if(!IsSideNameInList)
            {
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** Error: sidename=%12s is not in the list!!!",sidename.c_str());
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"*** AsFem exit!                            ***\n");
                PetscSynchronizedPrintf(PETSC_COMM_WORLD,"**********************************************\n");
                PetscFinalize();
                abort();
            }
        }
    }
}