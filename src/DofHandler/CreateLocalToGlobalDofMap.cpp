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
// create dofmap in AsFem

#include "DofHandler/DofHandler.h"

bool DofHandler::CreateLocalToGlobalDofMap(Mesh &mesh, int ndofspernode)
{
    pair<string,vector<int>> temp;
    HasDofMap=false;
    GlobalDofMap.clear();

    nDims=mesh.GetDims();
    nNodes=mesh.GetNodesNum();
    nElmts=mesh.GetElmtsNum();
    nNodesPerElmts=mesh.GetNodesNumPerElmt();

    nDofsPerNode=ndofspernode;
    nDofs=nNodes*nDofsPerNode;
    nDofsPerElmt=nNodesPerElmts*nDofsPerNode;
    nDofsPerBCElmt=nNodesPerBCElmt*nDofsPerNode;

    GlobalDofMap.resize(nElmts*nDofsPerElmt,0);

    int e,i,j,k,iInd;
    for(e=1;e<=nElmts;e++)
    {
        for(i=1;i<=nNodesPerElmts;i++)
        {
            k=mesh.GetIthConnJthIndex(e,i);
            for(j=1;j<=nDofsPerNode;j++)
            {
                iInd=(k-1)*nDofsPerNode+j;
                GlobalDofMap[(e-1)*nDofsPerElmt+(i-1)*nDofsPerNode+j-1]=iInd;
            }
        }
    }
    GlobalDofMap.resize(GlobalDofMap.size());
    HasDofMap=true;

    // Now to create boundary's dof map
    nNodesPerBCElmt=mesh.GetNodesNumPerBCElmt();
    nDofsPerBCElmt=nNodesPerBCElmt*nDofsPerNode;
    if(mesh.GetDims()==1)
    {
        LeftBCDofs.resize(mesh.GetSideBCElmtNum("left")*nNodesPerBCElmt*nDofsPerNode,0);
        for(e=1;e<=mesh.GetSideBCElmtNum("left");e++)
        {
            for(i=1;i<=nNodesPerBCElmt;i++)
            {
                j=mesh.GetSideBCIthConnJthIndex("left",e,i);
                for(k=1;k<=nDofsPerNode;k++)
                {
                    iInd=(j-1)*nDofsPerNode+k;
                    LeftBCDofs[(e-1)*nNodesPerBCElmt*nDofsPerNode+(i-1)*nDofsPerNode+k-1]=iInd;
                }
            }
        }

        //*** For right side
        RightBCDofs.resize(mesh.GetSideBCElmtNum("right")*nNodesPerBCElmt*nDofsPerNode,0);
        for(e=1;e<=mesh.GetSideBCElmtNum("right");e++)
        {
            for(i=1;i<=nNodesPerBCElmt;i++)
            {
                j=mesh.GetSideBCIthConnJthIndex("right",e,i);
                for(k=1;k<=nDofsPerNode;k++)
                {
                    iInd=(j-1)*nDofsPerNode+k;
                    RightBCDofs[(e-1)*nNodesPerBCElmt*nDofsPerNode+(i-1)*nDofsPerNode+k-1]=iInd;
                }
            }
        }

        temp=make_pair("left",LeftBCDofs);
        GlobalBCDofMap.push_back(temp);

        temp=make_pair("right",RightBCDofs);
        GlobalBCDofMap.push_back(temp);

        HasBCDofMap=true;
    }
    else if(mesh.GetDims()==2)
    {
        // for left side
        LeftBCDofs.resize(mesh.GetSideBCElmtNum("left")*nNodesPerBCElmt*nDofsPerNode,0);
        for(e=1;e<=mesh.GetSideBCElmtNum("left");e++)
        {
            for(i=1;i<=nNodesPerBCElmt;i++)
            {
                j=mesh.GetSideBCIthConnJthIndex("left",e,i);
                for(k=1;k<=nDofsPerNode;k++)
                {
                    iInd=(j-1)*nDofsPerNode+k;
                    LeftBCDofs[(e-1)*nNodesPerBCElmt*nDofsPerNode+(i-1)*nDofsPerNode+k-1]=iInd;
                }
            }
        }
        //*** For right side
        RightBCDofs.resize(mesh.GetSideBCElmtNum("right")*nNodesPerBCElmt*nDofsPerNode,0);
        for(e=1;e<=mesh.GetSideBCElmtNum("right");e++)
        {
            for(i=1;i<=nNodesPerBCElmt;i++)
            {
                j=mesh.GetSideBCIthConnJthIndex("right",e,i);
                for(k=1;k<=nDofsPerNode;k++)
                {
                    iInd=(j-1)*nDofsPerNode+k;
                    RightBCDofs[(e-1)*nNodesPerBCElmt*nDofsPerNode+(i-1)*nDofsPerNode+k-1]=iInd;
                }
            }
        }
        // for bottom side
        BottomBCDofs.resize(mesh.GetSideBCElmtNum("bottom")*nNodesPerBCElmt*nDofsPerNode,0);
        for(e=1;e<=mesh.GetSideBCElmtNum("bottom");e++)
        {
            for(i=1;i<=nNodesPerBCElmt;i++)
            {
                j=mesh.GetSideBCIthConnJthIndex("bottom",e,i);
                for(k=1;k<=nDofsPerNode;k++)
                {
                    iInd=(j-1)*nDofsPerNode+k;
                    BottomBCDofs[(e-1)*nNodesPerBCElmt*nDofsPerNode+(i-1)*nDofsPerNode+k-1]=iInd;
                }
            }
        }
        // for top side
        TopBCDofs.resize(mesh.GetSideBCElmtNum("top")*nNodesPerBCElmt*nDofsPerNode,0);
        for(e=1;e<=mesh.GetSideBCElmtNum("top");e++)
        {
            for(i=1;i<=nNodesPerBCElmt;i++)
            {
                j=mesh.GetSideBCIthConnJthIndex("top",e,i);
                for(k=1;k<=nDofsPerNode;k++)
                {
                    iInd=(j-1)*nDofsPerNode+k;
                    TopBCDofs[(e-1)*nNodesPerBCElmt*nDofsPerNode+(i-1)*nDofsPerNode+k-1]=iInd;
                }
            }
        }



        temp=make_pair("left",LeftBCDofs);
        GlobalBCDofMap.push_back(temp);

        temp=make_pair("right",RightBCDofs);
        GlobalBCDofMap.push_back(temp);

        temp=make_pair("bottom",BottomBCDofs);
        GlobalBCDofMap.push_back(temp);

        temp=make_pair("top",TopBCDofs);
        GlobalBCDofMap.push_back(temp);

        HasBCDofMap=true;
    }
    else if(mesh.GetDims()==3)
    {
        // for left side
        LeftBCDofs.resize(mesh.GetSideBCElmtNum("left")*nNodesPerBCElmt*nDofsPerNode,0);
        for(e=1;e<=mesh.GetSideBCElmtNum("left");e++)
        {
            for(i=1;i<=nNodesPerBCElmt;i++)
            {
                j=mesh.GetSideBCIthConnJthIndex("left",e,i);
                for(k=1;k<=nDofsPerNode;k++)
                {
                    iInd=(j-1)*nDofsPerNode+k;
                    LeftBCDofs[(e-1)*nNodesPerBCElmt*nDofsPerNode+(i-1)*nDofsPerNode+k-1]=iInd;
                }
            }
        }
        //*** For right side
        RightBCDofs.resize(mesh.GetSideBCElmtNum("right")*nNodesPerBCElmt*nDofsPerNode,0);
        for(e=1;e<=mesh.GetSideBCElmtNum("right");e++)
        {
            for(i=1;i<=nNodesPerBCElmt;i++)
            {
                j=mesh.GetSideBCIthConnJthIndex("right",e,i);
                for(k=1;k<=nDofsPerNode;k++)
                {
                    iInd=(j-1)*nDofsPerNode+k;
                    RightBCDofs[(e-1)*nNodesPerBCElmt*nDofsPerNode+(i-1)*nDofsPerNode+k-1]=iInd;
                }
            }
        }
        // for bottom side
        BottomBCDofs.resize(mesh.GetSideBCElmtNum("bottom")*nNodesPerBCElmt*nDofsPerNode,0);
        for(e=1;e<=mesh.GetSideBCElmtNum("bottom");e++)
        {
            for(i=1;i<=nNodesPerBCElmt;i++)
            {
                j=mesh.GetSideBCIthConnJthIndex("bottom",e,i);
                for(k=1;k<=nDofsPerNode;k++)
                {
                    iInd=(j-1)*nDofsPerNode+k;
                    BottomBCDofs[(e-1)*nNodesPerBCElmt*nDofsPerNode+(i-1)*nDofsPerNode+k-1]=iInd;
                }
            }
        }
        // for top side
        TopBCDofs.resize(mesh.GetSideBCElmtNum("top")*nNodesPerBCElmt*nDofsPerNode,0);
        for(e=1;e<=mesh.GetSideBCElmtNum("top");e++)
        {
            for(i=1;i<=nNodesPerBCElmt;i++)
            {
                j=mesh.GetSideBCIthConnJthIndex("top",e,i);
                for(k=1;k<=nDofsPerNode;k++)
                {
                    iInd=(j-1)*nDofsPerNode+k;
                    TopBCDofs[(e-1)*nNodesPerBCElmt*nDofsPerNode+(i-1)*nDofsPerNode+k-1]=iInd;
                }
            }
        }
        // for back side
        BackBCDofs.resize(mesh.GetSideBCElmtNum("back")*nNodesPerBCElmt*nDofsPerNode,0);
        for(e=1;e<=mesh.GetSideBCElmtNum("back");e++)
        {
            for(i=1;i<=nNodesPerBCElmt;i++)
            {
                j=mesh.GetSideBCIthConnJthIndex("back",e,i);
                for(k=1;k<=nDofsPerNode;k++)
                {
                    iInd=(j-1)*nDofsPerNode+k;
                    BackBCDofs[(e-1)*nNodesPerBCElmt*nDofsPerNode+(i-1)*nDofsPerNode+k-1]=iInd;
                }
            }
        }
        // for front side
        FrontBCDofs.resize(mesh.GetSideBCElmtNum("front")*nNodesPerBCElmt*nDofsPerNode,0);
        for(e=1;e<=mesh.GetSideBCElmtNum("front");e++)
        {
            for(i=1;i<=nNodesPerBCElmt;i++)
            {
                j=mesh.GetSideBCIthConnJthIndex("front",e,i);
                for(k=1;k<=nDofsPerNode;k++)
                {
                    iInd=(j-1)*nDofsPerNode+k;
                    FrontBCDofs[(e-1)*nNodesPerBCElmt*nDofsPerNode+(i-1)*nDofsPerNode+k-1]=iInd;
                }
            }
        }

        temp=make_pair("left",LeftBCDofs);
        GlobalBCDofMap.push_back(temp);

        temp=make_pair("right",RightBCDofs);
        GlobalBCDofMap.push_back(temp);

        temp=make_pair("bottom",BottomBCDofs);
        GlobalBCDofMap.push_back(temp);

        temp=make_pair("top",TopBCDofs);
        GlobalBCDofMap.push_back(temp);

        temp=make_pair("back",BackBCDofs);
        GlobalBCDofMap.push_back(temp);

        temp=make_pair("front",FrontBCDofs);
        GlobalBCDofMap.push_back(temp);

        HasBCDofMap=true;
    }
    return HasBCDofMap;
}