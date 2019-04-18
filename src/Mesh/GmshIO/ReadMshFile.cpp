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
// Created by walkandthinker on 18.04.19.
// read mesh for .msh file

#include "Mesh/GmshIO.h"
#include "MsgPrint/MsgPrintForInput.h"

bool GmshIO::ReadMshFile(vector<double> &NodeCoords,
                         vector<vector<int>> &Conn,
                         vector<pair<int,string>> &GmshPhyGroup)
{
    if(!HasMshFileName)
    {
        Msg_Gmsh_NoMshFile();
        return false;
    }
    ifstream in;

    in.open(MshFileName.c_str(),ios::in);
    string line,str,str0;

    int phydim,phyid;
    string phyname;

    ElmtMaxDim=0;

    while (!in.eof())
    {
        getline(in,line);
        if(line.find("$PhysicalNames")!=string::npos)
        {
            nPhysics=0;
            MaxPhyDim=-10;MinPhyDim=10;
            in>>nPhysics;
            getline(in,line);
            GmshPhyGroup.clear();
            GmshPhyGroup.resize(nPhysics);
            for(int i=0;i<nPhysics;i++)
            {
                //phyical dimension,physical id,physical name
                getline(in,line);
                istringstream s_stream(line);
                s_stream >> phydim >> phyid >> phyname;

                if(phydim>MaxPhyDim) MaxPhyDim=phydim;
                if(phydim<MinPhyDim) MinPhyDim=phydim;

                //remove ""
                phyname.erase(remove(phyname.begin(),phyname.end(),'"'), phyname.end());
                GmshPhyGroup[phyid-1]=make_pair(phydim,phyname);

            }

            getline(in,line);// read '$EndPhysicalNames'
        }
        else if(line.find("$NOD")!=string::npos ||
                line.find("$NOE")!=string::npos||
                line.find("$Nodes")!=string::npos)
        {
            // Read node's coordinates
            nNodes=0;
            in>>nNodes;
            NodeCoords.clear();
            NodeCoords.resize(4*nNodes,0.0);

            int nodeid;
            double x,y,z;

            for(int i=0;i<nNodes;i++)
            {
                in>>nodeid>>x>>y>>z;

                NodeCoords[(nodeid-1)*4+0]=1.0;
                NodeCoords[(nodeid-1)*4+1]=x;
                NodeCoords[(nodeid-1)*4+2]=y;
                NodeCoords[(nodeid-1)*4+3]=z;
            }

            getline(in,line);// read '$EndNodes'
        }
        else if(line.find("$Elements")!=string::npos||
                line.find("$ELM")!=string::npos)
        {
            in>>nElmts;

            if(nElmts<1) return false;

            Conn.clear();
            //elm-number elm-type number-of-tags .. node-list
            //                    first  tag: physical id
            //                    second tag: geometry id

            int elmtid,geoid,dim,elmttype,ntags,nodeid;
            int nNodesPerElmt;
            ElmtMaxDim=0;
            vector<int> tempvec;
            for(int e=0;e<nElmts;e++)
            {
                in>>elmtid>>elmttype>>ntags;
                nNodesPerElmt=GetNodesNumViaElmtType(elmttype);
                dim=GetElmtDimViaElmtType(elmttype);
                if(dim>ElmtMaxDim) ElmtMaxDim=dim;


                //Conn[e+1].push_back(nNodesPerElmt);
                tempvec.clear();
                tempvec.push_back(nNodesPerElmt);
                in>>phyid>>geoid;
                for(int i=0;i<nNodesPerElmt;i++)
                {
                    in>>nodeid;
                    tempvec.push_back(nodeid);
                }
                Conn.push_back(tempvec);
            }
            getline(in,line);
        }
    }
}

