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
                         vector<int> &ElmtVTKCellType,
                         vector<string> &ElmtTypeName,
                         vector<int> &PhyIDIndex,
                         map<string,vector<int>> &MeshNameSet,
                         map<int,vector<int>> &MeshIdSet,
                         vector<pair<int,string>> &GmshPhyGroup)
{
    if(!HasMshFileName)
    {
        Msg_Gmsh_NoMshFile();
        return false;
    }
    ifstream in;

    in.open(MshFileName.c_str(),ios::in);
    if(!in.is_open())
    {
        Msg_Gmsh_NoMshFile();
        return false;
    }
    string line,str,str0;

    int phydim,phyid;
    string phyname;

    ElmtMaxDim=0;

    BulkElmtTypeName="";
    while (!in.eof())
    {
        getline(in,line);
        if(line.find("$PhysicalNames")!=string::npos)
        {
            nPhysics=0;
            MaxPhyDim=-10;MinPhyDim=10;
            in>>nPhysics;// waring: even nPhysics=1, in elmts section, it can be greater
                         // than 1, so do not resize it, just clear and push bakc!!!
            getline(in,line);
            GmshPhyGroup.clear();
            PhyIDIndex.resize(nPhysics,-1);
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
                //GmshPhyGroup[phyid]=make_pair(phydim,phyname);
                GmshPhyGroup.push_back(make_pair(phydim,phyname));

                MeshNameSet[phyname].clear();
                MeshIdSet[phyid].clear();
                PhyIDIndex[i]=phyid;
            }
            MeshNameSet["all"].clear();
            MeshIdSet[0].clear();
            //GmshPhyGroup[0]=make_pair(MaxPhyDim,"all");

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

                if(x<Xmin) Xmin=x;
                if(x>Xmax) Xmax=x;

                if(y<Ymin) Ymin=y;
                if(y>Ymax) Ymax=y;

                if(z<Zmin) Zmin=z;
                if(z>Zmax) Zmax=z;
            }

            getline(in,line);// read '$EndNodes'
        }
        else if(line.find("$Elements")!=string::npos||
                line.find("$ELM")!=string::npos)
        {
            nElmts=0;
            nVolumeElmts=0;nSurfaceElmts=0;nLineElmts=0;nPointElmts=0;
            nBulkElmts=0;

            ElmtVTKCellType.clear();
            ElmtTypeName.clear();

            in>>nElmts;
            if(nElmts<1) return false;

            ElmtTypeName.resize(nElmts);
            ElmtVTKCellType.resize(nElmts,-1);

            Conn.clear();
            //elm-number elm-type number-of-tags .. node-list
            //                    first  tag: physical id
            //                    second tag: geometry id

            int elmtid,geoid,dim,phydim,elmttype,ntags,nodeid;
            int nNodesPerElmt;
            string elmttypename,phyname;
            int ind;
            int elmtvtktype;
            ElmtMaxDim=0;
            vector<int> tempvec;
            bool IsPhyIDInList;

            for(int e=0;e<nElmts;e++)
            {
                in>>elmtid>>elmttype>>ntags;
                nNodesPerElmt=GetNodesNumViaElmtType(elmttype);
                dim=GetElmtDimViaElmtType(elmttype);
                elmttypename=GetElmtNameViaElmtType(elmttype);
                elmtvtktype=GetVTKCellTypeViaElmtType(elmttype);

                ElmtTypeName[elmtid-1]=elmttypename;
                ElmtVTKCellType[elmtid-1]=elmtvtktype;

                if(dim>ElmtMaxDim) ElmtMaxDim=dim;
                if(dim>MaxPhyDim) MaxPhyDim=dim;
                if(dim<MinPhyDim) MinPhyDim=dim;

                if(dim==0) nPointElmts+=1;
                if(dim==1) nLineElmts+=1;
                if(dim==2) nSurfaceElmts+=1;
                if(dim==3) nVolumeElmts+=1;

                if(dim==MaxPhyDim)
                {
                    BulkElmtTypeName=GetElmtNameViaElmtType(elmttype);
                }


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
                // check if current phyid already is in the  phylist
                // if not exist, add it to the GmshPhyGroup
                IsPhyIDInList=false;
                for(int i=0;i<PhyIDIndex.size();i++)
                {
                    if(phyid==PhyIDIndex[i])
                    {
                        IsPhyIDInList=true;
                        phydim=GmshPhyGroup[PhyIDIndex[i]-1].first;
                        break;
                    }
                }

                if(!IsPhyIDInList)
                {
                    // element's physical id isn't in PhyGroup list
                    // then add it to the list
                    PhyIDIndex.push_back(phyid);
                    phyname=to_string(phyid);
                    //GmshPhyGroup[phyid]=make_pair(dim,phyname);
                    GmshPhyGroup.push_back(make_pair(dim,phyname));
                }

                ind=-1;
                for(int i=0;i<PhyIDIndex.size();i++)
                {
                    if(PhyIDIndex[i]==phyid)
                    {
                        ind=i;
                        break;
                    }
                }

                if(dim>GmshPhyGroup[ind].first) GmshPhyGroup[ind].first=dim;

                MeshIdSet[phyid].push_back(elmtid);

                MeshNameSet[GmshPhyGroup[ind].second].push_back(elmtid);

                if(dim==MaxPhyDim)
                {
                    // store all the bulk element(max dim) to "default" for phyname
                    //                                     to  0        for phyid
                    //MeshIdSet[0].push_back(elmtid);
                    //MeshNameSet["all"].push_back(elmtid);
                    nBulkElmts+=1;
                }
            }
            GmshPhyGroup.push_back(make_pair(MaxPhyDim,"all"));
            PhyIDIndex.push_back(0);
            getline(in,line);
        }
    }
}

