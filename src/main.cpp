#include <iostream>
#include <vector>
#include "Eigen/Eigen"



#include "Welcome.h"

#include "InputSystem/InputSystem.h"
#include "MsgPrint/MsgPrintForInput.h"

#include "Mesh/GmshIO.h"

using namespace std;
using namespace Eigen;

int main(int argc,char *argv[])
{
    Welcome(1.0);

    vector<double> NodeCoords;
    vector<vector<int>> Conn;
    vector<pair<int,string>> PhyGroup;
    map<string,vector<int>> MeshNameSet;
    map<int,vector<int>> MeshIdSet;
    vector<int> ElmtVTKCellType;
    vector<string> ElmtTypeName;

    GmshIO gmshIo("rect.msh");
    gmshIo.ReadMshFile(NodeCoords,Conn,ElmtVTKCellType,ElmtTypeName,MeshNameSet,MeshIdSet,PhyGroup);

    return 0;
}
