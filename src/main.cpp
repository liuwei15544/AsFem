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
    vector<pair<int,string>> GmshPhyGroup;

    GmshIO gmshIo;
    gmshIo.SetMshFileName("rect.msh");
    gmshIo.ReadMshFile(NodeCoords,Conn,GmshPhyGroup);

    for(int i=0;i<GmshPhyGroup.size();i++)
    {
        cout<<GmshPhyGroup[i].first<<"==>"<<GmshPhyGroup[i].second<<endl;
    }

    InputSystem inputSystem(argc,argv);



    return 0;
}
