#include <iostream>
#include <vector>
#include "Eigen/Eigen"



#include "Welcome.h"

#include "InputSystem/InputSystem.h"
#include "MsgPrint/MsgPrintForInput.h"

#include "Mesh/Mesh.h"

using namespace std;
using namespace Eigen;

int main(int argc,char *argv[])
{
    Welcome(1.0);

    Mesh mesh;

    mesh.SetMshFileName("rect.msh");
    mesh.ReadMesh("gmsh");
    mesh.PrintMeshInfo();

    mesh.SaveMeshAsVTU();

    return 0;
}
