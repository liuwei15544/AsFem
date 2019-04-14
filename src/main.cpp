#include <iostream>
#include "Eigen/Eigen"



#include "Welcome.h"

#include "InputSystem/InputSystem.h"
#include "MsgPrint/MsgPrintForInput.h"


using namespace std;
using namespace Eigen;

int main(int argc,char *argv[])
{
    Welcome(1.0);

    InputSystem inputSystem(argc,argv);

    return 0;
}
