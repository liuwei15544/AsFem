#include <iostream>

#include "Eigen/Eigen"
#include "Welcome.h"

#include "MsgPrint/MsgPrintForInput.h"
#include "MsgPrint/MsgPrintForProgram.h"

using namespace std;
using namespace Eigen;

int main(int args,char *argv[])
{
    Welcome(1.0);
    Msg_InputFileNameWrong("test.i");
    Msg_ExitProgram();
    return 0;
}
