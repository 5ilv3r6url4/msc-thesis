#include "gui/main_interface.h"
#include <QtGui>

int main(int argc, char** argv)
{
    // read command lines arguments
    QApplication application(argc,argv);

    // instantiate the interface
    Main_Interface m_interface;
    m_interface.show();

    // run main loop
    return application.exec();
}
