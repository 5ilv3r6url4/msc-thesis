QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.12.6

LIBS += -stdlib=libc++
QMAKE_CXXFLAGS += -std=c++11
QMAKE_CXXFLAGS += -stdlib=libc++

INCLUDEPATH += glew/include/

INCLUDEPATH += /usr/local/Cellar/superlu43/4.3_1/include/superlu
INCLUDEPATH += /usr/local/Cellar/eigen/3.3.2/include/eigen3/Eigen/src
INCLUDEPATH += /Library/Frameworks/QGLViewer.framework/Versions/2/Headers/
INCLUDEPATH += /usr/local/include/

LIBS += -L/usr/local/lib/ -lnlopt
LIBS += -L/usr/local/Cellar/superlu43/4.3_1/lib
LIBS += -lsuperlu

QMAKE_LFLAGS += -F/Library/Frameworks/
LIBS += -framework QGLViewer
LIBS += -framework Accelerate
