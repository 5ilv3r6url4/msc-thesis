LIBS += -L/usr/lib/lapack/
LIBS += -L/usr/local/lib
LIBS += -L/usr/lib

INCLUDEPATH += /usr/local/include/eigen3
INCLUDEPATH += /usr/include/superlu
INCLUDEPATH += /usr/include/suitesparse

LIBS += -lQGLViewer
LIBS += -lGLU
LIBS += -lsuperlu
LIBS += -lcholmod
