LIBS += opengl32.lib glu32.lib

LIBS += -L$$PWD/libs/glew/lib/x32/ -lglew32
INCLUDEPATH += $$PWD/libs/glew/inc

LIBS += -L$$PWD/libs/nlopt/lib/ -lnlopt
INCLUDEPATH += $$PWD/libs/nlopt/inc

CONFIG(release, debug|release): LIBS += -L$$PWD/libs/qglviewer/lib -lQGLViewer2
CONFIG(debug, debug|release): LIBS += -L$$PWD/libs/qglviewer/lib -lQGLViewerd2
INCLUDEPATH += $$PWD/libs/qglviewer/inc

glew_dll.files = $$PWD/libs/glew/bin/x32/glew32.dll
CONFIG(release, debug|release): glew_dll.path = $$PWD/build/release
CONFIG(debug, debug|release): glew_dll.path = $$PWD/build/debug

nlopt_dll.files = $$PWD/libs/nlopt/bin/nlopt.dll
CONFIG(release, debug|release): nlopt_dll.path = $$PWD/build/release
CONFIG(debug, debug|release): nlopt_dll.path = $$PWD/build/debug

CONFIG(release, debug|release): qglviewer_dll.files = $$PWD/libs/qglviewer/bin/QGLViewer2.dll
CONFIG(debug, debug|release): qglviewer_dll.files = $$PWD/libs/qglviewer/bin/QGLViewerd2.dll
CONFIG(release, debug|release): qglviewer_dll.path = $$PWD/build/release
CONFIG(debug, debug|release): qglviewer_dll.path = $$PWD/build/debug

INSTALLS += glew_dll nlopt_dll qglviewer_dll
