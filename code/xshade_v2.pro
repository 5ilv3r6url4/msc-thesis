win32 {
    !contains(QMAKE_TARGET.arch, x86_64) {
        include(xshade_v2_win32.pro)
    }
    else {
        include(xshade_v2_win64.pro)
    }
}

unix:!macx {
    include(xshade_v2_linux.pro)
}

macx {
    include(xshade_v2_mac.pro)
}

TARGET = xshade_v2
TEMPLATE = app

QT += core gui xml opengl widgets printsupport

CONFIG -= debug_and_release debug_and_release_target

PRO_PATH = $$PWD

SOURCES += $$PWD/$$files(src/*.cpp, true)
HEADERS += $$PWD/$$files(inc/*.h, true) \
    inc/utils/util_filepath_defs.h
FORMS += $$PWD/$$files(ui/*.ui, true)
RESOURCES += $$PWD/$$files(res/*.qrc, true)

INCLUDEPATH += inc/
INCLUDEPATH += libs/boost/

CONFIG(debug,debug|release): DESTDIR = $$PWD/build/debug
CONFIG(release,debug|release): DESTDIR = $$PWD/build/release
MOC_DIR = $$DESTDIR/moc
OBJECTS_DIR = $$DESTDIR/obj
RCC_DIR = $$DESTDIR/rcc
UI_DIR = $$DESTDIR/uci

DEFINES += QT_DEPRECATED_WARNINGS
