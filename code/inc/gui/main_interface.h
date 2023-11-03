#ifndef INTERFACE_H
#define INTERFACE_H

#include "ui_main_interface.h"

#include "gui/canvas.h"
#include "gui/dialog_new_file.h"
#include "gui/widget_selection_info.h"

#include "controller/editor.h"

#include <QSignalMapper>
#include <QFileDialog>
#include <QPdfWriter>
#include <QDateTime>
#include <QLayout>
#include <QDialogButtonBox>
#include <QLineEdit>

class Main_Interface : public QMainWindow,
                       public Ui::Main_Interface
{
    Q_OBJECT

public slots:
    void slot_update_canvas();
    void slot_adjust_size(int width, int height);

    void slot_set_canvas_tool(int mode);
    void slot_set_canvas_flag(int flag);
    void slot_file_actions(int act);

    void slot_set_active_layer(int layer);

    void slot_generate_report();

signals:
    void request_new_file(QString, int, int);

public:
    Main_Interface(QWidget *parent = 0);
    ~Main_Interface();

private:
    void dialog_box_new_file();

    Canvas* m_canvas;
    Editor* m_editor;

    Widget_Selection_Info* m_widget_selection_info;
    Dialog_New_File* m_dialog_new_file;

    QSignalMapper* m_signal_mapper_toolbox;
    QSignalMapper* m_signal_mapper_canvas_flags;
    QSignalMapper* m_signal_mapper_file_actions;
};

#endif //INTERFACE_H
