#include "gui/main_interface.h"

Main_Interface::Main_Interface(QWidget *parent)
    : QMainWindow(parent)
{
    setupUi(this);

    m_editor = new Editor();
    m_canvas = new Canvas(m_editor, this->widget_canvas_area);

    m_widget_selection_info = new Widget_Selection_Info();
    m_widget_selection_info->clear_all_forms_info();
    this->widget_central->layout()->addWidget(m_widget_selection_info);

    connect(m_editor,                SIGNAL(broadcast_polybezier_a_info(polybezier_info_form_s&)),
            m_widget_selection_info, SLOT(write_polybez_a_form_info(polybezier_info_form_s&)));
    connect(m_editor,                SIGNAL(broadcast_polybezier_b_info(polybezier_info_form_s&)),
            m_widget_selection_info, SLOT(write_polybez_b_form_info(polybezier_info_form_s&)));
    connect(m_editor,                SIGNAL(broadcast_intersection_info(intersection_info_form_s&)),
            m_widget_selection_info, SLOT(write_intersection_form_info(intersection_info_form_s&)));
    connect(m_editor,                SIGNAL(signal_deselect()),
            m_widget_selection_info, SLOT(clear_all_forms_info()));

    connect(m_widget_selection_info, SIGNAL(signal_spnbox_polybez_a_segment(int)),
            m_editor,                SLOT(response_polybez_a_segment(int)));
    connect(m_widget_selection_info, SIGNAL(signal_spnbox_polybez_b_segment(int)),
            m_editor,                SLOT(response_polybez_b_segment(int)));

    connect(m_widget_selection_info, SIGNAL(request_selected_polybez_curve_type_update(curve_type_t)),
            m_editor,                SLOT(response_selected_polybez_curve_type_update(curve_type_t)));
    connect(m_widget_selection_info, SIGNAL(request_selected_polybez_plane_normal_update(vector_t)),
            m_editor,                SLOT(response_selected_polybez_plane_normal_update(vector_t)));
    connect(m_widget_selection_info, SIGNAL(request_selected_polybez_plane_offset_update(double)),
            m_editor,                SLOT(response_selected_polybez_plane_offset_update(double)));
    connect(m_widget_selection_info, SIGNAL(request_selected_polybez_ctrl_point_update(int, int, point_t)),
            m_editor,                SLOT(response_selected_polybez_ctrl_point_update(int, int, point_t)));

    m_dialog_new_file = new Dialog_New_File();
    m_dialog_new_file->hide();

    m_signal_mapper_toolbox = new QSignalMapper(this);
    m_signal_mapper_canvas_flags = new QSignalMapper(this);
    m_signal_mapper_file_actions = new QSignalMapper(this);

    QAction *act_draw_cross = new QAction("cross", this);
    act_draw_cross->setCheckable(true);
    QAction *act_draw_silhouette = new QAction("silhouette", this);
    act_draw_silhouette->setCheckable(true);
    QAction *act_draw_boundary = new QAction("boundary", this);
    act_draw_boundary->setCheckable(true);

    QActionGroup *group_palette = new QActionGroup(this);
    group_palette->setExclusive(true);
    group_palette->addAction(act_draw_cross);
    group_palette->addAction(act_draw_silhouette);
    group_palette->addAction(act_draw_boundary);

    QMenu *menu_palette = new QMenu(this);
    menu_palette->addAction(act_draw_cross);
    menu_palette->addAction(act_draw_silhouette);
    menu_palette->addAction(act_draw_boundary);
    this->btn_palette->setMenu(menu_palette);
    act_draw_cross->setChecked(true);

    connect(menu_palette, SIGNAL(triggered(QAction*)), btn_palette, SLOT(setCurrentAction(QAction*)));

    connect(act_draw_cross,             SIGNAL(triggered()), m_signal_mapper_toolbox, SLOT(map()));
    connect(act_draw_silhouette,        SIGNAL(triggered()), m_signal_mapper_toolbox, SLOT(map()));
    connect(act_draw_boundary,          SIGNAL(triggered()), m_signal_mapper_toolbox, SLOT(map()));
    connect(this->btn_new_file,         SIGNAL(clicked()), m_signal_mapper_toolbox, SLOT(map()));
    connect(this->btn_save_xml,         SIGNAL(clicked()), m_signal_mapper_toolbox, SLOT(map()));
    connect(this->btn_snapshot,         SIGNAL(clicked()), m_signal_mapper_toolbox, SLOT(map()));
    connect(this->btn_new_layer,        SIGNAL(clicked()), m_signal_mapper_toolbox, SLOT(map()));
    connect(this->btn_draw,             SIGNAL(clicked()), m_signal_mapper_toolbox, SLOT(map()));
    connect(this->btn_draw_polygon,     SIGNAL(clicked()), m_signal_mapper_toolbox, SLOT(map()));
    connect(this->btn_select,           SIGNAL(clicked()), m_signal_mapper_toolbox, SLOT(map()));
    connect(this->btn_flip_normal,      SIGNAL(clicked()), m_signal_mapper_toolbox, SLOT(map()));
    connect(this->btn_flop_normal,      SIGNAL(clicked()), m_signal_mapper_toolbox, SLOT(map()));
    connect(this->btn_marquee_rect,     SIGNAL(clicked()), m_signal_mapper_toolbox, SLOT(map()));
    connect(this->btn_inspector,        SIGNAL(clicked()), m_signal_mapper_toolbox, SLOT(map()));
    connect(this->btn_add_intersect,    SIGNAL(clicked()), m_signal_mapper_toolbox, SLOT(map()));
    connect(this->btn_remove_intersect, SIGNAL(clicked()), m_signal_mapper_toolbox, SLOT(map()));

    m_signal_mapper_toolbox->setMapping(this->btn_draw,         0);
    m_signal_mapper_toolbox->setMapping(this->btn_select,       1);
    m_signal_mapper_toolbox->setMapping(this->btn_flip_normal,  2);
    m_signal_mapper_toolbox->setMapping(act_draw_cross,         3);
    m_signal_mapper_toolbox->setMapping(act_draw_silhouette,    4);
    m_signal_mapper_toolbox->setMapping(act_draw_boundary,      5);
    m_signal_mapper_toolbox->setMapping(this->btn_new_file,     6);
    m_signal_mapper_toolbox->setMapping(this->btn_remove_intersect, 7);
    m_signal_mapper_toolbox->setMapping(this->btn_inspector, 8);
    m_signal_mapper_toolbox->setMapping(this->btn_flip_normal, 9);
    m_signal_mapper_toolbox->setMapping(this->btn_flop_normal, 10);
    m_signal_mapper_toolbox->setMapping(this->btn_draw_polygon, 11);
    m_signal_mapper_toolbox->setMapping(this->btn_snapshot, 12);

    connect(m_signal_mapper_toolbox, SIGNAL(mapped(const int &)), this, SLOT(slot_set_canvas_tool(int)));

    connect(this->chkbox_show_curves,       SIGNAL(clicked()), m_signal_mapper_canvas_flags, SLOT(map()));
    connect(this->chkbox_show_planes,       SIGNAL(clicked()), m_signal_mapper_canvas_flags, SLOT(map()));
    connect(this->chkbox_show_normals,      SIGNAL(clicked()), m_signal_mapper_canvas_flags, SLOT(map()));
    connect(this->chkbox_show_fields,       SIGNAL(clicked()), m_signal_mapper_canvas_flags, SLOT(map()));
    connect(this->chkbox_show_labels,       SIGNAL(clicked()), m_signal_mapper_canvas_flags, SLOT(map()));
    connect(this->chkbox_show_all_layers,   SIGNAL(clicked()), m_signal_mapper_canvas_flags, SLOT(map()));
    connect(this->chkbox_show_patches,      SIGNAL(clicked()), m_signal_mapper_canvas_flags, SLOT(map()));
    connect(this->chkbox_show_isecs,        SIGNAL(clicked()), m_signal_mapper_canvas_flags, SLOT(map()));
    connect(this->radbtn_show_2D,           SIGNAL(clicked()), m_signal_mapper_canvas_flags, SLOT(map()));
    connect(this->radbtn_show_3D,           SIGNAL(clicked()), m_signal_mapper_canvas_flags, SLOT(map()));

    m_signal_mapper_canvas_flags->setMapping(this->chkbox_show_curves,       0);
    m_signal_mapper_canvas_flags->setMapping(this->chkbox_show_planes,       1);
    m_signal_mapper_canvas_flags->setMapping(this->chkbox_show_normals,      2);
    m_signal_mapper_canvas_flags->setMapping(this->chkbox_show_fields,       3);
    m_signal_mapper_canvas_flags->setMapping(this->chkbox_show_labels,       4);
    m_signal_mapper_canvas_flags->setMapping(this->chkbox_show_all_layers,   5);
    m_signal_mapper_canvas_flags->setMapping(this->chkbox_show_patches,      6);
    m_signal_mapper_canvas_flags->setMapping(this->chkbox_show_isecs        ,7);
    m_signal_mapper_canvas_flags->setMapping(this->radbtn_show_2D,           8);
    m_signal_mapper_canvas_flags->setMapping(this->radbtn_show_3D,           9);

    connect(m_signal_mapper_canvas_flags, SIGNAL(mapped(const int &)), this, SLOT(slot_set_canvas_flag(int)));

    connect(this->act_save_xml,         SIGNAL(triggered()), m_signal_mapper_file_actions, SLOT(map()));
    connect(this->act_load_xml,         SIGNAL(triggered()), m_signal_mapper_file_actions, SLOT(map()));
    connect(this->act_load_maya,        SIGNAL(triggered()), m_signal_mapper_file_actions, SLOT(map()));
    connect(this->act_export_obj,       SIGNAL(triggered()), m_signal_mapper_file_actions, SLOT(map()));
    connect(this->act_load_matlab,      SIGNAL(triggered()), m_signal_mapper_file_actions, SLOT(map()));

    m_signal_mapper_file_actions->setMapping(this->act_save_xml,    0);
    m_signal_mapper_file_actions->setMapping(this->act_load_xml,    1);
    m_signal_mapper_file_actions->setMapping(this->act_load_maya,   2);
    m_signal_mapper_file_actions->setMapping(this->act_export_obj,  3);
    m_signal_mapper_file_actions->setMapping(this->act_load_matlab, 4);

    connect(m_signal_mapper_file_actions, SIGNAL(mapped(const int &)), this, SLOT(slot_file_actions(int)));

    connect(this->spnbox_layers, SIGNAL(valueChanged(int)),     this, SLOT(slot_set_active_layer(int)));
    connect(this->btn_nlopt_solve, SIGNAL(clicked()),           m_editor, SLOT(slot_nlopt_solve()));
    connect(this->btn_nlopt_solve_all, SIGNAL(clicked()),       m_editor, SLOT(slot_nlopt_solve_all()));
    connect(this->btn_gen_report, SIGNAL(clicked()),            this, SLOT(slot_generate_report()));

    connect(m_canvas, SIGNAL(request_process_mouse_event(mouse_action_t, mouse_button_t, point_t, point_t)),
            m_editor, SLOT( slot_process_mouse_event(mouse_action_t, mouse_button_t, point_t, point_t)));

    connect(m_dialog_new_file,  SIGNAL(request_new_file(QString, int, int)),
            m_editor,           SLOT(slot_new_file(QString, int, int)));

    connect(m_canvas,  SIGNAL(request_adjust_size(int, int)),
            this,      SLOT(slot_adjust_size(int, int)));
}

Main_Interface::~Main_Interface()
{
    delete m_canvas;
    delete m_signal_mapper_toolbox;
    delete m_signal_mapper_canvas_flags;
    delete m_signal_mapper_file_actions;
    delete m_dialog_new_file;
}

void Main_Interface::slot_update_canvas() { m_canvas->update(); }
void Main_Interface::slot_adjust_size(int width, int height) { QMainWindow::resize(width, height); }

void Main_Interface::slot_set_canvas_tool(int mode)
{
    switch (mode) {
    // TODO setters
    case 0: { m_editor->m_active_tool = TOOL_DRAW_FREE; break; }
    case 1: { m_editor->m_active_tool = TOOL_SELECT; break; }
    case 2: { m_editor->m_active_tool = TOOL_INTERSECTION_NORMAL_FLIP; break; }
    case 3: { m_editor->m_curve_type = CURVE_CROSS; break; }
    case 4: { m_editor->m_curve_type = CURVE_SILHOUETTE; break; }
    case 5: { m_editor->m_curve_type = CURVE_BOUNDARY; break; }
    case 6: { m_dialog_new_file->show(); break; }
    case 7: { m_editor->m_active_tool = TOOL_REMOVE_INTERSECTION; break; }
    case 8: { m_editor->m_active_tool = TOOL_INSPECTOR; break; }
    case 9: { m_editor->m_active_tool = TOOL_INTERSECTION_NORMAL_FLIP; break; }
    case 10: { m_editor->m_active_tool = TOOL_INTERSECTION_NORMAL_FLOP; break; }
    case 11: { m_editor->m_active_tool = TOOL_DRAW_POLYGON; break; }
    case 12: { m_canvas->saveSnapshot(false, true); break; }
    default: { break; }
    }
}

void Main_Interface::slot_set_canvas_flag(int flag)
{
    switch (flag) {
    case 0: { m_canvas->set_draw_flag_show_curves(this->chkbox_show_curves->checkState()); break; }
    case 1: { m_canvas->set_draw_flag_show_planes(this->chkbox_show_planes->checkState()); break; }
    case 2: { m_canvas->set_draw_flag_show_normals(this->chkbox_show_normals->checkState()); break; }
    case 3: { m_canvas->set_draw_flag_show_fields(this->chkbox_show_fields->checkState()); break; }
    case 4: { m_canvas->set_draw_flag_show_labels(this->chkbox_show_labels->checkState()); break; }
    case 5: { m_canvas->set_draw_flag_show_all_layers(this->chkbox_show_all_layers->checkState()); break; }
    case 6: { m_canvas->set_draw_flag_show_patches(this->chkbox_show_patches->checkState()); break; }
    case 7: { m_canvas->set_draw_flag_show_isecs(this->chkbox_show_isecs->checkState()); break; }
    case 8: { m_canvas->set_canvas_display_mode(DISPLAY_2D); break; }
    case 9: { m_canvas->set_canvas_display_mode(DISPLAY_3D); break; }
    default: { break; }
    }
    m_canvas->update();
}

void Main_Interface::slot_file_actions(int act)
{
    switch (act) {
    case 0: {
        QString xml_filename = QFileDialog::getSaveFileName(this, "Choose a filename to save under", "../matsketching/paper/examples/", "XML (*.xml)");
        if (xml_filename.isEmpty()) {
            std::cout << "[ ERROR ] unable to save xml." << std::endl;
            return;
        }
        m_editor->save_file_xml(xml_filename);
        break;
    }
    case 1: {
        QString xml_filepath = QFileDialog::getOpenFileName(this, "Select an xml file", "../matsketching/paper/examples/", "Curve Networks (*.xml)");
        if (!m_editor->load_file_xml(xml_filepath)) {
            std::cout << "[ ERROR ] unable to load xml." << std::endl;
            return;
        }
        this->chkbox_show_normals->setChecked(false);
        break;
    }
    case 2: {
        QString maya_filepath = QFileDialog::getOpenFileName(this, "Select a maya file", "../matsketching/paper/examples/", "Splines (*.xml)");
        if (!m_editor->load_file_maya(maya_filepath)) {
            std::cout << "[ ERROR ] unable to load maya." << std::endl;
            return;
        }
        break;
    }
    case 3: { m_editor->export_obj(); break; }
    default: { break; }
    }
}


void Main_Interface::slot_set_active_layer(int layer) { m_editor->set_active_layer(layer); }

void Main_Interface::slot_generate_report()
{
    QFont font_header = QFont("QFont::SansSerif", 5, QFont::Bold);
    QFont font_body = QFont("QFont::SansSerif", 3);
    QFont font_body_small = QFont("QFont::SansSerif", 2);

    QString path_debug = "../../report/";
    const QString summary_debug_filepath(path_debug + QDateTime::currentDateTime().toString("yyyy-MM-dd") + ".pdf");

    QPdfWriter pdf_writer(summary_debug_filepath);

    QPageSize pdf_writer_pagesize = QPageSize(QSizeF(2, 6), QPageSize::Inch, QString(), QPageSize::SizeMatchPolicy::ExactMatch);
    QPageLayout::Orientation pdf_writer_orientation = QPageLayout::Orientation::Portrait;
    QMarginsF pdf_writer_margins = QMarginsF(0, 0, 0, 0);
    QPageLayout pdf_writer_layout = QPageLayout(pdf_writer_pagesize, pdf_writer_orientation, pdf_writer_margins);

    pdf_writer.setPageLayout(pdf_writer_layout);
    pdf_writer.setResolution(1000);

    QSize pdf_pagesize = pdf_writer_pagesize.sizePixels(1000);
    QPainter painter(&pdf_writer);

    QRectF box_img_nlopt;
    QRectF box_img_matlab;
    QRectF box_nlopt_summary = QRectF(0, pdf_pagesize.height()/4, pdf_pagesize.width()/2, 2 * pdf_pagesize.height()/3);
    QRectF box_matlab_summary = QRectF(pdf_pagesize.width()/2, pdf_pagesize.height()/4, pdf_pagesize.width()/2, 2 * pdf_pagesize.height()/3);

    QRectF box_layer_img_nlopt;
    QRectF box_layer_img_matlab;
    QRectF box_layer_nlopt_summary = QRectF(0, pdf_pagesize.height()/4, pdf_pagesize.width()/2, 2 * pdf_pagesize.height()/3);
    QRectF box_layer_matlab_summary = QRectF(pdf_pagesize.width()/2, pdf_pagesize.height()/4, pdf_pagesize.width()/2, 2 * pdf_pagesize.height()/3);

    QDir directory_debug(path_debug);
    QStringList debug_subdirectory_names = directory_debug.entryList(QDir::AllDirs | QDir::NoDotAndDotDot);

    foreach(QString debug_subdirectory_name, debug_subdirectory_names) {

        const QString debug_filepath(path_debug + debug_subdirectory_name + "-" + QDateTime::currentDateTime().toString("yyyy-MM-dd") + ".pdf");
        QPdfWriter pdf_writer_layers(debug_filepath);
        pdf_writer_layers.setPageLayout(pdf_writer_layout);
        pdf_writer_layers.setResolution(1000);
        QPainter painter_layers(&pdf_writer_layers);

        QString debug_best_guesses_nlopt_filepath = path_debug + debug_subdirectory_name + "/" + debug_subdirectory_name + "-NLOPT_BEST_GUESSES.txt";
        QString debug_best_guesses_matlab_filepath = path_debug + debug_subdirectory_name + "/" + debug_subdirectory_name + "-MATLAB_BEST_GUESSES.txt";

        QFile debug_best_guesses_nlopt_file(debug_best_guesses_nlopt_filepath);
        QFile debug_best_guesses_matlab_file(debug_best_guesses_matlab_filepath);

        if (debug_best_guesses_nlopt_file.open(QFile::ReadOnly | QFile::Text) && debug_best_guesses_matlab_file.open(QFile::ReadOnly | QFile::Text)) {
            QTextStream summary_layer_nlopt_filestream(&debug_best_guesses_nlopt_file);
            QTextStream summary_layer_matlab_filestream(&debug_best_guesses_matlab_file);
            QString best_guesses_nlopt_str = summary_layer_nlopt_filestream.readAll();
            QString best_guesses_matlab_str = summary_layer_matlab_filestream.readAll();
            QStringList best_guesses_nlopt = best_guesses_nlopt_str.split("\n");
            QStringList best_guesses_matlab = best_guesses_matlab_str.split("\n");
            for (int i = 0; i < best_guesses_nlopt.size() - 1; ++i) {
                QString debug_file_layers_nlopt = debug_subdirectory_name + "-" + "L" + best_guesses_nlopt[i].split("-")[0] + "-G" + best_guesses_nlopt[i].split("-")[1];
                QString debug_file_layers_matlab = debug_subdirectory_name + "-" + "L" + best_guesses_matlab[i].split("-")[0] + "-G" + best_guesses_matlab[i].split("-")[1];
                QString nlopt_img_path = debug_subdirectory_name + "-" + "L" + best_guesses_nlopt[i].split("-")[0] + "-G" + best_guesses_nlopt[i].split("-")[1] + "-NLOT_SNAPSHOT.jpg";
                QString matlab_img_path = debug_subdirectory_name + "-" + "L" + best_guesses_matlab[i].split("-")[0] + "-G" + best_guesses_matlab[i].split("-")[1] + "-MATLAB.jpg";

                // ------------------
                // IMAGE PRINTING

                QString img_nlopt_layer_filepath = path_debug + debug_subdirectory_name + "/" + nlopt_img_path;
                QString img_matlab_layer_filepath = path_debug + debug_subdirectory_name + "/" + matlab_img_path;
                QImage img_layer_nlopt(img_nlopt_layer_filepath);
                QImage img_layer_matlab(img_matlab_layer_filepath);

                QSize img_layer_nlopt_size = img_layer_nlopt.size();
                QSize img_layer_matlab_size = img_layer_matlab.size();

                if (img_layer_nlopt_size != img_layer_matlab_size) {
                    std::cout << ">> nlopt and matlab sizes don't match for ddebug image [" + debug_subdirectory_name.toStdString() + "]" << std::endl;        }

                QSize img_layer_size = img_layer_nlopt_size.width() > img_layer_matlab_size.width()? img_layer_nlopt_size : img_layer_matlab_size;
                if (img_layer_size.width() > img_layer_size.height()) {
                    img_layer_size.scale(pdf_pagesize.width()/2,0, Qt::AspectRatioMode::KeepAspectRatioByExpanding);
                    if (img_layer_size.height() > pdf_pagesize.height()/3) {
                        img_layer_size.scale(0, pdf_pagesize.height()/3, Qt::AspectRatioMode::KeepAspectRatioByExpanding);
                    }
                    int y_offset = (pdf_pagesize.height()/3 - img_layer_size.height())/2;
                    box_layer_img_nlopt = QRectF(0, y_offset, img_layer_size.width(), img_layer_size.height());
                    box_layer_img_matlab = QRectF(pdf_pagesize.width()/2, y_offset, img_layer_size.width(), img_layer_size.height());
                }
                else {
                    img_layer_size.scale(pdf_pagesize.height()/3,0, Qt::AspectRatioMode::KeepAspectRatioByExpanding);
                    if (img_layer_size.width() > pdf_pagesize.width()/2) {
                        img_layer_size.scale(pdf_pagesize.width()/2, 0, Qt::AspectRatioMode::KeepAspectRatioByExpanding);
                    }
                    int x_offset = (pdf_pagesize.width()/2 - img_layer_size.width())/2;
                    box_layer_img_nlopt = QRectF(x_offset, 150, img_layer_size.width(), img_layer_size.height());
                    box_layer_img_matlab = QRectF(x_offset + pdf_pagesize.width()/2, 150, img_layer_size.width(), img_layer_size.height());
                }

                painter_layers.drawImage(box_layer_img_nlopt, img_layer_nlopt);
                painter_layers.drawImage(box_layer_img_matlab, img_layer_matlab);

                // ------------------
                // TEXT PRINTING
                painter_layers.setFont(font_body);
                painter_layers.drawText(20, 120, "NLOPT-" + debug_file_layers_nlopt);
                painter_layers.drawText(pdf_pagesize.width()/2, 120, "MATLAB-" + debug_file_layers_matlab);

                QFile summary_layer_nlopt_file( path_debug + debug_subdirectory_name + "/_" + debug_file_layers_nlopt + "-NLOPT_SUMMARY.txt");
                QFile summary_layer_matlab_file( path_debug + debug_subdirectory_name + "/_" + debug_file_layers_matlab+ "-MATLAB_SUMMARY.txt");
                if (summary_layer_nlopt_file.open(QFile::ReadOnly | QFile::Text) && summary_layer_matlab_file.open(QFile::ReadOnly | QFile::Text)) {
                    QTextStream summary_layer_nlopt_filestream(&summary_layer_nlopt_file);
                    QTextStream summary_layer_matlab_filestream(&summary_layer_matlab_file);
                    painter_layers.setFont(font_body_small);
                    painter_layers.drawText(box_layer_matlab_summary, summary_layer_matlab_filestream.readAll());
                    painter_layers.drawText(box_layer_nlopt_summary, summary_layer_nlopt_filestream.readAll());
                }
                if (i >= best_guesses_nlopt.size() - 1) break;
                pdf_writer_layers.newPage();
            }
        }


        // ------------------
        // IMAGE PRINTING

        QString img_nlopt_filepath = path_debug + debug_subdirectory_name + "/" + debug_subdirectory_name + "-NLOPT_ALL.jpg";
        QString img_matlab_filepath = path_debug + debug_subdirectory_name + "/" + debug_subdirectory_name + "-MATLAB_ALL.jpg";
        QImage img_nlopt(img_nlopt_filepath);
        QImage img_matlab(img_matlab_filepath);

        QSize img_nlopt_size = img_nlopt.size();
        QSize img_matlab_size = img_matlab.size();

        if (img_nlopt_size != img_matlab_size) {
            std::cout << ">> nlopt and matlab sizes don't match for debug image [" + debug_subdirectory_name.toStdString() + "]" << std::endl;        }

        QSize img_size = img_nlopt_size.width() > img_matlab_size.width()? img_nlopt_size : img_matlab_size;
        if (img_size.width() > img_size.height()) {
            img_size.scale(pdf_pagesize.width()/2,0, Qt::AspectRatioMode::KeepAspectRatioByExpanding);
            if (img_size.height() > pdf_pagesize.height()/3) {
                img_size.scale(0, pdf_pagesize.height()/3, Qt::AspectRatioMode::KeepAspectRatioByExpanding);
            }
            int y_offset = (pdf_pagesize.height()/3 - img_size.height())/2;
            box_img_nlopt = QRectF(0, y_offset, img_size.width(), img_size.height());
            box_img_matlab = QRectF(pdf_pagesize.width()/2, y_offset, img_size.width(), img_size.height());
        }
        else {
            img_size.scale(pdf_pagesize.height()/3,0, Qt::AspectRatioMode::KeepAspectRatioByExpanding);
            if (img_size.width() > pdf_pagesize.width()/2) {
                img_size.scale(pdf_pagesize.width()/2, 0, Qt::AspectRatioMode::KeepAspectRatioByExpanding);
            }
            int x_offset = (pdf_pagesize.width()/2 - img_size.width())/2;
            box_img_nlopt = QRectF(x_offset, 150, img_size.width(), img_size.height());
            box_img_matlab = QRectF(x_offset + pdf_pagesize.width()/2, 150, img_size.width(), img_size.height());
        }

        painter.drawImage(box_img_nlopt, img_nlopt);
        painter.drawImage(box_img_matlab, img_matlab);

        // ------------------
        // TEXT PRINTING
        painter.setFont(font_header);
        painter.drawText(20, 120, debug_subdirectory_name);

        QFile summary_nlopt_file( path_debug + debug_subdirectory_name + "/" + debug_subdirectory_name + "-NLOPT_SUMMARY.txt");
        QFile summary_matlab_file( path_debug + debug_subdirectory_name + "/" + debug_subdirectory_name + "-MATLAB_SUMMARY.txt");
        if (summary_nlopt_file.open(QFile::ReadOnly | QFile::Text) && summary_matlab_file.open(QFile::ReadOnly | QFile::Text)) {
            QTextStream summary_nlopt_filestream(&summary_nlopt_file);
            QTextStream summary_matlab_filestream(&summary_matlab_file);
            painter.setFont(font_body);
            painter.drawText(box_matlab_summary, summary_matlab_filestream.readAll());
            painter.drawText(box_nlopt_summary, summary_nlopt_filestream.readAll());
        }
        if (debug_subdirectory_name == "TOOTHPASTE") break;
        pdf_writer.newPage();
    }
}
