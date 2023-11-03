#ifndef EDITOR_H
#define EDITOR_H

#include <qglviewer.h>

#include "model/curve_network.h"

#include "gui/widget_selection_info.h"

#include "utils/defs/util_data_defs.h"
#include "utils/defs/util_gui_defs.h"
#include "utils/gui/util_gui_forms.h"

#include <QWidget>
#include <QImage>

class Canvas;

class Editor : public QObject
{
    Q_OBJECT
    friend Canvas;

public slots:
    // TODO refactor slots to adhere to naming convention of tool/menu/file action
    void slot_new_file(QString filename, int image_width, int image_height);
    void slot_nlopt_solve();
    void slot_nlopt_solve_all();
    void slot_graph_cut();
    void slot_process_mouse_event(mouse_action_t event, mouse_button_t button, point_t screen_point, point_t canvas_point);
    void slot_snapshot(QString filepath, bool show_all, bool show_normals, bool show_fields);

    void response_polybez_a_segment(int segment);
    void response_polybez_b_segment(int segment);
    void response_selected_polybez_curve_type_update(curve_type_t curve_type);
    void response_selected_polybez_plane_normal_update(vector_t normal);
    void response_selected_polybez_plane_offset_update(double offset);
    void response_selected_polybez_ctrl_point_update(int segment, int ctrl_point_idx, point_t ctrl_point);
    void response_selected_intersection_position_update(point_t position);
    void response_selected_intersection_tangents_update(boost::optional<vector_t>, boost::optional<vector_t>);

signals:
    void broadcast_polybezier_a_info(polybezier_info_form_s& out_form);
    void broadcast_polybezier_b_info(polybezier_info_form_s& out_form);
    void broadcast_intersection_info(intersection_info_form_s& out_form);
    void signal_deselect();

public:
    Editor();
    ~Editor();

    const Curve_Network* get_network_model(int layer);
    void register_canvas(Canvas* canvas);

    void tool_draw_free(mouse_action_t event, mouse_button_t button, point_t screen_point, point_t canvas_point);
    void tool_draw_polygon(mouse_action_t event, mouse_button_t button, point_t screen_point, point_t canvas_point);
    void tool_select(mouse_action_t event, mouse_button_t button, point_t screen_point, point_t canvas_point);
    void tool_inspector(mouse_action_t event, mouse_button_t button, point_t screen_point, point_t canvas_point);
    void tool_remove_intersection(mouse_action_t event, mouse_button_t button, point_t screen_point, point_t canvas_point);
    void tool_intersection_normal_flip(mouse_action_t event, mouse_button_t button, point_t screen_point, point_t canvas_point);
    void tool_intersection_normal_flop(mouse_action_t event, mouse_button_t button, point_t screen_point, point_t canvas_point);
    void tool_move_ctrl_point();
    bool save_file_xml(QString filename) { return true; };
    bool load_file_xml(QString filepath);
    bool load_file_maya(QString filename){ return true; };

    void export_obj() {}
    void set_active_layer(int layer);

    // ----------------------
    // NETWORKS PROPERTIES
    // ----------------------
    int m_num_layers;
    int m_current_layer;

    // ----------------------
    // FILE PROPERTIES
    // ----------------------
    QString m_filename;

    // ----------------------
    // DRAWING DATA
    // ----------------------
    QVector<point_t> m_canvas_raw_path;
    QVector<point_t> m_canvas_draw_path;

    // ----------------------
    // EDITOR STATES
    // ----------------------
    curve_type_t m_curve_type;
    editor_tool_t m_active_tool;


private:
    // ----------------------
    // EDITOR PROPERTIES
    // ----------------------
    Canvas* m_canvas;
    QVector<Curve_Network*> m_curve_networks;

    polybezier_info_form_s package_polybez_info(polybez_ptr_t polybez_ptr, boost::optional<int> current_segment = boost::none);
    intersection_info_form_s package_intersection_info(isec_ptr_t isec_ptr);

    struct current_selection_t {
        boost::optional<boost::uuids::uuid> _polybez_uuid = boost::none;
        boost::optional<boost::uuids::uuid> _isec_uuid    = boost::none;
        int _ctrl_point_idx = -1;

        void select_polybez(boost::uuids::uuid uuid) { _polybez_uuid = uuid; _isec_uuid = boost::none; _ctrl_point_idx = -1; }
        void select_isec(boost::uuids::uuid uuid) { _polybez_uuid = boost::none; _isec_uuid = uuid; _ctrl_point_idx = -1; }
        void select_ctrl_point(int idx) { _isec_uuid = boost::none; _ctrl_point_idx = idx; }
        void unselect() { _polybez_uuid = boost::none; _isec_uuid = boost::none; _ctrl_point_idx = -1; }
    } m_current_selection;
};

#endif // EDITOR_H
