#include "controller/editor.h"

#include "gui/canvas.h"

#include "utils/math/util_bezier_fit.h"

Editor::Editor()
{    
    m_curve_networks.resize(0);
    m_current_layer = -1;
    m_num_layers = -1;
}

Editor::~Editor() { }

const Curve_Network* Editor::get_network_model(int layer)
{
    return m_curve_networks[layer];
}

void Editor::register_canvas(Canvas *canvas)
{
    m_canvas = canvas;
}

void Editor::response_polybez_a_segment(int segment)
{
    curve_graph_vertex_data_map_t polybez_ptr_map = boost::get(boost::vertex_data_pointer, m_curve_networks[m_current_layer]->m_curve_network_graph);
    if (m_current_selection._isec_uuid == boost::none) {
        if (m_current_selection._polybez_uuid == boost::none) {
            std::cout << "[ ERROR ] no polybezier selected." << std::endl;
            return;
        }
        boost::optional<curve_graph_vertex_desc_t> polybez_desc = Curve_Network::find_vertex_by_uuid(m_current_selection._polybez_uuid.value(), m_curve_networks[m_current_layer]->m_curve_network_graph);
        if (polybez_desc == boost::none) {
            std::cout << "[ ERROR ] attempt to manually edit polybezier that is not selected." << std::endl;
            return;
        }

        polybez_ptr_t polybez_ptr = polybez_ptr_map[polybez_desc.value()];
        polybezier_info_form_s updated_form = package_polybez_info(polybez_ptr, segment);
        emit broadcast_polybezier_a_info(updated_form);
        m_canvas->update();
    }

    else if (m_current_selection._polybez_uuid == boost::none) {
        if (m_current_selection._isec_uuid == boost::none) {
            std::cout << "[ ERROR ] no intersection selected." << std::endl;
            return;
        }
        //curve_graph_edge_data_map_t isec_ptr_map = boost::get(boost::edge_data_pointer, m_curve_networks[m_current_layer]->m_curve_network_graph);
        boost::optional<curve_graph_edge_desc_t> isec_desc = Curve_Network::find_edge_by_uuid(m_current_selection._isec_uuid.value(), m_curve_networks[m_current_layer]->m_curve_network_graph);
        if (isec_desc == boost::none) {
            std::cout << "[ ERROR ] attempt to manually edit intersection that is not selected." << std::endl;
            return;
        }

        //        boost::optional<curve_graph_vertex_desc_t> polybez_desc = Curve_Network::find_vertex_by_uuid(isec_ptr_map[isec_desc.value()]->m_polybez_a_id, m_curve_networks[m_current_layer]->m_curve_network_graph);
        //        if (polybez_desc == boost::none) {
        //            std::cout << "[ ERROR ] attempt to manually edit polybezier that is not selected." << std::endl;
        //            return;
        //        }
        //        polybez_ptr_t polybez_ptr = polybez_ptr_map[polybez_desc.value()];
        polybez_ptr_t polybez_ptr = polybez_ptr_map[isec_desc->m_target];
        polybezier_info_form_s updated_form = package_polybez_info(polybez_ptr, segment);
        emit broadcast_polybezier_a_info(updated_form);
        m_canvas->update();
    }
}

void Editor::response_polybez_b_segment(int segment)
{
    curve_graph_vertex_data_map_t polybez_ptr_map = boost::get(boost::vertex_data_pointer, m_curve_networks[m_current_layer]->m_curve_network_graph);
    if (m_current_selection._isec_uuid == boost::none) {
        std::cout << "[ ERROR ] no intersection selected." << std::endl;
        return;
    }
    curve_graph_edge_data_map_t isec_ptr_map = boost::get(boost::edge_data_pointer, m_curve_networks[m_current_layer]->m_curve_network_graph);
    boost::optional<curve_graph_edge_desc_t> isec_desc = Curve_Network::find_edge_by_uuid(m_current_selection._isec_uuid.value(), m_curve_networks[m_current_layer]->m_curve_network_graph);
    if (isec_desc == boost::none) {
        std::cout << "[ ERROR ] attempt to manually edit intersection that is not selected." << std::endl;
        return;
    }
    //    isec_ptr_t isec_ptr = isec_ptr_map[isec_desc.value()];

    //    boost::optional<curve_graph_vertex_desc_t> polybez_desc = Curve_Network::find_vertex_by_uuid(isec_ptr->m_polybez_b_id, m_curve_networks[m_current_layer]->m_curve_network_graph);
    //    if (polybez_desc == boost::none) {
    //        std::cout << "[ ERROR ] attempt to manually edit polybezier that is not selected." << std::endl;
    //        return;
    //    }
    //    polybez_ptr_t polybez_ptr = polybez_ptr_map[polybez_desc.value()];
    polybez_ptr_t polybez_ptr = polybez_ptr_map[isec_desc->m_source];
    polybezier_info_form_s updated_form = package_polybez_info(polybez_ptr, segment);
    emit broadcast_polybezier_b_info(updated_form);
    m_canvas->update();
}

void Editor::response_selected_polybez_curve_type_update(curve_type_t curve_type)
{
    if (m_current_selection._polybez_uuid == boost::none) {
        std::cout << "[ ERROR ] no polybezier selected." << std::endl;
        return;
    }
    boost::optional<curve_graph_vertex_desc_t> polybez_desc = Curve_Network::find_vertex_by_uuid(m_current_selection._polybez_uuid.value(), m_curve_networks[m_current_layer]->m_curve_network_graph);
    if (polybez_desc == boost::none) {
        std::cout << "[ ERROR ] attempt to manually edit polybezier that is not selected." << std::endl;
        return;
    }
    curve_graph_vertex_data_map_t polybez_ptr_map = boost::get(boost::vertex_data_pointer, m_curve_networks[m_current_layer]->m_curve_network_graph);
    polybez_ptr_t polybez_ptr = polybez_ptr_map[polybez_desc.value()];
    polybez_ptr->m_curve_type = curve_type;

    polybezier_info_form_s updated_form = package_polybez_info(polybez_ptr);
    emit broadcast_polybezier_a_info(updated_form);
    m_canvas->update();
}

// TODO private data with setters and getters?
void Editor::response_selected_polybez_plane_normal_update(vector_t normal)
{
    if (m_current_selection._polybez_uuid == boost::none) {
        std::cout << "[ ERROR ] no polybezier selected." << std::endl;
        return;
    }
    boost::optional<curve_graph_vertex_desc_t> polybez_desc = Curve_Network::find_vertex_by_uuid(m_current_selection._polybez_uuid.value(), m_curve_networks[m_current_layer]->m_curve_network_graph);
    if (polybez_desc == boost::none) {
        std::cout << "[ ERROR ] attempt to manually edit polybezier that is not selected." << std::endl;
        return;
    }
    curve_graph_vertex_data_map_t polybez_ptr_map = boost::get(boost::vertex_data_pointer, m_curve_networks[m_current_layer]->m_curve_network_graph);
    polybez_ptr_t polybez_ptr = polybez_ptr_map[polybez_desc.value()];
    polybez_ptr->m_plane_normal = normal;

    polybezier_info_form_s updated_form = package_polybez_info(polybez_ptr);
    emit broadcast_polybezier_a_info(updated_form);
    m_canvas->update();
}

void Editor::response_selected_polybez_plane_offset_update(double offset)
{
    if (m_current_selection._polybez_uuid == boost::none) {
        std::cout << "[ ERROR ] no polybezier selected." << std::endl;
        return;
    }
    boost::optional<curve_graph_vertex_desc_t> polybez_desc = Curve_Network::find_vertex_by_uuid(m_current_selection._polybez_uuid.value(), m_curve_networks[m_current_layer]->m_curve_network_graph);
    if (polybez_desc == boost::none) {
        std::cout << "[ ERROR ] attempt to manually edit polybezier that is not selected." << std::endl;
        return;
    }
    curve_graph_vertex_data_map_t polybez_ptr_map = boost::get(boost::vertex_data_pointer, m_curve_networks[m_current_layer]->m_curve_network_graph);
    polybez_ptr_t polybez_ptr = polybez_ptr_map[polybez_desc.value()];
    polybez_ptr->m_plane_offset = offset;

    polybezier_info_form_s updated_form = package_polybez_info(polybez_ptr);
    emit broadcast_polybezier_a_info(updated_form);
    m_canvas->update();
}

void Editor::response_selected_polybez_ctrl_point_update(int segment, int ctrl_point_idx, point_t ctrl_point)
{
    if (m_current_selection._polybez_uuid == boost::none) {
        std::cout << "[ ERROR ] no polybezier selected." << std::endl;
        return;
    }
    boost::optional<curve_graph_vertex_desc_t> polybez_desc = Curve_Network::find_vertex_by_uuid(m_current_selection._polybez_uuid.value(), m_curve_networks[m_current_layer]->m_curve_network_graph);
    if (polybez_desc == boost::none) {
        std::cout << "[ ERROR ] attempt to manually edit polybezier that is not selected." << std::endl;
        return;
    }
    curve_graph_vertex_data_map_t polybez_ptr_map = boost::get(boost::vertex_data_pointer, m_curve_networks[m_current_layer]->m_curve_network_graph);
    polybez_ptr_t polybez_ptr = polybez_ptr_map[polybez_desc.value()];

    point_t  old_ctrl_point = polybez_ptr->m_bezier_curves[segment].ctrl_points[ctrl_point_idx];
    vector_t delta_ctrl_point = ctrl_point - old_ctrl_point;

    if ( ctrl_point_idx == 0 && segment == 0 ) {
        polybez_ptr->m_bezier_curves[segment].ctrl_points[0] = ctrl_point;
        polybez_ptr->m_bezier_curves[segment].ctrl_points[1] += delta_ctrl_point;
    }

    else if ( ctrl_point_idx == 1 && segment == 0) {
        polybez_ptr->m_bezier_curves[segment].ctrl_points[1] = ctrl_point;
    }

    else if ( ctrl_point_idx == 3 && segment == polybez_ptr->m_bezier_curves.size() - 1) {
        polybez_ptr->m_bezier_curves[segment].ctrl_points[3] = ctrl_point;
        polybez_ptr->m_bezier_curves[segment].ctrl_points[2] += delta_ctrl_point;
    }

    else if ( ctrl_point_idx == 2 && segment == polybez_ptr->m_bezier_curves.size() - 1) {
        polybez_ptr->m_bezier_curves[segment].ctrl_points[2] = ctrl_point;
    }

    else {
        switch ( ctrl_point_idx )
        {
        case 0: {
            polybez_ptr->m_bezier_curves[segment].ctrl_points[0] = ctrl_point;
            polybez_ptr->m_bezier_curves[segment].ctrl_points[1] += delta_ctrl_point;
            polybez_ptr->m_bezier_curves[segment - 1].ctrl_points[3] = ctrl_point;
            polybez_ptr->m_bezier_curves[segment - 1].ctrl_points[2] += delta_ctrl_point;
            break;
        }
        case 1: {
            polybez_ptr->m_bezier_curves[segment].ctrl_points[1] = ctrl_point;
            polybez_ptr->m_bezier_curves[segment - 1].ctrl_points[2] -= delta_ctrl_point;
            break;
        }
        case 2: {
            polybez_ptr->m_bezier_curves[segment].ctrl_points[2] = ctrl_point;
            polybez_ptr->m_bezier_curves[segment + 1].ctrl_points[1] -= delta_ctrl_point;
            break;
        }
        case 3: {
            polybez_ptr->m_bezier_curves[segment].ctrl_points[3] = ctrl_point;
            polybez_ptr->m_bezier_curves[segment].ctrl_points[2] += delta_ctrl_point;
            polybez_ptr->m_bezier_curves[segment + 1].ctrl_points[0] = ctrl_point;
            polybez_ptr->m_bezier_curves[segment + 1].ctrl_points[1] += delta_ctrl_point;
            break;
        }
        default: {
            break;
        }
        }
    }
    polybezier_info_form_s updated_form = package_polybez_info(polybez_ptr, segment);
    emit broadcast_polybezier_a_info(updated_form);
    m_canvas->update();
}

void Editor::response_selected_intersection_position_update(point_t position) { }
void Editor::response_selected_intersection_tangents_update(boost::optional<vector_t>, boost::optional<vector_t>) { }

polybezier_info_form_s Editor::package_polybez_info(polybez_ptr_t polybez_ptr, boost::optional<int> current_segment)
{
    polybezier_info_form_s out_form;
    out_form.uuid = polybez_ptr->m_uuid;
    out_form.curve_type = polybez_ptr->m_curve_type;
    out_form.plane_normal = polybez_ptr->m_plane_normal;
    out_form.plane_offset = polybez_ptr->m_plane_offset;
    out_form.num_segments = polybez_ptr->m_bezier_curves.size();
    out_form.curr_segment = current_segment.value_or(0);
    out_form.ctrl_points[0] = polybez_ptr->m_bezier_curves[current_segment.value_or(0)].ctrl_points[0];
    out_form.ctrl_points[1] = polybez_ptr->m_bezier_curves[current_segment.value_or(0)].ctrl_points[1];
    out_form.ctrl_points[2] = polybez_ptr->m_bezier_curves[current_segment.value_or(0)].ctrl_points[2];
    out_form.ctrl_points[3] = polybez_ptr->m_bezier_curves[current_segment.value_or(0)].ctrl_points[3];
    return out_form;
}

intersection_info_form_s Editor::package_intersection_info(isec_ptr_t isec_ptr)
{
    intersection_info_form_s out_form;
    out_form.uuid = isec_ptr->m_uuid;
    out_form.uuid_polybez_a = isec_ptr->m_polybez_a_uuid;
    out_form.uuid_polybez_b = isec_ptr->m_polybez_b_uuid;
    out_form.position = isec_ptr->m_position;
    out_form.tval_a = isec_ptr->m_tval_a;
    out_form.tval_b = isec_ptr->m_tval_b;
    out_form.normal = isec_ptr->m_normal;
    out_form.tangent_a = isec_ptr->m_tangent_a;
    out_form.tangent_b = isec_ptr->m_tangent_b;
    return out_form;
}

void Editor::set_active_layer(int layer)
{
    if (layer < m_num_layers) {
        m_current_layer = layer;
    }
    else {
        m_current_layer = m_num_layers - 1;
    }
    m_canvas->update();
}

void Editor::slot_new_file(QString filename, int image_width, int image_height)
{
    // -------------------------------
    // clear the currently stored data
    if (m_curve_networks.size() != 0) {
        qDeleteAll(m_curve_networks);
        m_curve_networks.clear();
    }

    // ---------------------
    // create new empty data
    m_curve_networks.resize(1);
    Curve_Network* curve_network = new Curve_Network();
    m_curve_networks[0] = curve_network;

    m_curve_type = CURVE_CROSS;
    m_current_layer = 0;
    m_num_layers = 1;

    m_filename = filename;

    // -----------------
    // update new canvas
    m_canvas->new_drawing_canvas(image_width, image_height);
    connect(m_curve_networks[0], SIGNAL(request_snapshot(QString, bool, bool, bool)), this, SLOT(slot_snapshot(QString, bool, bool, bool)));
}

bool Editor::load_file_xml(QString filepath)
{
    /*------------------------------------------------------------------
     |
     |  I make the supposition that the xml file is correctly filled-in.
     |  I'm not checking that each node is where it's supposed to be.
     |
     *------------------------------------------------------------------*/

    // -----------------------------------------
    // check that filepath is correct and exists
    if (filepath.isEmpty()) { return false; }
    if (!QFile::exists (filepath)){
        std::cout << "[ ERROR ] file " << filepath.toStdString() <<" does not exist!" << std::endl;
        return false;
    }

    // -------------------------------
    // open in read-only and load data
    QFile xml_file(filepath);
    if (!xml_file.open (QIODevice::ReadOnly)) { return false; }

    QDomDocument doc ("CurveSetXML");
    if (!doc.setContent(&xml_file)) {
        xml_file.close();
        return false;
    }

    // ---------------------
    // save file information
    QFileInfo fileInfo(xml_file);
    QString filename = fileInfo.fileName();
    m_filename = filename.split(".")[0];
    xml_file.close();

    // -------------------------------
    // clear the currently stored data
    if (m_curve_networks.size() != 0) {
        qDeleteAll(m_curve_networks);
        m_curve_networks.clear();
    }

    // ---------------------
    // create new empty data
    m_curve_networks.resize(1);
    m_current_layer = 0;
    m_num_layers = 1;
    Curve_Network* curve_network = new Curve_Network();
    m_curve_networks[0] = curve_network;

    // --------------------------------------------------
    // root node attributes ( image_width, image_height )
    const QDomElement root = doc.documentElement();
    int image_width = (root.hasAttribute("image_width"))? root.attribute("image_width").toInt() : 0;
    int image_height = (root.hasAttribute("image_height"))? root.attribute("image_height").toInt() : 0;

//    double height_factor = 1.0;
//    double width_factor = 1.0;

//    if (image_height > 800) {
//        height_factor = 800.0 / (double)image_height;
//        width_factor = 800.0 * ((double)image_width / (double)image_height) / (double)image_width;
//    }

//    if (image_width > 950) {
//        width_factor = 950.0 / (double)image_width;
//        height_factor = 950.0 * ((double)image_height / (double)image_width) / (double)image_height;
//    }

    // -------------------------------------------------------------------
    // parse all curves, each curve has attribute (layer, curve type, ...)
    // and children (control_points_set)
    //
    // read attributes. For example for an integer XXX
    // splines[selected].setXXX(curve_node.attribute("XXX").toInt());
    for(QDomElement curve_node = root.firstChildElement("curve"); !curve_node.isNull(); curve_node = curve_node.nextSiblingElement()) {

        // ----------------------
        // basic curve components
        bool is_closed_curve;
        int layer;
        QVector<point_t> ctrl_points;
        curve_type_t curve_type;

        // parse layer and resize if needed
        layer = ( curve_node.hasAttribute("layer") ) ? curve_node.attribute("layer").toInt():0;
        if (layer >= m_num_layers) {
            m_num_layers = layer + 1;
            m_curve_networks.resize(m_num_layers);
        }

        if (m_curve_networks[layer] == NULL) {
            m_curve_networks[layer] = new Curve_Network();
        }

        // parse curve type
        int ct = ( curve_node.hasAttribute("curve_type") ) ? curve_node.attribute("curve_type").toInt():0;
        switch (ct)
        {
        case 0:  { curve_type = CURVE_CROSS; break;      }
        case 1:  { curve_type = CURVE_BOUNDARY; break;   }
        case 2:  { curve_type = CURVE_SILHOUETTE; break; }
        default: { curve_type = CURVE_CROSS; break;      }
        }

        // parse control points
        QDomElement ctrl_points_dom_element = curve_node.firstChildElement("control_points_set");
        if(!ctrl_points_dom_element.isNull()) {
            for(QDomElement ctrl_point = ctrl_points_dom_element.firstChildElement("control_point");
                !ctrl_point.isNull();
                ctrl_point = ctrl_point.nextSiblingElement()) {
                double x = 0.0;
                double y = 0.0;
                y = (1.0 - (ctrl_point.attribute("x").toDouble()) / (0.5 * (double)image_height));// * height_factor;
                x = ((ctrl_point.attribute("y").toDouble()) / (0.5 * (double)image_width) - 1.0);// * width_factor;
                ctrl_points.push_back(point_t(x, y));
            }
        }

        // calculate closed curve boolean
        double dist = (ctrl_points[0] - ctrl_points[ctrl_points.size() - 1]).L2_norm2D();
        is_closed_curve = (dist < 1.0e-2)? true : false;

        // add spline to correct network
        m_curve_networks[layer]->add_polybezier(ctrl_points, curve_type, is_closed_curve);
    }

    m_current_layer = 0;

    // -----------------
    // update new canvas
    //m_canvas->new_drawing_canvas(image_width * width_factor, image_height * height_factor);
    m_canvas->new_drawing_canvas(image_width, image_height);
    m_curve_type = CURVE_CROSS;

    for (int l = 0; l < m_num_layers; ++l) {
        connect(m_curve_networks[l], SIGNAL(request_snapshot(QString, bool, bool, bool)), this, SLOT(slot_snapshot(QString, bool, bool, bool)));
    }
    return true;
} // TODO: display filename in window // TODO: include tool class in header // TODO: multi layer view

void Editor::slot_nlopt_solve()
{
    double time_tracker = -1.0;
    int nxhairs = -1;
    int nxcurves = -1;
    int nxpatches = -1;
    int npatches = -1;
    m_curve_networks[m_current_layer]->run_nlopt_solve(time_tracker, nxhairs, nxcurves, nxpatches, npatches);
    m_canvas->update();
}

void Editor::slot_nlopt_solve_all()
{
    QString directory_path = DEBUG_DIRECTORY + m_filename.toUpper() + "/";
    QDir dir(directory_path);
    if (!dir.exists())
        dir.mkpath(".");

    m_canvas->set_draw_flag_show_all_layers(true);
    m_canvas->set_draw_flag_show_normals(false);
    m_canvas->update();
    QString file_path = directory_path + "input.png";
    m_canvas->saveSnapshot(file_path, true);

    m_canvas->set_draw_flag_show_normals(true);
    m_canvas->update();
    file_path = directory_path + "input_normals.png";
    m_canvas->saveSnapshot(file_path, true);

    m_canvas->set_draw_flag_show_all_layers(false);

    for (int l = 0; l < m_num_layers; ++l) {
        set_active_layer(l);
        m_canvas->set_draw_flag_show_normals(false);
        m_canvas->update();
        file_path = directory_path + "nlopt_input_layer_" + QString::number(l) + ".png";
        m_canvas->saveSnapshot(file_path, true);

        m_canvas->set_draw_flag_show_normals(true);
        m_canvas->update();
        file_path = directory_path + "nlopt_input_normals_layer_" + QString::number(l) + ".png";
        m_canvas->saveSnapshot(file_path, true);
    }

    double time_tracker = 0.0;
    int nxhairs = 0;
    int nxcurves = 0;
    int nxpatches = 0;
    int npatches = 0;

    for (int l = 0; l < m_num_layers; ++l) {
        set_active_layer(l);
        m_curve_networks[l]->run_nlopt_solve(time_tracker, nxhairs, nxcurves, nxpatches, npatches);
        std::cout << time_tracker << std::endl;
        m_canvas->update();

        m_canvas->set_draw_flag_show_normals(false);
        m_canvas->set_draw_flag_show_fields(true);
        m_canvas->update();
        file_path = directory_path + "nlopt_solution_fields_layer_" + QString::number(l) + ".png";
        m_canvas->saveSnapshot(file_path, true);

        m_canvas->set_draw_flag_show_normals(true);
        m_canvas->set_draw_flag_show_fields(false);
        m_canvas->update();
        file_path = directory_path + "nlopt_solution_normals_layer_" + QString::number(l) + ".png";
        m_canvas->saveSnapshot(file_path, true);
    }

    m_canvas->set_draw_flag_show_all_layers(true);

    m_canvas->set_draw_flag_show_normals(true);
    m_canvas->set_draw_flag_show_fields(false);
    m_canvas->update();
    file_path = directory_path + "nlopt_solution_normals_all.png";
    m_canvas->saveSnapshot(file_path, true);

    m_canvas->set_draw_flag_show_normals(false);
    m_canvas->set_draw_flag_show_fields(true);
    m_canvas->update();
    file_path = directory_path + "nlopt_solution_fields_all.png";
    m_canvas->saveSnapshot(file_path, true);

    file_path = directory_path + "nlopt_total_solve_duration.txt";
    QFile file(file_path);
    if (file.open(QIODevice::ReadWrite)) {
        QTextStream stream(&file);
        stream << time_tracker << endl;
    }
    file.close();

    QString file_stats_path = directory_path + m_filename + "-curve_network_statistics.txt";
    QFile file_stats(file_stats_path);
    if (file_stats.open(QIODevice::ReadWrite)) {
        QTextStream stream(&file_stats);
        stream << "num xhairs: " << nxhairs << endl;
        stream << "num xcurves: " << nxcurves << endl;
        stream << "num xpatches: " << nxpatches << endl;
        stream << "num patches: " << npatches << endl;
    }
    file_stats.close();

    std::cout << m_filename.toStdString() << ": " << time_tracker << std::endl;
}

void Editor::slot_graph_cut()
{
    m_curve_networks[m_current_layer]->run_graph_cut();
    m_canvas->update();
}

void Editor::slot_process_mouse_event(mouse_action_t event, mouse_button_t button, point_t screen_point, point_t canvas_point)
{
    // TODO need actual flag for active canvas vs inactive
    if (m_current_layer == -1) {
        return;
    }

    // TODO this keeps repeating with all mouse events, fix it.
    switch (m_active_tool)
    {
    case TOOL_DRAW_FREE: { tool_draw_free(event, button, screen_point, canvas_point); break; }
    case TOOL_DRAW_POLYGON: { tool_draw_polygon(event, button, screen_point, canvas_point); break; }
    case TOOL_SELECT: { tool_select(event, button, screen_point, canvas_point); break; }
    case TOOL_INTERSECTION_NORMAL_FLIP: { tool_intersection_normal_flip(event, button, screen_point, canvas_point); break; }
    case TOOL_INTERSECTION_NORMAL_FLOP: { tool_intersection_normal_flop(event, button, screen_point, canvas_point); break; }
    case TOOL_REMOVE_INTERSECTION: { tool_remove_intersection(event, button, screen_point, canvas_point); break; }
    case TOOL_INSPECTOR: { tool_inspector(event, button, screen_point, canvas_point); break; }
    default: {
        break;
    }
    }
}

void Editor::tool_draw_free(mouse_action_t event, mouse_button_t button, point_t screen_point, point_t canvas_point)
{
    switch (event)
    {
    case MOUSE_PRESS: {
        m_canvas_raw_path.push_back(screen_point);
        m_canvas_draw_path.push_back(canvas_point);
        break;
    }
    case MOUSE_MOVE: {
        m_canvas_raw_path.push_back(screen_point);
        m_canvas_draw_path.push_back(canvas_point);
        break;
    }
    case MOUSE_RELEASE: {
        if (button == LEFT_BUTTON) {
            m_canvas_raw_path.push_back(screen_point);
            m_canvas_draw_path.push_back(canvas_point);

            QVector<point_t> raw_bezier_ctrl_points;
            QVector<point_t> bezier_ctrl_points;
            Bezier_Fit bezier_fit;
            bool is_closed_curve = bezier_fit.fitBezier(m_canvas_raw_path, raw_bezier_ctrl_points);
            for(int point = 0; point < raw_bezier_ctrl_points.size(); point++) {
                point_t bezier_ctrl_point = m_canvas->calculate_canvas_point(raw_bezier_ctrl_points[point]);
                bezier_ctrl_points.push_back(bezier_ctrl_point);
            }
            m_curve_networks[m_current_layer]->add_polybezier(bezier_ctrl_points, m_curve_type, is_closed_curve);
        }

        m_canvas_raw_path.clear();
        m_canvas_draw_path.clear();
        break;
    }
    default: {
        break;
    }
    }
}

void Editor::tool_draw_polygon(mouse_action_t event, mouse_button_t button, point_t screen_point, point_t canvas_point)
{
    switch (event)
    {
    case MOUSE_PRESS: {
        m_canvas_raw_path.push_back(screen_point);
        m_canvas_draw_path.push_back(canvas_point);
        break;
    }
    case MOUSE_MOVE: {
        point_t p0_raw = m_canvas_raw_path[0];
        point_t p0_draw = m_canvas_draw_path[0];
        m_canvas_raw_path.clear();
        m_canvas_draw_path.clear();

        vector_t v_raw = screen_point - p0_raw;
        vector_t v_draw = canvas_point - p0_draw;

        for (double p = 0.0; p <= 10.0; ++p) {
            point_t p_raw = p0_raw + v_raw * (p / 10.0);
            point_t p_draw = p0_draw + v_draw * (p / 10.0);
            m_canvas_raw_path.push_back(p_raw);
            m_canvas_draw_path.push_back(p_draw);
        }
        break;
    }
    case MOUSE_RELEASE: {
        if (button == LEFT_BUTTON) {
            m_canvas_raw_path.push_back(screen_point);
            m_canvas_draw_path.push_back(canvas_point);

            QVector<point_t> raw_bezier_ctrl_points;
            QVector<point_t> bezier_ctrl_points;
            Bezier_Fit bezier_fit;
            bool is_closed_curve = bezier_fit.fitBezier(m_canvas_raw_path, raw_bezier_ctrl_points);
            for(int point = 0; point < raw_bezier_ctrl_points.size(); point++) {
                point_t bezier_ctrl_point = m_canvas->calculate_canvas_point(raw_bezier_ctrl_points[point]);
                bezier_ctrl_points.push_back(bezier_ctrl_point);
            }
            m_curve_networks[m_current_layer]->add_polybezier(bezier_ctrl_points, m_curve_type, is_closed_curve);
        }

        m_canvas_raw_path.clear();
        m_canvas_draw_path.clear();
        break;
    }
    default: {
        break;
    }
    }
}

void Editor::tool_select(mouse_action_t event, mouse_button_t button, point_t screen_point, point_t canvas_point)
{

}

void Editor::tool_remove_intersection(mouse_action_t event, mouse_button_t button, point_t screen_point, point_t canvas_point)
{
    GLubyte pixel[4];
    m_canvas->draw_picking_pixel(screen_point, pixel);
    if (pixel[3] == 2) { // intersection
        boost::remove_edge(pixel[0], pixel[1], m_curve_networks[m_current_layer]->m_curve_network_graph);
    }
}

void Editor::tool_intersection_normal_flip(mouse_action_t event, mouse_button_t button, point_t screen_point, point_t canvas_point)
{
    if (event != MOUSE_RELEASE) {
        return;
    }

    GLubyte pixel[4];
    m_canvas->draw_picking_pixel(screen_point, pixel);
    if (pixel[3] == 2) { // intersection
        std::pair<curve_graph_edge_desc_t, bool> isec_desc_return = boost::edge(pixel[0], pixel[1], m_curve_networks[m_current_layer]->m_curve_network_graph);
        if (isec_desc_return.second == false) {
            std::cout << "[ ERROR ] selection of intersection that does not exist in database." << std::endl;
            return;
        }

        curve_graph_edge_data_map_t isec_ptr_map = boost::get(boost::edge_data_pointer, m_curve_networks[m_current_layer]->m_curve_network_graph);
        boost::shared_ptr<Intersection> isec = isec_ptr_map[isec_desc_return.first];
        isec->flip_tangent_frame();
    }
    m_canvas->update();
}

void Editor::tool_intersection_normal_flop(mouse_action_t event, mouse_button_t button, point_t screen_point, point_t canvas_point)
{
    if (event != MOUSE_RELEASE) {
        return;
    }

    GLubyte pixel[4];
    m_canvas->draw_picking_pixel(screen_point, pixel);
    if (pixel[3] == 2) { // intersection
        std::pair<curve_graph_edge_desc_t, bool> isec_desc_return = boost::edge(pixel[0], pixel[1], m_curve_networks[m_current_layer]->m_curve_network_graph);
        if (isec_desc_return.second == false) {
            std::cout << "[ ERROR ] selection of intersection that does not exist in database." << std::endl;
            return;
        }

        curve_graph_edge_data_map_t isec_ptr_map = boost::get(boost::edge_data_pointer, m_curve_networks[m_current_layer]->m_curve_network_graph);
        boost::shared_ptr<Intersection> isec = isec_ptr_map[isec_desc_return.first];
        isec->flip_tangent_frame();
    }
    m_canvas->update();
}

void Editor::tool_inspector(mouse_action_t event, mouse_button_t button, point_t screen_point, point_t canvas_point)
{
    m_current_selection.unselect();
    emit signal_deselect();

    GLubyte pixel[4];
    m_canvas->draw_picking_pixel(screen_point, pixel);

    switch (pixel[3])
    {
    case 1: { // polybez
        curve_graph_vertex_desc_t polybez_desc = boost::vertex(pixel[0], m_curve_networks[m_current_layer]->m_curve_network_graph);
        curve_graph_vertex_data_map_t polybez_ptr_map = boost::get(boost::vertex_data_pointer, m_curve_networks[m_current_layer]->m_curve_network_graph);

        polybez_ptr_t polybez_ptr = polybez_ptr_map[polybez_desc];
        m_current_selection.select_polybez(polybez_ptr->m_uuid);

        polybezier_info_form_s out_form = package_polybez_info(polybez_ptr);
        emit broadcast_polybezier_a_info(out_form);
        break;
    }
    case 2: { // intersection
        std::pair<curve_graph_edge_desc_t, bool> isec_desc_return = boost::edge(pixel[0], pixel[1], m_curve_networks[m_current_layer]->m_curve_network_graph);
        if (isec_desc_return.second == false) {
            std::cout << "[ ERROR ] selection of intersection that does not exist in database." << std::endl;
            return;
        }

        curve_graph_edge_data_map_t isec_ptr_map = boost::get(boost::edge_data_pointer, m_curve_networks[m_current_layer]->m_curve_network_graph);
        boost::shared_ptr<Intersection> isec_ptr = isec_ptr_map[isec_desc_return.first];
        m_current_selection.select_isec(isec_ptr->m_uuid);

        intersection_info_form_s out_form_isec = package_intersection_info(isec_ptr);
        emit broadcast_intersection_info(out_form_isec);

        curve_graph_vertex_data_map_t polybez_ptr_map = boost::get(boost::vertex_data_pointer, m_curve_networks[m_current_layer]->m_curve_network_graph);

        boost::optional<curve_graph_vertex_desc_t> polybez_a_desc = Curve_Network::find_vertex_by_uuid(isec_ptr->m_polybez_a_uuid, m_curve_networks[m_current_layer]->m_curve_network_graph);
        boost::optional<curve_graph_vertex_desc_t> polybez_b_desc = Curve_Network::find_vertex_by_uuid(isec_ptr->m_polybez_b_uuid, m_curve_networks[m_current_layer]->m_curve_network_graph);

        boost::shared_ptr<Polybezier> polybez_a = polybez_ptr_map[polybez_a_desc.value()];
        boost::shared_ptr<Polybezier> polybez_b = polybez_ptr_map[polybez_b_desc.value()];

        polybezier_info_form_s out_form_polybez_a = package_polybez_info(polybez_a);
        polybezier_info_form_s out_form_polybez_b = package_polybez_info(polybez_b);

        emit broadcast_polybezier_a_info(out_form_polybez_a);
        emit broadcast_polybezier_b_info(out_form_polybez_b);

        break;
    }
    default: {
        break;
    }
    }
}

void Editor::slot_snapshot(QString filepath, bool show_all, bool show_normals, bool show_fields)
{
    m_canvas->set_draw_flag_show_all_layers(show_all);
    m_canvas->set_draw_flag_show_normals(show_normals);
    m_canvas->set_draw_flag_show_fields(show_fields);
    m_canvas->update();

    if (show_all) { m_canvas->saveSnapshot(filepath + "-ALL-" + m_filename + ".png", true); }
    else { m_canvas->saveSnapshot(filepath + "-L" + QString::number(m_current_layer) + "-" + m_filename + ".png", true); }
}
