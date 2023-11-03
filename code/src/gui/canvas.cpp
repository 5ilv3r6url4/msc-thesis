#include <GL/glew.h>

#include "gui/canvas.h"

#include <stdio.h>

#include <iostream>
#include <QFileDialog>
#include <QTextStream>

#define TWO_M_PI	(2.*M_PI)

/*------------------------------------------ CANVAS ---
 |
 |  Canvas is where all the drawing is done and user
 |  interactions with the drawing area is triggered.
 |  Just creating a canvas is insufficient to begin
 |  drawing and interacting, a new file needs to be
 |  created in the editor.
 |
 *-----------------------------------------------------*/
Canvas::Canvas(Editor* editor, QWidget *parent)
    : QGLViewer(parent)
{
    // --------------------------
    // link the editor and canvas
    editor->register_canvas(this);
    m_editor = editor;

    // ---------------
    // canvas defaults
    // properties
    m_canvas_zoom               = 1.0;
    m_canvas_center             = point_t(0.0, 0.0);

    // state
    m_canvas_display_mode       = DISPLAY_2D;
    m_canvas_is_panning         = false;

    // draw flags
    m_draw_flag_show_curves     = true;
    m_draw_flag_show_isecs      = true;
    m_draw_flag_show_planes     = false;
    m_draw_flag_show_normals    = false;
    m_draw_flag_show_fields     = false;
    m_draw_flag_show_patches    = false;
    m_draw_flag_show_labels     = false;
    m_draw_flag_show_all_layers = false;

    // sp. render data
    m_draw_curve_thickness      = 1;
    m_draw_point_size           = 5;
    m_light_dir_x               = 0.5f;
    m_light_dir_y               = 1.0f;
    m_light_dir_z               = 1.0f;
}

Canvas::~Canvas()
{
    glDeleteTextures(1,&m_buffer_picking);
    glDeleteFramebuffersEXT(1,&m_buffer_picking_fbo);
}

/*-------------------------------------------- INIT ---*
 |
 |  QGLViewer virtual.
 |
 *-----------------------------------------------------*/
void Canvas::init() {
    // ---------
    // init glew
    glewInit();
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);

    // ---------------------------------------------
    // generate size 0 empty buffers for safe delete
    glGenTextures(1, &m_buffer_picking);
    glBindTexture( GL_TEXTURE_RECTANGLE_ARB, m_buffer_picking);
    glTexImage2D( GL_TEXTURE_RECTANGLE_ARB, 0, GL_RGBA32F_ARB, 0, 0, 0, GL_RGBA, GL_FLOAT, NULL);
    glTexParameterf(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
    glTexParameterf(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);

    glGenFramebuffersEXT(1, &m_buffer_picking_fbo);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, m_buffer_picking_fbo);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_RECTANGLE_ARB, m_buffer_picking, 0);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}

/*-------------------------------------------- DRAW ---
 |
 |  Central point for all draw calls.
 |
 *-----------------------------------------------------*/
void Canvas::draw()
{
    // TODO need actual flag for active canvas vs inactive
    if (m_editor->m_current_layer == -1) {
        glClearColor(charcoal[0], charcoal[1], charcoal[2], charcoal[3]);
        glClear(GL_COLOR_BUFFER_BIT);
        glViewport(0, 0, width(), height());
    }

    else {
        glClearColor(white[0], white[1], white[2], white[3]);
        glClear(GL_COLOR_BUFFER_BIT);
        glViewport(0, 0, width(), height());
        if (!m_editor->m_curve_networks.isEmpty()) {
            switch(m_canvas_display_mode)
            {
            case DISPLAY_2D:
                draw_2D();
                break;
            case DISPLAY_3D:
                draw_3D();
                break;
            default:
                draw_2D();
                break;
            }
        }
    }
}
// TODO disable all buttons when no file open in canvas

/*--------------------------------------- DRAW_3D -----
 |
 |  Central point for all 3D draw calls. TODO
 |
 *-----------------------------------------------------*/
void Canvas::draw_3D() { }

/*--------------------------------------- DRAW_2D -----
 |
 |  Central point for all 2D draw calls.
 |
 *-----------------------------------------------------*/
void Canvas::draw_2D()
{
    glMatrixMode(GL_PROJECTION); glLoadIdentity();
    glMatrixMode(GL_MODELVIEW) ; glLoadIdentity();

    // -------------------
    // enable antialiasing
    // -------------------
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_POLYGON_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    int curr_layer = m_editor->m_current_layer;

    // ------------------------------
    // draw curve being drawn by user
    // ------------------------------
    if ( m_editor->m_canvas_draw_path.size() > 0 ) { draw_2D_freehand(); }

    for (int l = 0; l < m_editor->m_num_layers; ++l) {
        // ------------------------------
        // opacity toggle
        // ------------------------------
        float opacity = 0.25f;
        if (l == curr_layer) { opacity = 1.0f; }

        // ------------------------------
        // draw polybez curves in network
        // ------------------------------
        if ( m_draw_flag_show_curves ) { draw_2D_curves(l, opacity); }

        // ------------------------------
        // opacity toggle
        // ------------------------------
        if (m_draw_flag_show_all_layers || l == curr_layer) { opacity = 1.0f; }
        else { opacity = 0.0f; }

        // -----------------------------
        // draw intersections in network
        // -----------------------------
        if ( m_draw_flag_show_isecs ) { draw_2D_intersections(l, opacity); }

        // ------------------------------------
        // draw polybez curve planes in network
        // ------------------------------------
        if ( m_draw_flag_show_planes ) { draw_2D_planes(l, opacity); }

        // ------------------------------------
        // draw intersection normals in network
        // ------------------------------------
        if ( m_draw_flag_show_normals ) { draw_2D_tangent_frames(l, opacity); }

        // ------------------------------------------------
        // draw normal fields along cross curves in network
        // ------------------------------------------------
        if ( m_draw_flag_show_fields ) { draw_2D_fields(l, opacity); }

        // ------------------------------------
        // draw disjoint patches in network
        // ------------------------------------
        if ( m_draw_flag_show_patches ) { draw_2D_patches(l, opacity); }
    }

    // TODO doubles
    glMatrixMode(GL_MODELVIEW); glLoadIdentity();
    glColor3f(1.0,1.0,1.0);

    //disable antialiasing
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_POLYGON_SMOOTH);
    glDisable(GL_BLEND);
}

void Canvas::draw_2D_freehand()
{
    glMatrixMode(GL_MODELVIEW);
    glScalef(m_canvas_zoom, m_canvas_zoom, 0.f);
    glTranslatef(m_canvas_center.x, m_canvas_center.y, 0.0f);

    glColor3fv(gray);
    glLineWidth(m_draw_curve_thickness);
    glBegin(GL_LINE_STRIP);

    for(int p = 0; p < m_editor->m_canvas_draw_path.size(); p++){
        double x = m_editor->m_canvas_draw_path[p].x;
        double y = m_editor->m_canvas_draw_path[p].y;
        glVertex3f(x, y, 0.f);
    }
    glEnd();
}

void Canvas::draw_2D_curves(int layer, float opacity)
{
    const Curve_Network* curve_network_readonly = m_editor->get_network_model(layer);
    curve_graph_t curve_network_graph = curve_network_readonly->m_curve_network_graph;
    curve_graph_vertex_data_map_t polybez_ptr_map = boost::get(boost::vertex_data_pointer, curve_network_graph);

    if (m_draw_flag_show_all_layers) { opacity = 1.0f; }

    glMatrixMode(GL_MODELVIEW);
    glScalef(m_canvas_zoom, m_canvas_zoom, 0.f);
    glTranslatef(m_canvas_center.x, m_canvas_center.y, 0.0f);

    glLineWidth(m_draw_curve_thickness);
    for (auto polybez_desc : boost::make_iterator_range(boost::vertices(curve_network_graph))) {

        // colour coded curve types
        switch (polybez_ptr_map[polybez_desc]->m_curve_type)
        {
        case CURVE_CROSS: { glColor4f(orange[0], orange[1], orange[2], opacity); break; }
        case CURVE_BOUNDARY: { glColor4f(green[0], green[1], green[2], opacity); break; }
        case CURVE_SILHOUETTE: { glColor4f(blue[0], blue[1], blue[2], opacity);  break; }
        default: { glColor4fv(black); break; }
        }

        for (int bz = 0; bz < polybez_ptr_map[polybez_desc]->m_bezier_curves.size(); ++bz) {
            glBegin(GL_LINE_STRIP);
            for(double t = 0.f; t <= 1.f; t += T_STEP){
                point_t P = polybez_ptr_map[polybez_desc]->compute_point_at(bz, t);
                glVertex3f(P.x, P.y, 0.f);
            }
            glEnd();
        }
    }
}

void Canvas::draw_2D_intersections(int layer, float opacity)
{
    const Curve_Network* curve_network_readonly = m_editor->get_network_model(layer);
    curve_graph_t curve_network_graph = curve_network_readonly->m_curve_network_graph;
    curve_graph_edge_data_map_t isec_ptr_map = boost::get(boost::edge_data_pointer, curve_network_graph);

    glMatrixMode(GL_MODELVIEW);
    glScalef(m_canvas_zoom, m_canvas_zoom, 0.f);
    glTranslatef(m_canvas_center.x, m_canvas_center.y, 0.0f);

    glColor4f(blue[0], blue[1], blue[2], opacity);
    glPointSize(m_draw_point_size);
    glBegin(GL_POINTS);

    for (auto isec_desc : boost::make_iterator_range(boost::edges(curve_network_graph))) {
        point_t isec_pos = isec_ptr_map[isec_desc]->m_position;
        glVertex3f(isec_pos.x, isec_pos.y, 0.0f);
    }
    glEnd();
}

void Canvas::draw_2D_tangent_frames(int layer, float opacity)
{
    const Curve_Network* curve_network_readonly = m_editor->get_network_model(layer);
    curve_graph_t curve_network_graph = curve_network_readonly->m_curve_network_graph;
    curve_graph_vertex_data_map_t polybez_ptr_map = boost::get(boost::vertex_data_pointer, curve_network_graph);
    curve_graph_edge_data_map_t isec_ptr_map = boost::get(boost::edge_data_pointer, curve_network_graph);

    glMatrixMode(GL_MODELVIEW);
    glScalef(m_canvas_zoom, m_canvas_zoom, 0.f);
    glTranslatef(m_canvas_center.x, m_canvas_center.y, 0.0f);

    glLineWidth(m_draw_curve_thickness);

    for (auto isec : boost::make_iterator_range(boost::edges(curve_network_graph))) {
        if (polybez_ptr_map[isec.m_target]->m_curve_type != CURVE_CROSS || polybez_ptr_map[isec.m_source]->m_curve_type != CURVE_CROSS)
            continue;

        point_t  isec_pos = isec_ptr_map[isec]->m_position;
        vector_t isec_tan_a = isec_ptr_map[isec]->m_tangent_a * 0.1;
        vector_t isec_tan_b = isec_ptr_map[isec]->m_tangent_b * 0.1;
        vector_t isec_normal = isec_ptr_map[isec]->m_normal * 0.2;

        if (isec_tan_a.z > 0.0) { glColor4f(red[0], red[1], red[2], opacity); }
        else { glColor4f(blush[0], blush[1], blush[2], opacity); }
        glBegin(GL_LINE_STRIP);
        glVertex3f(isec_pos.x, isec_pos.y, 0.f);
        glVertex3f(isec_pos.x + isec_tan_a.x, isec_pos.y + isec_tan_a.y, 0.f);
        glEnd();

        if (isec_tan_b.z > 0.0) { glColor4f(red[0], red[1], red[2], opacity); }
        else { glColor4f(blush[0], blush[1], blush[2], opacity); }
        glBegin(GL_LINE_STRIP);
        glVertex3f(isec_pos.x, isec_pos.y, 0.f);
        glVertex3f(isec_pos.x + isec_tan_b.x, isec_pos.y + isec_tan_b.y, 0.f);
        glEnd();

        glColor4f(maroon[0], maroon[1], maroon[2], opacity);
        glBegin(GL_LINE_STRIP);
        glVertex3f(isec_pos.x, isec_pos.y, 0.f);
        glVertex3f(isec_pos.x + isec_normal.x, isec_pos.y + isec_normal.y, 0.f);
        glEnd();
    }
}

void Canvas::draw_2D_fields(int layer, float opacity)
{
    const Curve_Network* curve_network_readonly = m_editor->get_network_model(layer);
    curve_graph_t curve_network_graph = curve_network_readonly->m_curve_network_graph;
    curve_graph_vertex_data_map_t polybez_ptr_map = boost::get(boost::vertex_data_pointer, curve_network_graph);

    glMatrixMode(GL_MODELVIEW);
    glScalef(m_canvas_zoom, m_canvas_zoom, 0.f);
    glTranslatef(m_canvas_center.x, m_canvas_center.y, 0.0f);

    glLineWidth(m_draw_curve_thickness);

    for (auto polybez_desc : boost::make_iterator_range(boost::vertices(curve_network_graph))) {
        if (polybez_ptr_map[polybez_desc]->m_curve_type != CURVE_CROSS)
            continue;

        for (int i = 0; i < polybez_ptr_map[polybez_desc]->m_integrated_normals.size(); ++i) {
            point_t p = polybez_ptr_map[polybez_desc]->m_integrated_normals[i].position;
            vector_t n = (polybez_ptr_map[polybez_desc]->m_integrated_normals[i].normal) * 0.1;
            vector_t nc = {n.x, n.y, n.z};
            nc.normalize();

            glBegin(GL_LINE_STRIP);
            glColor4f((float)nc.x, (float)nc.y, (float)nc.z, opacity);
            glVertex3f(p.x, p.y, 0.f);
            glVertex3f(p.x + n.x, p.y + n.y, 0.f);
            glEnd();
        }
    }
}

void Canvas::draw_2D_patches(int layer, float opacity)
{
    const Curve_Network* curve_network_readonly = m_editor->get_network_model(layer);
    QVector<curve_graph_t> curve_network_subgraphs = curve_network_readonly->m_curve_network_subgraphs;

    int offset = 0;
    if (m_draw_flag_show_all_layers) {
        opacity = 1.0f;
        if (layer != 0) {
            for (int l = 0; l < layer; ++l) {
                const Curve_Network* curve_network_readonly_prev = m_editor->get_network_model(l);
                QVector<curve_graph_t> curve_network_subgraphs_prev = curve_network_readonly_prev->m_curve_network_subgraphs;
                offset += curve_network_subgraphs_prev.size();
            }
        }
    }

    glMatrixMode(GL_MODELVIEW);
    glScalef(m_canvas_zoom, m_canvas_zoom, 0.f);
    glTranslatef(m_canvas_center.x, m_canvas_center.y, 0.0f);

    glLineWidth(m_draw_curve_thickness);

    for (int i = 0; i < curve_network_subgraphs.size(); ++i) {
        curve_graph_t curve_network_subgraph = curve_network_subgraphs[i];
        curve_graph_vertex_data_map_t polybez_ptr_map = boost::get(boost::vertex_data_pointer, curve_network_subgraph);
        for (auto polybez_desc : boost::make_iterator_range(boost::vertices(curve_network_subgraph))) {
            // colour coded curve types
            switch (i + offset)
            {
            case 0:  { glColor4f(purple[0], purple[1], purple[2], opacity);         break; }
            case 1:  { glColor4f(lime[0], lime[1], lime[2], opacity);               break; }
            case 2:  { glColor4f(navy[0], navy[1], navy[2], opacity);               break; }
            case 3:  { glColor4f(cyan[0], cyan[1], cyan[2], opacity);               break; }
            case 4:  { glColor4f(orange[0], orange[1], orange[2], opacity);         break; }
            case 5:  { glColor4f(pink[0], pink[1], pink[2], opacity);               break; }
            case 6:  { glColor4f(lavendar[0], lavendar[1], lavendar[2], opacity);   break; }
            case 7:  { glColor4f(maroon[0], maroon[1], maroon[2], opacity);         break; }
            case 8:  { glColor4f(magenta[0], magenta[1], magenta[2], opacity);      break; }
            case 9:  { glColor4f(yellow[0], yellow[1], yellow[2], opacity);         break; }
            case 10: { glColor4f(teal[0], teal[1], teal[2], opacity);               break; }
            case 11: { glColor4f(olive[0], olive[1], olive[2], opacity);            break; }
            default: { glColor4f(black[0], black[1], black[2], opacity);            break; }
            }

            for (int bz = 0; bz < polybez_ptr_map[polybez_desc]->m_bezier_curves.size(); ++bz) {
                glBegin(GL_LINE_STRIP);
                for(double t = 0.f; t <= 1.f; t += T_STEP){
                    point_t P = polybez_ptr_map[polybez_desc]->compute_point_at(bz, t);
                    glVertex3f(P.x, P.y, 0.f);
                }
                glEnd();
            }
        }
    }
}

void Canvas::draw_2D_labels(int layer, float opacity)
{
    glMatrixMode(GL_MODELVIEW);
    glScalef(m_canvas_zoom, m_canvas_zoom, 0.f);
    glTranslatef(m_canvas_center.x, m_canvas_center.y, 0.0f);
    glLineWidth(m_draw_curve_thickness);

    const Curve_Network* curve_network_readonly = m_editor->get_network_model(layer);
    curve_graph_t curve_network_graph = curve_network_readonly->m_curve_network_graph;
    curve_graph_vertex_data_map_t polybez_ptr_map = boost::get(boost::vertex_data_pointer, curve_network_graph);

    for (auto polybez : boost::make_iterator_range(boost::vertices(curve_network_graph))) {
        QFont serifFont("Times", 12, QFont::Bold);
        QString label = QString::number(polybez);
        double xview = polybez_ptr_map[polybez]->m_bezier_curves[0].ctrl_points[0].x;
        double yview = polybez_ptr_map[polybez]->m_bezier_curves[0].ctrl_points[0].y;
        int xscreen = (int)((xview + 1.0)*(0.5 * width()));
        int yscreen = (int)((yview - 1.0)*(-0.5 * height()));
        glColor4f(black[0], black[1], black[2], opacity);
        drawText(xscreen, yscreen, label, serifFont);
    }
}

void Canvas::draw_2D_planes(int layer, float opacity)
{
    glMatrixMode(GL_MODELVIEW);
    glScalef(m_canvas_zoom, m_canvas_zoom, 0.f);
    glTranslatef(m_canvas_center.x, m_canvas_center.y, 0.0f);
    glLineWidth(m_draw_curve_thickness);

    // TODO
}

void Canvas::draw_2D_debug(int layer)
{
    glMatrixMode(GL_MODELVIEW);
    glScalef(m_canvas_zoom, m_canvas_zoom, 0.f);
    glTranslatef(m_canvas_center.x, m_canvas_center.y, 0.0f);
    glLineWidth(m_draw_curve_thickness);

    const Curve_Network* curve_network_readonly = m_editor->get_network_model(layer);
    curve_graph_t curve_network_graph = curve_network_readonly->m_curve_network_graph;
    QVector<curve_graph_t> curve_network_subgraphs = curve_network_readonly->m_curve_network_subgraphs;

    curve_graph_vertex_data_map_t polybez_ptr_map = boost::get(boost::vertex_data_pointer, curve_network_graph);

    // ADD AS NEEDED
}

/*----------------------------------------- DRAW_WITH_ID -----
 |
 |  Draws curves, isecs, and control points using
 |  RGBA to object index along RGB and object type along A.
 |
 *------------------------------------------------------------*/
void Canvas::draw_picking_pixel(point_t mouse_point, GLubyte *pixel)
{
    // TODO const curve_network_graph -- readonly
    const Curve_Network* curve_network_readonly = m_editor->get_network_model(m_editor->m_current_layer);
    curve_graph_t curve_network_graph = curve_network_readonly->m_curve_network_graph;

    curve_graph_vertex_data_map_t polybez_ptr_map = boost::get(boost::vertex_data_pointer, curve_network_graph);
    curve_graph_edge_data_map_t isec_ptr_map = boost::get(boost::edge_data_pointer, curve_network_graph);

    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, m_buffer_picking_fbo);

    glDisable(GL_BLEND);
    glDisable(GL_DITHER);
    glDisable(GL_FOG);
    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_1D);
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_TEXTURE_3D);
    glShadeModel(GL_FLAT);

    glMatrixMode(GL_PROJECTION); glLoadIdentity();
    glMatrixMode(GL_MODELVIEW) ; glLoadIdentity();

    glClearColor(0.f, 0.f, 0.f, 0.f);
    glClear(GL_COLOR_BUFFER_BIT);
    glViewport(0, 0, width(), height());

    glScalef(m_canvas_zoom, m_canvas_zoom, 0.f);
    glTranslatef(m_canvas_center.x, m_canvas_center.y, 0.0f);

    // ----------------------------------------------------------------
    // draw curves with their boost graph vertex index in the R channel
    // ----------------------------------------------------------------
    glLineWidth(m_draw_curve_thickness * 10); // TODO make picking sizes constant/static/define? - in global params?
    for (auto polybez_desc : boost::make_iterator_range(boost::vertices(curve_network_graph))) {
        glColor4ub(polybez_desc, 0, 0, 1);
        for (int bz = 0; bz < polybez_ptr_map[polybez_desc]->m_bezier_curves.size(); ++bz) {
            glBegin(GL_LINE_STRIP);
            for(double t = 0.f; t <= 1.f; t += T_STEP){
                point_t P = polybez_ptr_map[polybez_desc]->compute_point_at(bz, t);
                glVertex3f(P.x, P.y, 0.f);
            }
            glEnd();
        }
    }

    // ----------------------------------------------------------------------------
    // draw intersections with their boost graph edge index in the R and G channels
    // ----------------------------------------------------------------------------
    glPointSize(m_draw_point_size * 3);
    glBegin(GL_POINTS);
    for (auto isec_desc : boost::make_iterator_range(boost::edges(curve_network_graph))) {
        curve_graph_vertex_desc_t polybez_a_desc = boost::source(isec_desc, curve_network_graph);
        curve_graph_vertex_desc_t polybez_b_desc = boost::target(isec_desc, curve_network_graph);

        curve_graph_vertex_data_map_t polybez_ptr_map = boost::get(boost::vertex_data_pointer, curve_network_graph);
        curve_graph_edge_data_map_t isec_ptr_map = boost::get(boost::edge_data_pointer, curve_network_graph);

        if (isec_ptr_map[isec_desc]->m_polybez_a_uuid == polybez_ptr_map[polybez_a_desc]->m_uuid)
            //    std::cout << "yolo" << std::endl;
        {true; }

        glColor4ub(polybez_a_desc, polybez_b_desc, 0, 2);
        point_t isec_pos = isec_ptr_map[isec_desc]->m_position;
        glVertex3f(isec_pos.x, isec_pos.y, 0.0f);
    }
    glEnd();

    glReadPixels(mouse_point.x, height() - mouse_point.y - 1, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, (void *)pixel);

    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}

/*------------------------------------ KEYPRESSEVENT ---
 |
 |  Passes state of canvas and key to editor.
 |
 *-----------------------------------------------------*/
void Canvas::keyPressEvent(QKeyEvent *e) { }

/*---------------------------------- MOUSEPRESSEVENT ---
 |
 |  Passes mouse info in screen and in canvas to editor.
 |
 *-----------------------------------------------------*/
void Canvas::mousePressEvent(QMouseEvent* e)
{
    switch(m_canvas_display_mode)
    {
    case DISPLAY_2D: {
        point_t raw_point = point_t(e->x(), e->y());
        point_t draw_point = calculate_canvas_point(raw_point);
        mouse_button_t button = mouse_button_t(e->button());

        emit request_process_mouse_event(MOUSE_PRESS, button, raw_point, draw_point);
        break;
    }
    default: { QGLViewer::mousePressEvent(e); break; }
    }
    update();
}

/*----------------------------------- MOUSEMOVEEVENT ---
 |
 |  Passes mouse info in screen and in canvas to editor.
 |
 *-----------------------------------------------------*/
void Canvas::mouseMoveEvent(QMouseEvent* e)
{
    switch(m_canvas_display_mode)
    {
    case DISPLAY_2D: {
        point_t raw_point = point_t(e->x(), e->y());
        point_t draw_point = calculate_canvas_point(raw_point);
        mouse_button_t button = mouse_button_t(e->button());

        emit request_process_mouse_event(MOUSE_MOVE, button, raw_point, draw_point);
        break;
    }
    default: { QGLViewer::mousePressEvent(e); break; }
    }
    update();
}

/*-------------------------------- MOUSERELEASEEVENT ---
 |
 |  Passes mouse info in screen and in canvas to editor.
 |
 *-----------------------------------------------------*/
void Canvas::mouseReleaseEvent(QMouseEvent* e)
{
    switch(m_canvas_display_mode)
    {
    case DISPLAY_2D: {
        point_t raw_point = point_t(e->x(), e->y());
        point_t draw_point = calculate_canvas_point(raw_point);
        mouse_button_t button = mouse_button_t(e->button());

        emit request_process_mouse_event(MOUSE_RELEASE, button, raw_point, draw_point);
        break;
    }
    default: { QGLViewer::mousePressEvent(e); break; }
    }
    update();
}

/*--------------------------------------- HELPSTRING ---
 |
 |  QGLViewer help built-in.
 |
 *-----------------------------------------------------*/
QString Canvas::helpString() const
{
    QString text("<h2>V i e w e r</h2>");
    text += "No help yet... ";
    return text;
}

/*--------------------------- CALCULATE_CANVAS_POINT ---
 |
 |  Returns canvas coordinates from screen coordinates.
 |
 *-----------------------------------------------------*/
point_t Canvas::calculate_canvas_point(point_t point)
{
    double x = point.x / (0.5 * width()) - 1.0;
    double y = 1.0 - point.y / (0.5 * height());
    x /= m_canvas_zoom;
    y /= m_canvas_zoom;
    x -= m_canvas_center.x;
    y -= m_canvas_center.y;
    return point_t(x, y);
}

/*------------------------- SLOT_SET_CURVE_THICKNESS ---
 |
 |  Slot that updates the curve draw thickness.
 |
 *-----------------------------------------------------*/
void Canvas::slot_set_curve_thickness(int thickness)
{
    m_draw_curve_thickness = thickness;
    update();
}

/*------------------------- SLOT_SET_LIGHT_DIRECTION ---
 |
 |  Slot that sets the light direction for sp rendering.
 |
 *-----------------------------------------------------*/
void Canvas::slot_set_light_direction(float x, float y)
{
    //convert from image space to spherical space
    float phi = (0.5-0.5*x)*TWO_M_PI;
    float theta = y*M_PI;
    //convert from spherical space to xyz space
    float sintheta = sinf(theta);
    m_light_dir_x = sintheta * cosf(phi);
    m_light_dir_z = sintheta * sinf(phi);
    m_light_dir_y = cosf(theta);

    //float lightX = x*2.-1.;
    //float lightY = -(y*2.-1.);

    //diffShader.updateLightDir(lx,ly,lz);
    update();
}

/*---------------------------------------- NEW_FILE ---
 |
 |  New picking buffer, update dimensions, reset flags.
 |
 *-----------------------------------------------------*/
void Canvas::new_drawing_canvas(int image_width, int image_height)
{
    // --------------
    // delete buffers
    glDeleteTextures(1, &m_buffer_picking);
    glDeleteFramebuffersEXT(1, &m_buffer_picking_fbo);

    // --------------------
    // generate new buffers
    glGenTextures(1, &m_buffer_picking);
    glBindTexture( GL_TEXTURE_RECTANGLE_ARB, m_buffer_picking);
    glTexImage2D( GL_TEXTURE_RECTANGLE_ARB, 0, GL_RGBA32F_ARB, image_width, image_height, 0, GL_RGBA, GL_FLOAT, NULL);
    glTexParameterf(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
    glTexParameterf(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);

    glGenFramebuffersEXT(1, &m_buffer_picking_fbo);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, m_buffer_picking_fbo);
    glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_RECTANGLE_ARB, m_buffer_picking, 0);
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);

    // -----------------------
    // resize qglviewer widget
    resize(image_width, image_height);
    int center_x = (GRID_WIDTH / 2) - (image_width / 2);
    int center_y = (GRID_HEIGHT / 2) - (image_height / 2);
    this->move(center_x, center_y);

    // ----------------
    // restore defaults
    // properties
    m_canvas_zoom = 1.0;
    m_canvas_center = point_t(0.0, 0.0);

    // state
    m_canvas_display_mode       = DISPLAY_2D;
    m_canvas_is_panning         = false;

    // draw flags
    m_draw_flag_show_curves     = true;
    m_draw_flag_show_isecs      = true;
    m_draw_flag_show_planes     = false;
    m_draw_flag_show_normals    = false;
    m_draw_flag_show_fields     = false;
    m_draw_flag_show_patches    = false;
    m_draw_flag_show_labels     = false;
    m_draw_flag_show_all_layers = false;

    // sp. render data
    m_draw_curve_thickness  = 1;
    m_draw_point_size       = 5;
    m_light_dir_x           = 0.5f;
    m_light_dir_y           = 1.0f;
    m_light_dir_z           = 1.0f;

    update();
}
