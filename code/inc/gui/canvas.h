#include <qglviewer.h>

#include "model/curve_network.h"

#include "controller/editor.h"

#include "utils/defs/util_gui_defs.h"
#include "utils/defs/util_data_defs.h"
#include "utils/math/util_math_defs.h"
#include "utils/gui/util_color_chart.h"

#include <QKeyEvent>
#include <QMouseEvent>
#include <QWidget>
#include <QDebug>

class Canvas : public QGLViewer
{
    Q_OBJECT

public slots:
    void slot_set_curve_thickness(int thickness);
    void slot_set_light_direction(float x, float y);

signals:
    // TODO this doesnt need to be a signal
    void request_process_mouse_event(mouse_action_t, mouse_button_t, point_t, point_t);
    void request_adjust_size(int width, int height);

public:
    Canvas(Editor *editor, QWidget *parent = 0);
    ~Canvas();

    void new_drawing_canvas(int image_width, int image_height);

    // ----------------------
    // CANVAS PICKING
    // ----------------------
    void draw_picking_pixel(point_t mouse_point, GLubyte* pixel);

    // ----------------------
    // CANVAS HELPERS
    // ----------------------
    point_t calculate_canvas_point(point_t point);

    // ----------------------
    // CANVAS FLAG SETTERS
    // ----------------------
    void set_canvas_display_mode(display_mode_t mode) { m_canvas_display_mode = mode; }

    void set_draw_flag_show_curves(bool flag) { m_draw_flag_show_curves = flag; }
    void set_draw_flag_show_isecs(bool flag) { m_draw_flag_show_isecs = flag; }
    void set_draw_flag_show_planes(bool flag) { m_draw_flag_show_planes = flag; }
    void set_draw_flag_show_normals(bool flag) { m_draw_flag_show_normals = flag; }
    void set_draw_flag_show_fields(bool flag) { m_draw_flag_show_fields = flag; }
    void set_draw_flag_show_patches(bool flag) { m_draw_flag_show_patches = flag; }
    void set_draw_flag_show_labels(bool flag) { m_draw_flag_show_labels = flag; }
    void set_draw_flag_show_all_layers(bool flag) { m_draw_flag_show_all_layers = flag; }

    // TODO: move to private
    // ----------------------
    // CANVAS PROPERTIES
    // ----------------------
    double m_canvas_zoom;
    point_t m_canvas_center;
    
protected:
    virtual void draw();
    virtual QString helpString() const;
    virtual void init(); // deprecated
    virtual void drawLight(GLenum light,float scale = 1.0f) const { } // necessary to avoid compiling errors

private:
    virtual void keyPressEvent(QKeyEvent *e);
    virtual void mousePressEvent(QMouseEvent* e);
    virtual void mouseMoveEvent(QMouseEvent *e);
    virtual void mouseReleaseEvent(QMouseEvent *e);

    void draw_2D();
    void draw_3D();

    void draw_2D_freehand();
    void draw_2D_curves(int layer, float opacity);
    void draw_2D_intersections(int layer, float opacity);
    void draw_2D_fields(int layer, float opacity);
    void draw_2D_tangent_frames(int layer, float opacity);
    void draw_2D_planes(int layer, float opacity);
    void draw_2D_patches(int layer, float opacity);
    void draw_2D_labels(int layer, float opacity);
    void draw_2D_debug(int layer);

    // ----------------------
    // CANVAS STATES
    // ----------------------
    display_mode_t m_canvas_display_mode;
    bool m_canvas_is_panning;


    // ----------------------
    // CANVAS DRAWING FLAGS
    // ----------------------
    bool m_draw_flag_show_curves;
    bool m_draw_flag_show_isecs;
    bool m_draw_flag_show_planes;
    bool m_draw_flag_show_normals;
    bool m_draw_flag_show_fields;
    bool m_draw_flag_show_patches;
    bool m_draw_flag_show_labels;
    bool m_draw_flag_show_all_layers;


    // -----------------------
    // CANVAS SP. RENDER DATA
    // -----------------------
    int m_draw_curve_thickness;
    int m_draw_point_size;
    float m_light_dir_x;
    float m_light_dir_y;
    float m_light_dir_z;


    // ----------------------
    // CANVAS PICKING BUFFERS
    // ----------------------
    GLuint m_buffer_picking;
    GLuint m_buffer_picking_fbo;


    // ----------------------
    // TOOL CONTROLLER
    // ----------------------
    Editor* m_editor;
};
