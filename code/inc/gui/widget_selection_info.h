#ifndef WIDGET_SELECTION_INFO_H
#define WIDGET_SELECTION_INFO_H

#include "utils/gui/util_gui_forms.h"

// ---------------------------------
// QT incompatability with boost fix
// ---------------------------------
#ifndef Q_MOC_RUN
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/optional.hpp>
#endif

#include <QWidget>
#include <QComboBox>
#include <QDoubleValidator>

namespace Ui {
class Widget_Selection_Info;
}

class Widget_Selection_Info : public QWidget
{
    Q_OBJECT

public slots:    
    void write_polybez_a_form_info(polybezier_info_form_s& form);
    void write_polybez_b_form_info(polybezier_info_form_s& form);
    void write_intersection_form_info(intersection_info_form_s& form);

    void clear_all_forms_info();

signals:
    void signal_spnbox_polybez_a_segment(int);
    void signal_spnbox_polybez_b_segment(int);

    void request_selected_polybez_curve_type_update(curve_type_t);
    void request_selected_polybez_plane_normal_update(vector_t);
    void request_selected_polybez_plane_offset_update(double);
    void request_selected_polybez_ctrl_point_update(int, int, point_t);

    //void request_selected_intersection_position_update(point_t);
    //void request_selected_intersection_tangents_update(boost::optional<vector_t>, boost::optional<vector_t>);

public:
    explicit Widget_Selection_Info(QWidget *parent = nullptr);
    ~Widget_Selection_Info();

private:
    Ui::Widget_Selection_Info *ui;
    QDoubleValidator* m_input_validator_double;

    void send_polybez_form_info(int flag);
    void send_intersection_form_info() {}

    void clear_polybez_a_form_info();
    void clear_polybez_b_form_info();
    void clear_intersection_form_info();
};

#endif // WIDGET_SELECTION_INFO_H
