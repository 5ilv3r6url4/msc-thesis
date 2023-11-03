#include "ui_widget_selection_info.h"

#include "gui/widget_selection_info.h"

Widget_Selection_Info::Widget_Selection_Info(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Widget_Selection_Info)
{
    ui->setupUi(this);

    m_input_validator_double = new QDoubleValidator(-1, 1, 3);

    // ------------------------------------------------------------------------
    // polybezier a input formattings
    ui->input_polybez_a_normal_x->setValidator(m_input_validator_double);
    ui->input_polybez_a_normal_y->setValidator(m_input_validator_double);
    ui->input_polybez_a_normal_z->setValidator(m_input_validator_double);
    ui->input_polybez_a_ctrl_point_0_x->setValidator(m_input_validator_double);
    ui->input_polybez_a_ctrl_point_0_y->setValidator(m_input_validator_double);
    ui->input_polybez_a_ctrl_point_0_z->setValidator(m_input_validator_double);
    ui->input_polybez_a_ctrl_point_1_x->setValidator(m_input_validator_double);
    ui->input_polybez_a_ctrl_point_1_y->setValidator(m_input_validator_double);
    ui->input_polybez_a_ctrl_point_1_z->setValidator(m_input_validator_double);
    ui->input_polybez_a_ctrl_point_2_x->setValidator(m_input_validator_double);
    ui->input_polybez_a_ctrl_point_2_y->setValidator(m_input_validator_double);
    ui->input_polybez_a_ctrl_point_2_z->setValidator(m_input_validator_double);
    ui->input_polybez_a_ctrl_point_3_x->setValidator(m_input_validator_double);
    ui->input_polybez_a_ctrl_point_3_y->setValidator(m_input_validator_double);
    ui->input_polybez_a_ctrl_point_3_z->setValidator(m_input_validator_double);
    // ------------------------------------------------------------------------

    // ------------------------------------------------------------------------
    // polybezier b input formattings
    ui->input_polybez_b_normal_x->setValidator(m_input_validator_double);
    ui->input_polybez_b_normal_y->setValidator(m_input_validator_double);
    ui->input_polybez_b_normal_z->setValidator(m_input_validator_double);
    ui->input_polybez_b_ctrl_point_0_x->setValidator(m_input_validator_double);
    ui->input_polybez_b_ctrl_point_0_y->setValidator(m_input_validator_double);
    ui->input_polybez_b_ctrl_point_0_z->setValidator(m_input_validator_double);
    ui->input_polybez_b_ctrl_point_1_x->setValidator(m_input_validator_double);
    ui->input_polybez_b_ctrl_point_1_y->setValidator(m_input_validator_double);
    ui->input_polybez_b_ctrl_point_1_z->setValidator(m_input_validator_double);
    ui->input_polybez_b_ctrl_point_2_x->setValidator(m_input_validator_double);
    ui->input_polybez_b_ctrl_point_2_y->setValidator(m_input_validator_double);
    ui->input_polybez_b_ctrl_point_2_z->setValidator(m_input_validator_double);
    ui->input_polybez_b_ctrl_point_3_x->setValidator(m_input_validator_double);
    ui->input_polybez_b_ctrl_point_3_y->setValidator(m_input_validator_double);
    ui->input_polybez_b_ctrl_point_3_z->setValidator(m_input_validator_double);
    // ------------------------------------------------------------------------

    // ------------------------------------------------------------------------
    // intersection input formattings
    ui->input_isec_position_x->setValidator(m_input_validator_double);
    ui->input_isec_position_y->setValidator(m_input_validator_double);
    ui->input_isec_position_z->setValidator(m_input_validator_double);
    ui->input_isec_normal_x->setValidator(m_input_validator_double);
    ui->input_isec_normal_y->setValidator(m_input_validator_double);
    ui->input_isec_normal_z->setValidator(m_input_validator_double);
    ui->input_isec_tval_a->setValidator(m_input_validator_double);
    ui->input_isec_tval_b->setValidator(m_input_validator_double);
    ui->input_isec_tan_a_x->setValidator(m_input_validator_double);
    ui->input_isec_tan_a_y->setValidator(m_input_validator_double);
    ui->input_isec_tan_a_z->setValidator(m_input_validator_double);
    ui->input_isec_tan_b_x->setValidator(m_input_validator_double);
    ui->input_isec_tan_b_y->setValidator(m_input_validator_double);
    ui->input_isec_tan_b_z->setValidator(m_input_validator_double);
    // ------------------------------------------------------------------------

    // ------------------------------------------------------------------------
    // polybezier a signals & slots connections to listen for parameter edits
    connect(ui->cmbobox_polybez_a_curve_type, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this] ( ) { send_polybez_form_info(0); });
    connect(ui->input_polybez_a_normal_x, &QLineEdit::editingFinished, this, [this] ( ) { send_polybez_form_info(1); });
    connect(ui->input_polybez_a_normal_y, &QLineEdit::editingFinished, this, [this] ( ) { send_polybez_form_info(1); });
    connect(ui->input_polybez_a_normal_z, &QLineEdit::editingFinished, this, [this] ( ) { send_polybez_form_info(1); });
    // offset
    connect(ui->input_polybez_a_ctrl_point_0_x, &QLineEdit::editingFinished, this, [this] ( ) { send_polybez_form_info(3); });
    connect(ui->input_polybez_a_ctrl_point_0_y, &QLineEdit::editingFinished, this, [this] ( ) { send_polybez_form_info(3); });
    connect(ui->input_polybez_a_ctrl_point_0_z, &QLineEdit::editingFinished, this, [this] ( ) { send_polybez_form_info(3); });
    connect(ui->input_polybez_a_ctrl_point_1_x, &QLineEdit::editingFinished, this, [this] ( ) { send_polybez_form_info(4); });
    connect(ui->input_polybez_a_ctrl_point_1_y, &QLineEdit::editingFinished, this, [this] ( ) { send_polybez_form_info(4); });
    connect(ui->input_polybez_a_ctrl_point_1_z, &QLineEdit::editingFinished, this, [this] ( ) { send_polybez_form_info(4); });
    connect(ui->input_polybez_a_ctrl_point_2_x, &QLineEdit::editingFinished, this, [this] ( ) { send_polybez_form_info(5); });
    connect(ui->input_polybez_a_ctrl_point_2_y, &QLineEdit::editingFinished, this, [this] ( ) { send_polybez_form_info(5); });
    connect(ui->input_polybez_a_ctrl_point_2_z, &QLineEdit::editingFinished, this, [this] ( ) { send_polybez_form_info(5); });
    connect(ui->input_polybez_a_ctrl_point_3_x, &QLineEdit::editingFinished, this, [this] ( ) { send_polybez_form_info(6); });
    connect(ui->input_polybez_a_ctrl_point_3_y, &QLineEdit::editingFinished, this, [this] ( ) { send_polybez_form_info(6); });
    connect(ui->input_polybez_a_ctrl_point_3_z, &QLineEdit::editingFinished, this, [this] ( ) { send_polybez_form_info(6); });
    // ------------------------------------------------------------------------
}

Widget_Selection_Info::~Widget_Selection_Info()
{
    delete ui;
    delete m_input_validator_double;
}

void Widget_Selection_Info::write_polybez_a_form_info(polybezier_info_form_s& in_form)
{
    this->blockSignals(true);
    ui->spnbox_polybez_a_segment->disconnect();
    connect(ui->spnbox_polybez_a_segment, QOverload<int>::of(&QSpinBox::valueChanged), this, [this] (int segment) { emit signal_spnbox_polybez_a_segment(segment); });

    ui->label_polybez_a_id->setText(QString::fromStdString(boost::uuids::to_string(in_form.uuid)));

    ui->cmbobox_polybez_a_curve_type->setCurrentIndex(in_form.curve_type);

    ui->input_polybez_a_normal_x->setText(QString::number(in_form.plane_normal.x, 'f', 4));
    ui->input_polybez_a_normal_y->setText(QString::number(in_form.plane_normal.y, 'f', 4));
    ui->input_polybez_a_normal_z->setText(QString::number(in_form.plane_normal.z, 'f', 4));

    //ui->input_polybez_a_offset->setText(QString::number(plane_offset, 'f', 4));

    ui->input_polybez_a_ctrl_point_0_x->setText(QString::number(in_form.ctrl_points[0].x, 'f', 4));
    ui->input_polybez_a_ctrl_point_0_y->setText(QString::number(in_form.ctrl_points[0].y, 'f', 4));
    ui->input_polybez_a_ctrl_point_0_z->setText(QString::number(in_form.ctrl_points[0].z, 'f', 4));

    ui->input_polybez_a_ctrl_point_1_x->setText(QString::number(in_form.ctrl_points[1].x, 'f', 4));
    ui->input_polybez_a_ctrl_point_1_y->setText(QString::number(in_form.ctrl_points[1].y, 'f', 4));
    ui->input_polybez_a_ctrl_point_1_z->setText(QString::number(in_form.ctrl_points[1].z, 'f', 4));

    ui->input_polybez_a_ctrl_point_2_x->setText(QString::number(in_form.ctrl_points[2].x, 'f', 4));
    ui->input_polybez_a_ctrl_point_2_y->setText(QString::number(in_form.ctrl_points[2].y, 'f', 4));
    ui->input_polybez_a_ctrl_point_2_z->setText(QString::number(in_form.ctrl_points[2].z, 'f', 4));

    ui->input_polybez_a_ctrl_point_3_x->setText(QString::number(in_form.ctrl_points[3].x, 'f', 4));
    ui->input_polybez_a_ctrl_point_3_y->setText(QString::number(in_form.ctrl_points[3].y, 'f', 4));
    ui->input_polybez_a_ctrl_point_3_z->setText(QString::number(in_form.ctrl_points[3].z, 'f', 4));

    ui->spnbox_polybez_a_segment->setRange(0, in_form.num_segments - 1);
    ui->spnbox_polybez_a_segment->setValue(in_form.curr_segment);

    ui->group_polybez_a_info->setEnabled(true);
    this->blockSignals(false);
}

void Widget_Selection_Info::write_polybez_b_form_info(polybezier_info_form_s& in_form)
{
    this->blockSignals(true);
    ui->spnbox_polybez_b_segment->disconnect();
    connect(ui->spnbox_polybez_b_segment, QOverload<int>::of(&QSpinBox::valueChanged), this, [this] (int segment) { emit signal_spnbox_polybez_b_segment(segment); });

    ui->label_polybez_b_id->setText(QString::fromStdString(boost::uuids::to_string(in_form.uuid)));

    ui->cmbobox_polybez_b_curve_type->setCurrentIndex(in_form.curve_type);

    ui->input_polybez_b_normal_x->setText(QString::number(in_form.plane_normal.x, 'f', 4));
    ui->input_polybez_b_normal_y->setText(QString::number(in_form.plane_normal.y, 'f', 4));
    ui->input_polybez_b_normal_z->setText(QString::number(in_form.plane_normal.z, 'f', 4));

    //ui->input_polybez_b_offset->setText(QString::number(plane_offset, 'f', 4));

    ui->input_polybez_b_ctrl_point_0_x->setText(QString::number(in_form.ctrl_points[0].x, 'f', 4));
    ui->input_polybez_b_ctrl_point_0_y->setText(QString::number(in_form.ctrl_points[0].y, 'f', 4));
    ui->input_polybez_b_ctrl_point_0_z->setText(QString::number(in_form.ctrl_points[0].z, 'f', 4));

    ui->input_polybez_b_ctrl_point_1_x->setText(QString::number(in_form.ctrl_points[1].x, 'f', 4));
    ui->input_polybez_b_ctrl_point_1_y->setText(QString::number(in_form.ctrl_points[1].y, 'f', 4));
    ui->input_polybez_b_ctrl_point_1_z->setText(QString::number(in_form.ctrl_points[1].z, 'f', 4));

    ui->input_polybez_b_ctrl_point_2_x->setText(QString::number(in_form.ctrl_points[1].x, 'f', 4));
    ui->input_polybez_b_ctrl_point_2_y->setText(QString::number(in_form.ctrl_points[1].y, 'f', 4));
    ui->input_polybez_b_ctrl_point_2_z->setText(QString::number(in_form.ctrl_points[1].z, 'f', 4));

    ui->input_polybez_b_ctrl_point_3_x->setText(QString::number(in_form.ctrl_points[1].x, 'f', 4));
    ui->input_polybez_b_ctrl_point_3_y->setText(QString::number(in_form.ctrl_points[1].y, 'f', 4));
    ui->input_polybez_b_ctrl_point_3_z->setText(QString::number(in_form.ctrl_points[1].z, 'f', 4));

    ui->spnbox_polybez_b_segment->setRange(0, in_form.num_segments - 1);
    ui->spnbox_polybez_b_segment->setValue(in_form.curr_segment);

    ui->group_polybez_b_info->setEnabled(true);
    this->blockSignals(false);
}

void Widget_Selection_Info::write_intersection_form_info(intersection_info_form_s &form)
{
    this->blockSignals(true);
    ui->label_isec_id->setText(QString::fromStdString(boost::uuids::to_string(form.uuid)));

    ui->input_isec_position_x->setText(QString::number(form.position.x, 'f', 4));
    ui->input_isec_position_y->setText(QString::number(form.position.y, 'f', 4));
    ui->input_isec_position_z->setText(QString::number(form.position.z, 'f', 4));

    ui->input_isec_normal_x->setText(QString::number(form.normal.x, 'f', 4));
    ui->input_isec_normal_y->setText(QString::number(form.normal.y, 'f', 4));
    ui->input_isec_normal_z->setText(QString::number(form.normal.z, 'f', 4));

    ui->input_isec_tval_a->setText(QString::number(form.tval_a, 'f', 4));
    ui->input_isec_tval_b->setText(QString::number(form.tval_b, 'f', 4));

    ui->input_isec_tan_a_x->setText(QString::number(form.tangent_a.x, 'f', 4));
    ui->input_isec_tan_a_y->setText(QString::number(form.tangent_a.y, 'f', 4));
    ui->input_isec_tan_a_z->setText(QString::number(form.tangent_a.z, 'f', 4));
    ui->input_isec_tan_b_x->setText(QString::number(form.tangent_b.x, 'f', 4));
    ui->input_isec_tan_b_y->setText(QString::number(form.tangent_b.y, 'f', 4));
    ui->input_isec_tan_b_z->setText(QString::number(form.tangent_b.z, 'f', 4));

    ui->group_isec_info->setEnabled(true);
    this->blockSignals(false);
}

void Widget_Selection_Info::send_polybez_form_info(int flag)
{
    if (!ui->group_polybez_a_info->isEnabled()) {
        return;
    }

    switch (flag)
    {
    case 0: { // curve type changed
        emit request_selected_polybez_curve_type_update((curve_type_t)ui->cmbobox_polybez_a_curve_type->currentIndex());
        break;
    }
    case 1: { // plane normal changed
        vector_t normal = vector_t(ui->input_polybez_a_normal_x->text().toDouble(),
                                   ui->input_polybez_a_normal_y->text().toDouble(),
                                   ui->input_polybez_a_normal_z->text().toDouble());
        emit request_selected_polybez_plane_normal_update(normal);
        break;
    }
    case 2: { // plane offset changed
        emit request_selected_polybez_plane_offset_update(0.0);
        break;
    }
    case 3: { // ctrl point 0 changed
        point_t ctrl_point = point_t(ui->input_polybez_a_ctrl_point_0_x->text().toDouble(),
                                     ui->input_polybez_a_ctrl_point_0_y->text().toDouble(),
                                     ui->input_polybez_a_ctrl_point_0_z->text().toDouble());
        emit request_selected_polybez_ctrl_point_update(ui->spnbox_polybez_a_segment->value(), 0, ctrl_point);
        break;
    }
    case 4: { // ctrl point 1 changed
        point_t ctrl_point = point_t(ui->input_polybez_a_ctrl_point_1_x->text().toDouble(),
                                     ui->input_polybez_a_ctrl_point_1_y->text().toDouble(),
                                     ui->input_polybez_a_ctrl_point_1_z->text().toDouble());
        emit request_selected_polybez_ctrl_point_update(ui->spnbox_polybez_a_segment->value(), 1, ctrl_point);
        break;
    }
    case 5: { // ctrl point 2 changed
        point_t ctrl_point = point_t(ui->input_polybez_a_ctrl_point_2_x->text().toDouble(),
                                     ui->input_polybez_a_ctrl_point_2_y->text().toDouble(),
                                     ui->input_polybez_a_ctrl_point_2_z->text().toDouble());
        emit request_selected_polybez_ctrl_point_update(ui->spnbox_polybez_a_segment->value(), 2, ctrl_point);
        break;
    }
    case 6: { // ctrl point 3 changed
        point_t ctrl_point = point_t(ui->input_polybez_a_ctrl_point_3_x->text().toDouble(),
                                     ui->input_polybez_a_ctrl_point_3_y->text().toDouble(),
                                     ui->input_polybez_a_ctrl_point_3_z->text().toDouble());
        emit request_selected_polybez_ctrl_point_update(ui->spnbox_polybez_a_segment->value(), 3, ctrl_point);
        break;
    }
    default: {
        break;
    }
    }
}

void Widget_Selection_Info::clear_polybez_a_form_info( )
{
    this->blockSignals(true);

    ui->group_polybez_a_info->setDisabled(true);

    ui->label_polybez_a_id->setText("-1");

    ui->cmbobox_polybez_a_curve_type->setCurrentIndex(0);

    ui->input_polybez_a_normal_x->clear();
    ui->input_polybez_a_normal_y->clear();
    ui->input_polybez_a_normal_z->clear();

    ui->input_polybez_a_ctrl_point_0_x->clear();
    ui->input_polybez_a_ctrl_point_0_y->clear();
    ui->input_polybez_a_ctrl_point_0_z->clear();

    ui->input_polybez_a_ctrl_point_1_x->clear();
    ui->input_polybez_a_ctrl_point_1_y->clear();
    ui->input_polybez_a_ctrl_point_1_z->clear();

    ui->input_polybez_a_ctrl_point_2_x->clear();
    ui->input_polybez_a_ctrl_point_2_y->clear();
    ui->input_polybez_a_ctrl_point_2_z->clear();

    ui->input_polybez_a_ctrl_point_3_x->clear();
    ui->input_polybez_a_ctrl_point_3_y->clear();
    ui->input_polybez_a_ctrl_point_3_z->clear();

    ui->spnbox_polybez_a_segment->disconnect();

    ui->spnbox_polybez_a_segment->setRange(-1, -1);
    ui->spnbox_polybez_a_segment->setValue(-1);

    this->blockSignals(false);
}

void Widget_Selection_Info::clear_polybez_b_form_info( )
{
    this->blockSignals(true);

    ui->group_polybez_b_info->setDisabled(true);

    ui->label_polybez_b_id->setText("-1");

    ui->cmbobox_polybez_b_curve_type->setCurrentIndex(0);

    ui->input_polybez_b_normal_x->clear();
    ui->input_polybez_b_normal_y->clear();
    ui->input_polybez_b_normal_z->clear();

    ui->input_polybez_b_ctrl_point_0_x->clear();
    ui->input_polybez_b_ctrl_point_0_y->clear();
    ui->input_polybez_b_ctrl_point_0_z->clear();

    ui->input_polybez_b_ctrl_point_1_x->clear();
    ui->input_polybez_b_ctrl_point_1_y->clear();
    ui->input_polybez_b_ctrl_point_1_z->clear();

    ui->input_polybez_b_ctrl_point_2_x->clear();
    ui->input_polybez_b_ctrl_point_2_y->clear();
    ui->input_polybez_b_ctrl_point_2_z->clear();

    ui->input_polybez_b_ctrl_point_3_x->clear();
    ui->input_polybez_b_ctrl_point_3_y->clear();
    ui->input_polybez_b_ctrl_point_3_z->clear();

    ui->spnbox_polybez_b_segment->disconnect();

    ui->spnbox_polybez_b_segment->setRange(-1, -1);
    ui->spnbox_polybez_b_segment->setValue(-1);

    this->blockSignals(false);
}

void Widget_Selection_Info::clear_intersection_form_info()
{
    this->blockSignals(true);

    ui->group_isec_info->setDisabled(true);

    ui->label_isec_id->setText("-1");

    ui->input_isec_position_x->clear();
    ui->input_isec_position_y->clear();
    ui->input_isec_position_z->clear();

    ui->input_isec_normal_x->clear();
    ui->input_isec_normal_y->clear();
    ui->input_isec_normal_z->clear();

    ui->label_isec_tval_a_polybez_id->setText("-1");
    ui->label_isec_tval_b_polybez_id->setText("-1");

    ui->input_isec_tval_a->clear();
    ui->input_isec_tval_b->clear();

    ui->label_isec_tan_a_polybez_id->setText("-1");
    ui->label_isec_tan_b_polybez_id->setText("-1");

    ui->input_isec_tan_a_x->clear();
    ui->input_isec_tan_a_y->clear();
    ui->input_isec_tan_a_z->clear();
    ui->input_isec_tan_b_x->clear();
    ui->input_isec_tan_b_y->clear();
    ui->input_isec_tan_b_z->clear();

    this->blockSignals(false);
}

void Widget_Selection_Info::clear_all_forms_info()
{
    clear_polybez_a_form_info();
    clear_polybez_b_form_info();
    clear_intersection_form_info();
}
