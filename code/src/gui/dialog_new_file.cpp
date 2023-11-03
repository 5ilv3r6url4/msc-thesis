#include "gui/dialog_new_file.h"

Dialog_New_File::Dialog_New_File(QWidget *parent)
    : QWidget(parent)
{
    setWindowTitle(tr("Create New File"));

    m_button_ok     = new QPushButton(tr("&OK"), this);
    m_button_cancel = new QPushButton(tr("&Cancel"), this);

    connect(m_button_ok,        SIGNAL(clicked()), this, SLOT(slot_ok()));
    connect(m_button_cancel,    SIGNAL(clicked()), this, SLOT(slot_cancel()));

    m_input_filename    = new QLineEdit();
    m_input_width       = new QLineEdit();
    m_input_height      = new QLineEdit();

    m_input_filename->setPlaceholderText("...");
    m_input_width->setPlaceholderText("...");
    m_input_height->setPlaceholderText("...");

    m_input_width->setValidator(new QIntValidator(1, 999999, this));
    m_input_height->setValidator(new QIntValidator(1, 999999, this));

    QGridLayout* layout_main = new QGridLayout(this);

    m_label_filename    = new QLabel(tr("Filename:"));
    m_label_width       = new QLabel(tr("Width:"));
    m_label_height      = new QLabel(tr("Height:"));

    layout_main->addWidget(m_label_filename,    1, 0);
    layout_main->addWidget(m_label_width,       2, 0);
    layout_main->addWidget(m_label_height,      3, 0);

    layout_main->addWidget(m_input_filename,    1, 1);
    layout_main->addWidget(m_input_width,       2, 1);
    layout_main->addWidget(m_input_height,      3, 1);

    layout_main->addWidget(new QLabel(tr("")),  10, 1);
    layout_main->addWidget(m_button_cancel,     11, 0);
    layout_main->addWidget(m_button_ok,         11, 1);

    connect(new QShortcut(QKeySequence::Quit, this), &QShortcut::activated, qApp, &QApplication::quit);
}

Dialog_New_File::~Dialog_New_File()
{
    delete m_button_ok;
    delete m_button_cancel;
    delete m_label_filename;
    delete m_label_width;
    delete m_label_height;
    delete m_input_filename;
    delete m_input_width;
    delete m_input_height;
}

void Dialog_New_File::slot_ok() {

    bool temp_valid_input = true;

    if (m_input_filename->text().isEmpty()) {
        temp_valid_input = false;
        m_input_filename->setStyleSheet("background-color: rgba(255, 0, 0, 50);");
    } else { m_input_filename->setStyleSheet("background-color: white;"); }

    if (m_input_width->text().isEmpty()) {
        temp_valid_input = false;
        m_input_width->setStyleSheet("background-color: rgba(255, 0, 0, 50);");
    } else { m_input_width->setStyleSheet("background-color: white;"); }

    if (m_input_height->text().isEmpty()) {
        temp_valid_input = false;
        m_input_height->setStyleSheet("background-color: rgba(255, 0, 0, 50);");
    } else { m_input_height->setStyleSheet("background-color: white;"); }

    if (!temp_valid_input) {
        return;
    }

    emit request_new_file(m_input_filename->text(), m_input_width->text().toInt(), m_input_height->text().toInt());

    this->close();

    m_input_filename->clear();
    m_input_width->clear();
    m_input_height->clear();
}

void Dialog_New_File::slot_cancel()
{
    this->close();

    m_input_filename->clear();
    m_input_width->clear();
    m_input_height->clear();
}
