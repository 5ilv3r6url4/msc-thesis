#ifndef DIALOG_NEW_FILE_H
#define DIALOG_NEW_FILE_H

#include <QWidget>
#include <QShortcut>
#include <QApplication>
#include <QGridLayout>
#include <QPushButton>
#include <QLabel>
#include <QLineEdit>
#include <QIntValidator>

class Dialog_New_File : public QWidget
{
    Q_OBJECT

private slots:
    void slot_ok();
    void slot_cancel();

signals:
    void request_new_file(QString, int, int);

public:
    Dialog_New_File(QWidget *parent = 0);
    ~Dialog_New_File();

private:

    QPushButton*    m_button_ok;
    QPushButton*    m_button_cancel;

    QLabel*         m_label_filename;   QLineEdit*      m_input_filename;
    QLabel*         m_label_width;      QLineEdit*      m_input_width;
    QLabel*         m_label_height;     QLineEdit*      m_input_height;

};

#endif // DIALOG_NEW_FILE_H
