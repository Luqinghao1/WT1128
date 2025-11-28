#include "navbtn.h"
#include "ui_navbtn.h"
#include <QEvent>

NavBtn::NavBtn(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::NavBtn)
{
    ui->setupUi(this);
    ui->widget->installEventFilter(this);
    m_index = 0;

    // 设置按钮固定高度，防止拉伸
    this->setFixedHeight(80);
    this->setMinimumWidth(110);
}

NavBtn::~NavBtn()
{
    delete ui;
}

void NavBtn::setPicName(QString pic, QString Name)
{
    ui->labelPic->setStyleSheet(pic);
    ui->labelName->setText(Name);
}

QString NavBtn::getName()
{
    return ui->labelName->text();
}

void NavBtn::setIndex(int index)
{
    m_index = index;
}

void NavBtn::setNormalStyle()
{
    ui->widget->setStyleSheet("#widget{background-color: rgb(0,0,0,0);}");
}

void NavBtn::setClickedStyle()
{
    ui->widget->setStyleSheet("#widget{background-color: rgb(27,45,85,100);}");
}

int NavBtn::getIndex()
{
    return m_index;
}

bool NavBtn::eventFilter(QObject *watched, QEvent *event)
{
    if(event->type() == QEvent::MouseButtonPress)
    {
        emit sigClicked(ui->labelName->text());
        ui->widget->setStyleSheet("#widget{background-color: rgb(27,45,85,100);}");
    }
    return false;
}
