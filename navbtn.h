#ifndef NAVBTN_H
#define NAVBTN_H

#include <QWidget>

//左侧导航栏按钮
namespace Ui {
class NavBtn;
}

class NavBtn : public QWidget
{
    Q_OBJECT

public:
    explicit NavBtn(QWidget *parent = nullptr);
    ~NavBtn();
    void setPicName(QString pic,QString Name);
    QString getName();
    void setIndex(int index);
    void setNormalStyle();
    void setClickedStyle();
    int getIndex();
signals:
    void sigClicked(QString name);
protected:
    bool eventFilter(QObject *watched, QEvent *event);

private:
    Ui::NavBtn *ui;
    int m_index;
};

#endif // NAVBTN_H
