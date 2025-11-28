#ifndef MONITORWIDGET_H
#define MONITORWIDGET_H

#include <QWidget>
#include <QString>
#include "newprojectdialog.h" // 添加这个头文件引用

// 监控界面
namespace Ui {
class MonitorWidget;
}

class MonitorWidget : public QWidget
{
    Q_OBJECT

public:
    explicit MonitorWidget(QWidget *parent = nullptr);
    ~MonitorWidget();

    void init();
    void createNewProject(); // 添加这个函数声明

signals:
    // 新项目创建成功信号
    void newProjectCreated();

    // 文件读取成功信号，携带文件路径和类型
    void fileLoaded(const QString& filePath, const QString& fileType);

private slots:
    // 点击"新建"按钮的槽函数
    void onNewProjectClicked();

    // 点击"读取"按钮的槽函数
    void onLoadFileClicked();

private:
    Ui::MonitorWidget *ui;
};

#endif // MONITORWIDGET_H
