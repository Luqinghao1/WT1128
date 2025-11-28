#include "monitorwidget.h"
#include "ui_monitorwidget.h"
#include "newprojectdialog.h"
#include <QDebug>
#include <QFileDialog>
#include <QMessageBox>
#include <QApplication>
#include <QPalette>

MonitorWidget::MonitorWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::MonitorWidget)
{
    ui->setupUi(this);
    init();
}

MonitorWidget::~MonitorWidget()
{
    delete ui;
}

void MonitorWidget::init()
{
    qDebug() << "初始化监控界面...";

    // 设置一个空的顶部样式
    QString topPicStyle = "";
    QString topName = "  ";

    // 设置整体背景为透明
    this->setStyleSheet("background-color: transparent;");
    ui->widget_5->setStyleSheet("background-color: transparent;");

    // 调整按钮之间的间距
    ui->gridLayout_3->setHorizontalSpacing(30);
    ui->gridLayout_3->setVerticalSpacing(10);

    // 创建大字体
    QFont bigFont;
    bigFont.setPointSize(16);  // 设置为16点大小，通常默认是8-10点
    bigFont.setBold(true);

    // 创建背景颜色
    QColor backgroundColor(148, 226, 255);

    // 强制设置样式字符串 - 使用更具体的选择器
    QString forceStyle = QString(
        "MonitoStateW { "
        "background-color: rgb(148, 226, 255); "
        "border-radius: 10px; "
        "padding: 10px; "
        "} "
        "MonitoStateW * { "
        "background-color: transparent; "
        "} "
        "MonitoStateW:hover { "
        "background-color: rgb(120, 200, 240); "
        "} "
        "QLabel { "
        "color: #333333; "
        "font-weight: bold; "
        "margin-top: 5px; "
        "background-color: transparent; "
        "}"
        );

    // "新建"按钮设置
    QString centerPicStyle = "border-image: url(:/new/prefix1/Resource/Mon1.png);";
    QString bottomName = "新建";
    ui->MonitState1->setTextInfo(centerPicStyle, topPicStyle, topName, bottomName);
    ui->MonitState1->setFixedSize(128, 160);

    // 设置背景色 - 使用多种方法确保生效
    ui->MonitState1->setStyleSheet(forceStyle);
    ui->MonitState1->setAutoFillBackground(true);
    QPalette palette1 = ui->MonitState1->palette();
    palette1.setColor(QPalette::Window, backgroundColor);
    ui->MonitState1->setPalette(palette1);

    ui->MonitState1->setFont(bigFont);
    ui->MonitState1->setMouseTracking(true);
    connect(ui->MonitState1, SIGNAL(sigClicked()), this, SLOT(onNewProjectClicked()));
    qDebug() << "连接新建按钮信号到槽...";

    // "打开"按钮设置
    centerPicStyle = "border-image: url(:/new/prefix1/Resource/Mon2.png);";
    bottomName = "打开";
    ui->MonitState2->setTextInfo(centerPicStyle, topPicStyle, topName, bottomName);
    ui->MonitState2->setFixedSize(128, 160);

    // 设置背景色
    ui->MonitState2->setStyleSheet(forceStyle);
    ui->MonitState2->setAutoFillBackground(true);
    QPalette palette2 = ui->MonitState2->palette();
    palette2.setColor(QPalette::Window, backgroundColor);
    ui->MonitState2->setPalette(palette2);

    ui->MonitState2->setFont(bigFont);
    ui->MonitState2->setMouseTracking(true);

    // "读取"按钮设置
    centerPicStyle = "border-image: url(:/new/prefix1/Resource/Mon3.png);";
    bottomName = "读取";
    ui->MonitState3->setTextInfo(centerPicStyle, topPicStyle, topName, bottomName);
    ui->MonitState3->setFixedSize(128, 160);

    // 设置背景色
    ui->MonitState3->setStyleSheet(forceStyle);
    ui->MonitState3->setAutoFillBackground(true);
    QPalette palette3 = ui->MonitState3->palette();
    palette3.setColor(QPalette::Window, backgroundColor);
    ui->MonitState3->setPalette(palette3);

    ui->MonitState3->setFont(bigFont);
    ui->MonitState3->setMouseTracking(true);
    connect(ui->MonitState3, SIGNAL(sigClicked()), this, SLOT(onLoadFileClicked()));
    qDebug() << "连接读取按钮信号到槽...";

    // "退出"按钮设置 - 更改为退出而不是保存
    centerPicStyle = "border-image: url(:/new/prefix1/Resource/Mon4.png);";
    bottomName = "退出";
    ui->MonitState4->setTextInfo(centerPicStyle, topPicStyle, topName, bottomName);
    ui->MonitState4->setFixedSize(128, 160);

    // 设置背景色
    ui->MonitState4->setStyleSheet(forceStyle);
    ui->MonitState4->setAutoFillBackground(true);
    QPalette palette4 = ui->MonitState4->palette();
    palette4.setColor(QPalette::Window, backgroundColor);
    ui->MonitState4->setPalette(palette4);

    ui->MonitState4->setFont(bigFont);
    ui->MonitState4->setMouseTracking(true);

    // 连接退出按钮的点击信号
    connect(ui->MonitState4, &MonitoStateW::sigClicked, this, [=]() {
        qDebug() << "退出按钮被点击...";
        QApplication::quit();
    });

    // 调试输出 - 检查样式是否正确设置
    qDebug() << "MonitState1 样式:" << ui->MonitState1->styleSheet();
    qDebug() << "MonitState1 AutoFillBackground:" << ui->MonitState1->autoFillBackground();
}

void MonitorWidget::onNewProjectClicked()
{
    qDebug() << "新建按钮被点击，准备打开新建项目对话框...";

    // 创建并显示新建项目对话框
    NewProjectDialog* dialog = new NewProjectDialog(this);

    // 以模态方式显示对话框
    int result = dialog->exec();

    if (result == QDialog::Accepted) {
        // 用户点击了"完成"按钮
        qDebug() << "用户创建了新项目";
        // 这里可以添加后续处理代码
        emit newProjectCreated();
    } else {
        // 用户取消了操作
        qDebug() << "用户取消了新建项目";
    }

    // 清理对话框
    delete dialog;
}

void MonitorWidget::onLoadFileClicked()
{
    qDebug() << "读取按钮被点击，准备打开文件选择对话框...";

    // 设置文件过滤器，只显示Excel和TXT文件
    QString filter = tr("Excel Files (*.xlsx *.xls);;Text Files (*.txt);;All Files (*.*)");

    // 打开文件选择对话框
    QString filePath = QFileDialog::getOpenFileName(
        this,
        tr("选择要读取的文件"),
        QString(),  // 默认目录
        filter
        );

    // 如果用户取消了选择，filePath将为空
    if (filePath.isEmpty()) {
        qDebug() << "用户取消了文件选择";
        return;
    }

    qDebug() << "用户选择了文件:" << filePath;

    // 确定文件类型
    QString fileType;
    if (filePath.endsWith(".xlsx", Qt::CaseInsensitive) ||
        filePath.endsWith(".xls", Qt::CaseInsensitive)) {
        fileType = "excel";
    }
    else if (filePath.endsWith(".txt", Qt::CaseInsensitive)) {
        fileType = "txt";
    }
    else {
        fileType = "unknown";
    }

    // 发送文件加载信号，包含文件路径和类型
    emit fileLoaded(filePath, fileType);

    // 创建一个具有黑色文本的消息框
    QMessageBox msgBox;
    msgBox.setWindowTitle(tr("文件读取"));
    msgBox.setText(tr("文件已成功读取，正在准备显示数据..."));
    msgBox.setIcon(QMessageBox::Information);
    // 确保消息框中的文本为黑色
    msgBox.setStyleSheet("QLabel{color: black;}");
    msgBox.exec();
}
