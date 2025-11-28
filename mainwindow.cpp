#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "navbtn.h"
#include "monitorwidget.h"
#include "dataeditorwidget.h"
#include "modelmanager.h"
#include "plottingwidget.h"
#include "fittingwidget.h"
#include "settingswidget.h"

#include <QDateTime>
#include <QMessageBox>
#include <QDebug>
#include <QRegularExpression>
#include <QStandardItemModel>
#include <QTimer>
#include <QSpacerItem>
#include <QStackedWidget>
#include <cmath>
#include <QStatusBar>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    // 设置窗口标题
    this->setWindowTitle("陆相泥纹型及混积型页岩油压裂水平井非均匀产液机制与试井解释方法研究");

    // 由于标题较长，适当增加窗口最小宽度
    this->setMinimumWidth(1024);

    init();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::init()
{
    // 初始化导航栏按钮
    for(int i = 0 ; i<6;i++)  // 6个按钮（项目、数据、模型、图表、拟合、设置）
    {
        NavBtn* btn = new NavBtn(ui->widgetNav);
        btn->setMinimumWidth(110);
        btn->setIndex(i);
        // 设置字体颜色为黑色
        btn->setStyleSheet("color: black;");

        switch (i) {
        case 0:
            btn->setPicName("border-image: url(:/new/prefix1/Resource/X0.png);",tr("项目"));
            btn->setClickedStyle();
            ui->stackedWidget->setCurrentIndex(0);
            break;
        case 1:
            btn->setPicName("border-image: url(:/new/prefix1/Resource/X1.png);",tr("数据"));
            break;
        case 2:
            btn->setPicName("border-image: url(:/new/prefix1/Resource/X2.png);",tr("模型"));
            break;
        case 3:
            btn->setPicName("border-image: url(:/new/prefix1/Resource/X3.png);",tr("图表"));
            break;
        case 4:
            btn->setPicName("border-image: url(:/new/prefix1/Resource/X4.png);",tr("拟合"));
            break;
        case 5:
            btn->setPicName("border-image: url(:/new/prefix1/Resource/X5.png);",tr("设置"));
            break;
        default:
            break;
        }
        m_NavBtnMap.insert(btn->getName(),btn);
        ui->verticalLayoutNav->addWidget(btn);
        connect(btn,&NavBtn::sigClicked,[=](QString name)
                {
                    QMap<QString,NavBtn*>::Iterator item = m_NavBtnMap.begin();
                    while (item != m_NavBtnMap.end()) {
                        if(item.key() != name)
                        {
                            ((NavBtn*)(item.value()))->setNormalStyle();
                        }
                        item++;
                    }
                    int targetIndex = m_NavBtnMap.value(name)->getIndex();
                    ui->stackedWidget->setCurrentIndex(targetIndex);

                    // --- 特殊处理：切换页面时的数据传递 ---
                    if (name == tr("图表")) {
                        onTransferDataToPlotting();
                    }
                    else if (name == tr("拟合")) {
                        transferDataToFitting();
                    }
                });
    }

    // 添加垂直弹性空间，防止按钮拉伸
    QSpacerItem* verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);
    ui->verticalLayoutNav->addSpacerItem(verticalSpacer);

    // 初始化时间显示
    ui->labelTime->setText(QDateTime::currentDateTime().toString("yyyy-MM-dd hh-mm-ss").replace(" ","\n"));
    connect(&m_timer,&QTimer::timeout,[=]
            {
                ui->labelTime->setText(QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").replace(" ","\n"));
                ui->labelTime->setStyleSheet("color: black;");
            });
    m_timer.start(1000);

    // 1. 初始化监控界面
    m_MonitorWidget = new MonitorWidget(ui->pageMonitor);
    ui->verticalLayoutMonitor->addWidget(m_MonitorWidget);
    connect(m_MonitorWidget, &MonitorWidget::newProjectCreated, this, &MainWindow::onProjectCreated);
    connect(m_MonitorWidget, &MonitorWidget::fileLoaded, this, &MainWindow::onFileLoaded);

    // 2. 创建数据编辑器并将其添加到数据页面
    m_DataEditorWidget = new DataEditorWidget(ui->pageHand);
    ui->verticalLayoutHandle->addWidget(m_DataEditorWidget);
    connect(m_DataEditorWidget, &DataEditorWidget::fileChanged, this, &MainWindow::onFileLoaded);
    connect(m_DataEditorWidget, &DataEditorWidget::dataChanged, this, &MainWindow::onDataEditorDataChanged);

    // 3. 创建并初始化模型管理器
    m_ModelManager = new ModelManager(this);
    m_ModelManager->initializeModels(ui->pageParamter);
    connect(m_ModelManager, &ModelManager::calculationCompleted,
            this, &MainWindow::onModelCalculationCompleted);

    // 4. 创建并初始化绘图界面
    m_PlottingWidget = new PlottingWidget(ui->pageData);
    ui->verticalLayout_2->addWidget(m_PlottingWidget);
    connect(m_PlottingWidget, &PlottingWidget::analysisCompleted,
            this, &MainWindow::onPlotAnalysisCompleted);

    // 5. 创建并初始化拟合界面
    if (ui->pageFitting && ui->verticalLayoutFitting) {
        m_FittingWidget = new FittingWidget(ui->pageFitting);
        ui->verticalLayoutFitting->addWidget(m_FittingWidget);

        connect(m_FittingWidget, &FittingWidget::fittingCompleted,
                this, &MainWindow::onFittingCompleted);
        // 连接进度信号
        connect(m_FittingWidget, &FittingWidget::sigProgress,
                this, &MainWindow::onFittingProgressChanged);
    } else {
        qWarning() << "MainWindow: pageFitting或verticalLayoutFitting为空！无法创建拟合界面";
        m_FittingWidget = nullptr;
    }

    // 6. 创建并初始化设置界面
    m_SettingsWidget = new SettingsWidget(ui->pageAlarm);
    ui->verticalLayout_3->addWidget(m_SettingsWidget);

    connect(m_SettingsWidget, &SettingsWidget::systemSettingsChanged,
            this, &MainWindow::onSystemSettingsChanged);
    connect(m_SettingsWidget, &SettingsWidget::autoSaveIntervalChanged,
            this, &MainWindow::onAutoSaveIntervalChanged);
    connect(m_SettingsWidget, &SettingsWidget::backupSettingsChanged,
            this, &MainWindow::onBackupSettingsChanged);

    // 调用各子模块的初始化函数
    initMonitorForm();
    initDataEditorForm();
    initModelForm();
    initPlottingForm();
    initFittingForm();
}

void MainWindow::initMonitorForm()
{
    qDebug() << "初始化监控界面";
}

void MainWindow::initDataEditorForm()
{
    qDebug() << "初始化数据编辑器界面";
}

void MainWindow::initModelForm()
{
    if (m_ModelManager) {
        qDebug() << "模型界面初始化完成";
    }
}

void MainWindow::initPlottingForm()
{
    qDebug() << "初始化绘图界面";
}

void MainWindow::initFittingForm()
{
    // 注入 ModelManager 依赖，使 FittingWidget 能调用模型计算
    if (m_FittingWidget && m_ModelManager) {
        m_FittingWidget->setModelManager(m_ModelManager);
        qDebug() << "拟合界面初始化完成，依赖已注入";
    }
}

void MainWindow::onProjectCreated()
{
    qDebug() << "处理新项目创建";
    QMessageBox::information(this, tr("项目创建"), tr("新项目已创建成功！"));
    ui->stackedWidget->setCurrentIndex(2);
    updateNavigationState();
}

void MainWindow::onFileLoaded(const QString& filePath, const QString& fileType)
{
    qDebug() << "MainWindow收到文件加载/更换通知：" << filePath;
    ui->stackedWidget->setCurrentIndex(1);

    // 更新导航按钮状态
    QMap<QString,NavBtn*>::Iterator item = m_NavBtnMap.begin();
    while (item != m_NavBtnMap.end()) {
        ((NavBtn*)(item.value()))->setNormalStyle();
        if(item.key() == tr("数据")) {
            ((NavBtn*)(item.value()))->setClickedStyle();
        }
        item++;
    }

    if (m_DataEditorWidget && sender() != m_DataEditorWidget) {
        m_DataEditorWidget->loadData(filePath, fileType);
    }

    m_hasValidData = true;
    QTimer::singleShot(1000, this, &MainWindow::onDataReadyForPlotting);
}

void MainWindow::onPlotAnalysisCompleted(const QString &analysisType, const QMap<QString, double> &results)
{
    qDebug() << "绘图分析完成：" << analysisType;
}

void MainWindow::onDataReadyForPlotting()
{
    transferDataFromEditorToPlotting();
}

void MainWindow::onTransferDataToPlotting()
{
    if (!hasDataLoaded()) {
        return;
    }
    transferDataFromEditorToPlotting();
}

void MainWindow::onDataEditorDataChanged()
{
    // 如果当前在图表页或拟合页，数据变动时实时更新
    if (ui->stackedWidget->currentIndex() == 3) {
        transferDataFromEditorToPlotting();
    }
    m_hasValidData = hasDataLoaded();
}

void MainWindow::onModelCalculationCompleted(const QString &analysisType, const QMap<QString, double> &results)
{
    qDebug() << "模型计算完成：" << analysisType;
}

// -------------------------------------------------------------------------
// 拟合模块数据传递逻辑 (关键修改部分)
// -------------------------------------------------------------------------

void MainWindow::transferDataToFitting()
{
    if (!m_FittingWidget || !m_DataEditorWidget) return;

    qDebug() << "正在同步数据至拟合模块...";

    QStandardItemModel* model = m_DataEditorWidget->getDataModel();

    // [修改] 如果没有数据，传递空向量以清空图表，绝不生成演示数据
    if (!model || model->rowCount() == 0) {
        m_FittingWidget->setObservedData(QVector<double>(), QVector<double>(), QVector<double>());
        return;
    }

    QVector<double> tVec, pVec, dVec;
    double p_initial = 0.0;
    bool p_init_found = false;

    // 1. 扫描初始压力 (P0)
    // 假设第2列(index 1)是压力列
    for(int r=0; r<model->rowCount(); ++r) {
        QModelIndex idx = model->index(r, 1);
        if(idx.isValid()) {
            double p = idx.data().toDouble();
            if (std::abs(p) > 1e-6) { // 找到第一个非零压力作为初始压力
                p_initial = p;
                p_init_found = true;
                break;
            }
        }
    }

    // 2. 读取数据并转换为 Delta P
    for(int r=0; r<model->rowCount(); ++r) {
        double t = model->index(r, 0).data().toDouble(); // 第1列：时间
        double p_raw = model->index(r, 1).data().toDouble(); // 第2列：压力

        if (t > 0) {
            tVec.append(t);
            // [核心] 转换为压力差: Delta P = |P - P_initial|
            pVec.append(std::abs(p_raw - p_initial));
        }
    }

    // 3. 自动计算导数 (如果 DataEditor 没有提供)
    // 假设 DataEditor 目前没有计算导数，我们在这里计算 Bourdet 导数
    dVec.resize(tVec.size());
    if (tVec.size() > 2) {
        dVec[0] = 0;
        dVec[tVec.size()-1] = 0;

        for(int i=1; i<tVec.size()-1; ++i) {
            double lnt1 = std::log(tVec[i-1]);
            double lnt2 = std::log(tVec[i]);
            double lnt3 = std::log(tVec[i+1]);

            // 防止时间点重复导致除零
            if (std::abs(lnt2 - lnt1) < 1e-9 || std::abs(lnt3 - lnt2) < 1e-9) {
                dVec[i] = 0;
                continue;
            }

            // 对数时间的差分
            double d1 = (pVec[i] - pVec[i-1]) / (lnt2 - lnt1);
            double d2 = (pVec[i+1] - pVec[i]) / (lnt3 - lnt2);

            // Bourdet 加权平均
            double w1 = (lnt3 - lnt2) / (lnt3 - lnt1);
            double w2 = (lnt2 - lnt1) / (lnt3 - lnt1);

            dVec[i] = d1 * w1 + d2 * w2;
        }
    }

    // 将处理好的一套数据（时间、压差、导数）传递给拟合模块
    m_FittingWidget->setObservedData(tVec, pVec, dVec);
}

void MainWindow::onFittingCompleted(ModelManager::ModelType modelType, const QMap<QString, double> &parameters)
{
    QString typeName = ModelManager::getModelTypeName(modelType);
    QString msg = QString("拟合完成！\n模型: %1\n\n最优参数:\n").arg(typeName);

    QMapIterator<QString, double> i(parameters);
    while (i.hasNext()) {
        i.next();
        msg += QString("%1: %2\n").arg(i.key()).arg(i.value(), 0, 'f', 6);
    }
    QMessageBox::information(this, "拟合结果", msg);
}

void MainWindow::onFittingProgressChanged(int progress)
{
    if (this->statusBar()) {
        this->statusBar()->showMessage(QString("正在拟合... %1%").arg(progress));
        if(progress >= 100) this->statusBar()->showMessage("拟合完成", 5000);
    }
}

// -------------------------------------------------------------------------

void MainWindow::onSystemSettingsChanged()
{
    qDebug() << "系统设置已更改";
}

void MainWindow::onAutoSaveIntervalChanged(int interval)
{
    qDebug() << "自动保存间隔已更改为：" << interval << "分钟";
    static QTimer* autoSaveTimer = nullptr;
    if (!autoSaveTimer) {
        autoSaveTimer = new QTimer(this);
        connect(autoSaveTimer, &QTimer::timeout, this, [this]() {
            qDebug() << "执行自动保存...";
        });
    }
    autoSaveTimer->stop();
    if (interval > 0) {
        autoSaveTimer->start(interval * 60 * 1000);
    }
}

void MainWindow::onBackupSettingsChanged(bool enabled)
{
    qDebug() << "备份设置已更改：" << (enabled ? "启用" : "禁用");
}

void MainWindow::onPerformanceSettingsChanged()
{
    qDebug() << "性能设置已更改";
}

QStandardItemModel* MainWindow::getDataEditorModel() const
{
    if (!m_DataEditorWidget) return nullptr;
    return m_DataEditorWidget->getDataModel();
}

QString MainWindow::getCurrentFileName() const
{
    if (!m_DataEditorWidget) return QString();
    return m_DataEditorWidget->getCurrentFileName();
}

bool MainWindow::hasDataLoaded()
{
    if (!m_DataEditorWidget) return false;
    return m_DataEditorWidget->hasData();
}

void MainWindow::transferDataFromEditorToPlotting()
{
    if (!m_DataEditorWidget || !m_PlottingWidget) return;

    QStandardItemModel* model = m_DataEditorWidget->getDataModel();
    if (model && model->rowCount() > 0 && model->columnCount() > 0) {
        QString fileName = m_DataEditorWidget->getCurrentFileName();
        m_PlottingWidget->setTableDataFromModel(model, fileName);
        m_hasValidData = true;
    } else {
        // 如果是绘图界面，暂保留演示数据以展示功能，若需要也可删除
        WellTestData wellData = createDemoWellTestData();
        m_PlottingWidget->setWellTestData(wellData);
        m_hasValidData = true;
    }
}

WellTestData MainWindow::createDemoWellTestData()
{
    WellTestData wellData;
    wellData.wellName = "演示井-001";
    wellData.testType = "压力恢复试井";
    wellData.testDate = QDateTime::currentDateTime();

    int dataPoints = 150;
    for (int i = 0; i < dataPoints; ++i) {
        double time = 0.01 * std::pow(10, i * 4.0 / dataPoints);
        double pressure = 20.0;
        if (time < 0.1) {
            pressure += 3.0 * (1.0 - std::exp(-time * 10));
        }
        else if (time < 10) {
            pressure += 2.5 + 1.5 * std::log10(time);
        }
        else {
            pressure += 2.5 + 1.5 * std::log10(time) + 0.5 * std::log10(time / 10);
        }
        pressure += 0.05 * std::sin(i * 0.3) + 0.02 * (rand() % 100 - 50) / 50.0;
        wellData.time.append(time);
        wellData.pressure.append(pressure);
    }
    return wellData;
}

void MainWindow::updateNavigationState()
{
    QMap<QString,NavBtn*>::Iterator item = m_NavBtnMap.begin();
    while (item != m_NavBtnMap.end()) {
        ((NavBtn*)(item.value()))->setNormalStyle();
        if(item.key() == tr("模型")) {
            ((NavBtn*)(item.value()))->setClickedStyle();
        }
        item++;
    }
}
