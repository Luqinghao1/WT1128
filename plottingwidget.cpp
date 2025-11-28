#include "plottingwidget.h"
#include "ui_plottingwidget.h"
#include <QPaintEvent>
#include <QPainter>
#include <QApplication>
#include <QClipboard>
#include <QPrintPreviewDialog>
#include <QPrinter>
#include <QPixmap>
#include <QDateTime>
#include <QDebug>
#include <QtMath>
#include <QMenuBar>
#include <QToolBar>
#include <QStatusBar>

// 线型转换函数实现
Qt::PenStyle lineStyleToQt(LineStyle style)
{
    switch (style) {
    case LineStyle::Solid: return Qt::SolidLine;
    case LineStyle::Dash: return Qt::DashLine;
    case LineStyle::Dot: return Qt::DotLine;
    case LineStyle::DashDot: return Qt::DashDotLine;
    case LineStyle::DashDotDot: return Qt::DashDotDotLine;
    default: return Qt::SolidLine;
    }
}

QString lineStyleToString(LineStyle style)
{
    switch (style) {
    case LineStyle::Solid: return "实线";
    case LineStyle::Dash: return "虚线";
    case LineStyle::Dot: return "点线";
    case LineStyle::DashDot: return "点划线";
    case LineStyle::DashDotDot: return "双点划线";
    default: return "实线";
    }
}

LineStyle stringToLineStyle(const QString &str)
{
    if (str == "虚线") return LineStyle::Dash;
    if (str == "点线") return LineStyle::Dot;
    if (str == "点划线") return LineStyle::DashDot;
    if (str == "双点划线") return LineStyle::DashDotDot;
    return LineStyle::Solid;
}

// =======================
// PlottingWidget 类实现
// =======================

PlottingWidget::PlottingWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::PlottingWidget),
    m_hasTableData(false),
    m_isDragging(false),
    m_isSelecting(false),
    m_isPanning(false),
    m_isDraggingLegend(false),
    m_legendOffset(0, 0),
    m_zoomFactor(1.0),
    m_zoomFactorX(1.0),
    m_zoomFactorY(1.0),
    m_viewCenter(0, 0),
    m_panOffset(0, 0)
{
    ui->setupUi(this);
    initializeUI();
    setupDefaultSettings();
    setupConnections();

    setMouseTracking(true);
    ui->widget_plot->setMouseTracking(true);
    ui->widget_plot->installEventFilter(this);

    // 美化坐标标签
    m_coordinateLabel = new QLabel(this);
    m_coordinateLabel->setStyleSheet(
        "QLabel { "
        "   background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1, "
        "                              stop: 0 rgba(33, 150, 243, 0.95), "
        "                              stop: 1 rgba(25, 118, 210, 0.95)); "
        "   border: 2px solid #1976D2; "
        "   border-radius: 6px; "
        "   padding: 6px 12px; "
        "   font-size: 9pt; "
        "   font-weight: bold; "
        "   color: white; "
        "}"
        );
    m_coordinateLabel->hide();
}

PlottingWidget::~PlottingWidget()
{
    for (PlotWindow *window : m_plotWindows) {
        if (window) {
            window->close();
            window->deleteLater();
        }
    }
    for (DualPlotWindow *window : m_dualPlotWindows) {
        if (window) {
            window->close();
            window->deleteLater();
        }
    }
    delete ui;
}

bool PlottingWidget::eventFilter(QObject *obj, QEvent *event)
{
    if (obj == ui->widget_plot && event->type() == QEvent::Paint) {
        QPaintEvent *paintEvent = static_cast<QPaintEvent*>(event);
        paintPlotArea(paintEvent);
        return true;
    }
    return QWidget::eventFilter(obj, event);
}

void PlottingWidget::initializeUI()
{
    ui->splitter_main->setStretchFactor(0, 1);
    ui->splitter_main->setStretchFactor(1, 3);

    m_showGridCheck = ui->checkBox_showGrid;
    m_showLegendCheck = ui->checkBox_showLegend;
    m_gridColorBtn = ui->pushButton_gridColor;

    m_curvesListWidget = ui->listWidget_curves;

    // 美化网格颜色按钮
    m_gridColorBtn->setStyleSheet(
        "QPushButton { "
        "   background-color: white; "
        "   border: 2px solid #2196F3; "
        "   color: #2196F3; "
        "   border-radius: 6px; "
        "   padding: 8px 16px; "
        "   font-weight: bold; "
        "}"
        "QPushButton:hover { "
        "   background-color: #E3F2FD; "
        "   border-color: #1976D2; "
        "}"
        "QPushButton:pressed { "
        "   background-color: #BBDEFB; "
        "}"
        );

    setupContextMenu();

    ui->label_dataInfo->setText("📊 数据信息：未加载数据");

    updateControlsFromSettings();
}

void PlottingWidget::setupDefaultSettings()
{
    m_plotSettings.showGrid = true;
    m_plotSettings.logScaleX = false;
    m_plotSettings.logScaleY = false;
    m_plotSettings.backgroundColor = Qt::white;
    m_plotSettings.gridColor = QColor(224, 224, 224);
    m_plotSettings.textColor = Qt::black;
    m_plotSettings.lineWidth = 2;
    m_plotSettings.pointSize = 4;
    m_plotSettings.xAxisTitle = "时间 (小时)";
    m_plotSettings.yAxisTitle = "产量";
    m_plotSettings.plotTitle = "数据曲线";
    m_plotSettings.autoScale = true;
    m_plotSettings.xMin = 0;
    m_plotSettings.xMax = 100;
    m_plotSettings.yMin = 0;
    m_plotSettings.yMax = 100;
    m_plotSettings.showLegend = true;
    m_plotSettings.legendPosition = QPointF(0.8, 0.1);
    m_plotSettings.xAxisType = AxisType::Linear;
    m_plotSettings.yAxisType = AxisType::Linear;
}

void PlottingWidget::setupConnections()
{
    connect(ui->pushButton_addCurve, &QPushButton::clicked, this, &PlottingWidget::onAddCurve);
    connect(ui->pushButton_editCurve, &QPushButton::clicked, this, &PlottingWidget::onEditCurve);
    connect(ui->pushButton_removeCurve, &QPushButton::clicked, this, &PlottingWidget::onRemoveCurve);

    // 修改连接：使用新的压力产量联合绘图功能
    connect(ui->pushButton_pressureProdData, &QPushButton::clicked, this, &PlottingWidget::onPressureProdDataPlot);
    connect(ui->pushButton_pressureDerivative, &QPushButton::clicked, this, &PlottingWidget::onPressureDerivativePlot);

    connect(ui->listWidget_curves, &QListWidget::itemSelectionChanged, this, &PlottingWidget::onCurveSelectionChanged);
    connect(ui->listWidget_curves, &QListWidget::itemDoubleClicked, [this](QListWidgetItem* item) {
        if (!item) return;
        int index = item->data(Qt::UserRole).toInt();
        if (index >= 0 && index < m_curves.size()) {
            m_curves[index].visible = !m_curves[index].visible;
            updateCurvesList();
            updatePlot();
        }
    });

    connect(ui->checkBox_showGrid, &QCheckBox::toggled, this, &PlottingWidget::updatePlot);
    connect(ui->checkBox_showLegend, &QCheckBox::toggled, this, &PlottingWidget::updatePlot);

    connect(ui->pushButton_gridColor, &QPushButton::clicked,
            this, &PlottingWidget::onColorSettingsChanged);
}

void PlottingWidget::setupContextMenu()
{
    m_contextMenu = new QMenu(this);
    m_contextMenu->setStyleSheet(
        "QMenu { "
        "   background-color: white; "
        "   border: 2px solid #2196F3; "
        "   border-radius: 6px; "
        "   padding: 8px; "
        "}"
        "QMenu::item { "
        "   padding: 8px 24px; "
        "   color: #212121; "
        "   border-radius: 4px; "
        "   margin: 2px; "
        "}"
        "QMenu::item:selected { "
        "   background-color: #2196F3; "
        "   color: white; "
        "}"
        "QMenu::separator { "
        "   height: 2px; "
        "   background: #BBDEFB; "
        "   margin: 4px 8px; "
        "}"
        );

    m_dataMenu = m_contextMenu->addMenu("📍 数据标记");
    m_addMarkerAction = m_dataMenu->addAction("➕ 添加标记点");
    m_addAnnotationAction = m_dataMenu->addAction("📝 添加注释");
    m_dataMenu->addSeparator();
    m_removeLastMarkerAction = m_dataMenu->addAction("❌ 删除最后标记");
    m_removeAllMarkersAction = m_dataMenu->addAction("🗑️ 删除所有标记");
    m_removeAllAnnotationsAction = m_dataMenu->addAction("🗑️ 删除所有注释");

    m_contextMenu->addSeparator();

    m_zoomMenu = m_contextMenu->addMenu("🔍 缩放操作");
    m_zoomInAction = m_zoomMenu->addAction("➕ 放大 (+25%)");
    m_zoomOutAction = m_zoomMenu->addAction("➖ 缩小 (-25%)");
    m_zoomFitAction = m_zoomMenu->addAction("📐 适应窗口");
    m_resetZoomAction = m_zoomMenu->addAction("🔄 重置缩放");

    m_zoomMenu->addSeparator();
    m_zoomXInAction = m_zoomMenu->addAction("↔️ 横向放大");
    m_zoomXOutAction = m_zoomMenu->addAction("↔️ 横向缩小");
    m_zoomYInAction = m_zoomMenu->addAction("↕️ 纵向放大");
    m_zoomYOutAction = m_zoomMenu->addAction("↕️ 纵向缩小");

    connect(m_addMarkerAction, &QAction::triggered, this, &PlottingWidget::onMarkerAdded);
    connect(m_addAnnotationAction, &QAction::triggered, this, &PlottingWidget::onAnnotationAdded);
    connect(m_removeLastMarkerAction, &QAction::triggered, this, &PlottingWidget::onRemoveLastMarker);
    connect(m_removeAllMarkersAction, &QAction::triggered, this, &PlottingWidget::onRemoveAllMarkers);
    connect(m_removeAllAnnotationsAction, &QAction::triggered, this, &PlottingWidget::onRemoveAllAnnotations);
    connect(m_zoomInAction, &QAction::triggered, this, &PlottingWidget::zoomIn);
    connect(m_zoomOutAction, &QAction::triggered, this, &PlottingWidget::zoomOut);
    connect(m_zoomFitAction, &QAction::triggered, this, &PlottingWidget::zoomToFit);
    connect(m_resetZoomAction, &QAction::triggered, this, &PlottingWidget::resetZoom);

    // 连接单独缩放功能
    connect(m_zoomXInAction, &QAction::triggered, this, &PlottingWidget::zoomXIn);
    connect(m_zoomXOutAction, &QAction::triggered, this, &PlottingWidget::zoomXOut);
    connect(m_zoomYInAction, &QAction::triggered, this, &PlottingWidget::zoomYIn);
    connect(m_zoomYOutAction, &QAction::triggered, this, &PlottingWidget::zoomYOut);
}

void PlottingWidget::onAddCurve()
{
    if (!m_hasTableData || m_tableData.columns.isEmpty()) {
        QMessageBox::warning(this, "添加曲线", "请先加载数据！");
        return;
    }

    showDataSelectionDialog("自定义曲线");
}

void PlottingWidget::onPressureProdDataPlot()
{
    if (!m_hasTableData || m_tableData.columns.isEmpty()) {
        QMessageBox::warning(this, "压力产量数据", "请先加载数据！");
        return;
    }

    showPressureProdDataDialog();
}

void PlottingWidget::onPressureDerivativePlot()
{
    if (!m_hasTableData || m_tableData.columns.isEmpty()) {
        QMessageBox::warning(this, "压力导数", "请先加载数据！");
        return;
    }

    showDataSelectionDialog("压力导数");
}

void PlottingWidget::showPressureProdDataDialog()
{
    QDialog dialog(this);
    applyDialogStyle(&dialog);  // 应用美化样式
    dialog.setWindowTitle("📊 压力产量数据分析");
    dialog.setModal(true);
    dialog.resize(680, 750);

    QVBoxLayout *layout = new QVBoxLayout(&dialog);

    QLabel *titleLabel = new QLabel("压力产量联合数据分析");
    titleLabel->setStyleSheet("font-size: 14pt; font-weight: bold; color: #333333; margin: 10px;");
    titleLabel->setAlignment(Qt::AlignCenter);
    layout->addWidget(titleLabel);

    // 创建Tab控件
    QTabWidget *tabWidget = new QTabWidget();

    // 压力数据选项卡
    QWidget *pressureTab = new QWidget();
    QGridLayout *pressureLayout = new QGridLayout(pressureTab);

    pressureLayout->addWidget(new QLabel("时间数据列:"), 0, 0);
    QComboBox *pressureTimeCombo = new QComboBox();
    for (int i = 0; i < m_tableData.headers.size(); ++i) {
        pressureTimeCombo->addItem(QString("%1 (列%2)").arg(m_tableData.headers[i]).arg(i + 1));
    }
    pressureLayout->addWidget(pressureTimeCombo, 0, 1);

    pressureLayout->addWidget(new QLabel("时间轴类型:"), 0, 2);
    QComboBox *pressureTimeAxisTypeCombo = new QComboBox();
    pressureTimeAxisTypeCombo->addItems(QStringList() << "常规坐标系" << "对数坐标系");
    pressureLayout->addWidget(pressureTimeAxisTypeCombo, 0, 3);

    pressureLayout->addWidget(new QLabel("压力数据列:"), 1, 0);
    QComboBox *pressureDataCombo = new QComboBox();
    for (int i = 0; i < m_tableData.headers.size(); ++i) {
        pressureDataCombo->addItem(QString("%1 (列%2)").arg(m_tableData.headers[i]).arg(i + 1));
    }
    if (pressureDataCombo->count() > 1) {
        pressureDataCombo->setCurrentIndex(1);
    }
    pressureLayout->addWidget(pressureDataCombo, 1, 1);

    pressureLayout->addWidget(new QLabel("压力轴类型:"), 1, 2);
    QComboBox *pressureAxisTypeCombo = new QComboBox();
    pressureAxisTypeCombo->addItems(QStringList() << "常规坐标系" << "对数坐标系");
    pressureLayout->addWidget(pressureAxisTypeCombo, 1, 3);

    pressureLayout->addWidget(new QLabel("曲线名称:"), 2, 0);
    QLineEdit *pressureNameEdit = new QLineEdit("压力数据");
    pressureLayout->addWidget(pressureNameEdit, 2, 1);

    pressureLayout->addWidget(new QLabel("压力轴单位:"), 2, 2);
    QLineEdit *pressureUnitEdit = new QLineEdit("MPa");
    pressureLayout->addWidget(pressureUnitEdit, 2, 3);

    pressureLayout->addWidget(new QLabel("线宽:"), 3, 0);
    QSpinBox *pressureLineWidthSpin = new QSpinBox();
    pressureLineWidthSpin->setRange(1, 10);
    pressureLineWidthSpin->setValue(2);
    pressureLayout->addWidget(pressureLineWidthSpin, 3, 1);

    pressureLayout->addWidget(new QLabel("点大小:"), 3, 2);
    QSpinBox *pressurePointSizeSpin = new QSpinBox();
    pressurePointSizeSpin->setRange(1, 20);
    pressurePointSizeSpin->setValue(4);
    pressureLayout->addWidget(pressurePointSizeSpin, 3, 3);

    pressureLayout->setRowStretch(4, 1);
    tabWidget->addTab(pressureTab, "压力数据设置");

    // 产量数据选项卡
    QWidget *productionTab = new QWidget();
    QGridLayout *productionLayout = new QGridLayout(productionTab);

    productionLayout->addWidget(new QLabel("时间数据列:"), 0, 0);
    QComboBox *productionTimeCombo = new QComboBox();
    for (int i = 0; i < m_tableData.headers.size(); ++i) {
        productionTimeCombo->addItem(QString("%1 (列%2)").arg(m_tableData.headers[i]).arg(i + 1));
    }
    productionLayout->addWidget(productionTimeCombo, 0, 1);

    productionLayout->addWidget(new QLabel("产量数据列:"), 1, 0);
    QComboBox *productionDataCombo = new QComboBox();
    for (int i = 0; i < m_tableData.headers.size(); ++i) {
        productionDataCombo->addItem(QString("%1 (列%2)").arg(m_tableData.headers[i]).arg(i + 1));
    }
    if (productionDataCombo->count() > 2) {
        productionDataCombo->setCurrentIndex(2);
    }
    productionLayout->addWidget(productionDataCombo, 1, 1);

    productionLayout->addWidget(new QLabel("产量轴类型:"), 1, 2);
    QComboBox *productionAxisTypeCombo = new QComboBox();
    productionAxisTypeCombo->addItems(QStringList() << "常规坐标系" << "对数坐标系");
    productionLayout->addWidget(productionAxisTypeCombo, 1, 3);

    productionLayout->addWidget(new QLabel("曲线名称:"), 2, 0);
    QLineEdit *productionNameEdit = new QLineEdit("产量数据");
    productionLayout->addWidget(productionNameEdit, 2, 1);

    productionLayout->addWidget(new QLabel("产量轴单位:"), 2, 2);
    QLineEdit *productionUnitEdit = new QLineEdit("m³/d");
    productionLayout->addWidget(productionUnitEdit, 2, 3);

    productionLayout->addWidget(new QLabel("曲线类型:"), 3, 0);
    QComboBox *curveTypeCombo = new QComboBox();
    curveTypeCombo->addItems(QStringList() << "时间vs产量" << "时间段vs产量");
    productionLayout->addWidget(curveTypeCombo, 3, 1);

    productionLayout->addWidget(new QLabel("线宽:"), 4, 0);
    QSpinBox *productionLineWidthSpin = new QSpinBox();
    productionLineWidthSpin->setRange(1, 10);
    productionLineWidthSpin->setValue(2);
    productionLayout->addWidget(productionLineWidthSpin, 4, 1);

    productionLayout->addWidget(new QLabel("点大小:"), 4, 2);
    QSpinBox *productionPointSizeSpin = new QSpinBox();
    productionPointSizeSpin->setRange(1, 20);
    productionPointSizeSpin->setValue(4);
    productionLayout->addWidget(productionPointSizeSpin, 4, 3);

    productionLayout->setRowStretch(5, 1);
    tabWidget->addTab(productionTab, "产量数据设置");

    // 公共设置选项卡
    QWidget *commonTab = new QWidget();
    QGridLayout *commonLayout = new QGridLayout(commonTab);

    commonLayout->addWidget(new QLabel("时间轴标签:"), 0, 0);
    QLineEdit *timeLabelEdit = new QLineEdit("时间");
    commonLayout->addWidget(timeLabelEdit, 0, 1);

    commonLayout->addWidget(new QLabel("时间轴单位:"), 0, 2);
    QLineEdit *timeUnitEdit = new QLineEdit("小时");
    commonLayout->addWidget(timeUnitEdit, 0, 3);

    QCheckBox *syncTimeAxisCheck = new QCheckBox("同步时间数据列");
    syncTimeAxisCheck->setChecked(true);
    commonLayout->addWidget(syncTimeAxisCheck, 1, 0, 1, 2);

    connect(syncTimeAxisCheck, &QCheckBox::toggled, [=](bool checked) {
        if (checked) {
            productionTimeCombo->setCurrentIndex(pressureTimeCombo->currentIndex());
        }
        productionTimeCombo->setEnabled(!checked);
    });

    connect(pressureTimeCombo, QOverload<int>::of(&QComboBox::currentIndexChanged), [=](int index) {
        if (syncTimeAxisCheck->isChecked()) {
            productionTimeCombo->setCurrentIndex(index);
        }
    });

    QCheckBox *newWindowCheck = new QCheckBox("在新窗口中显示");
    newWindowCheck->setChecked(true);
    commonLayout->addWidget(newWindowCheck, 2, 0, 1, 2);

    commonLayout->setRowStretch(3, 1);
    tabWidget->addTab(commonTab, "公共设置");

    layout->addWidget(tabWidget);

    QHBoxLayout *buttonLayout = new QHBoxLayout();
    QPushButton *okButton = new QPushButton("确定");
    QPushButton *cancelButton = new QPushButton("取消");
    buttonLayout->addStretch();
    buttonLayout->addWidget(okButton);
    buttonLayout->addWidget(cancelButton);
    layout->addLayout(buttonLayout);

    connect(okButton, &QPushButton::clicked, &dialog, &QDialog::accept);
    connect(cancelButton, &QPushButton::clicked, &dialog, &QDialog::reject);

    if (dialog.exec() == QDialog::Accepted) {
        // 获取压力数据设置
        int pressureTimeIndex = pressureTimeCombo->currentIndex();
        int pressureDataIndex = pressureDataCombo->currentIndex();
        AxisType pressureTimeAxisType = (pressureTimeAxisTypeCombo->currentIndex() == 0) ? AxisType::Linear : AxisType::Logarithmic;
        AxisType pressureAxisType = (pressureAxisTypeCombo->currentIndex() == 0) ? AxisType::Linear : AxisType::Logarithmic;

        // 获取产量数据设置
        int productionTimeIndex = productionTimeCombo->currentIndex();
        int productionDataIndex = productionDataCombo->currentIndex();
        AxisType productionAxisType = (productionAxisTypeCombo->currentIndex() == 0) ? AxisType::Linear : AxisType::Logarithmic;
        QString curveTypeStr = curveTypeCombo->currentText();

        // 数据验证
        if (pressureTimeIndex == pressureDataIndex || productionTimeIndex == productionDataIndex) {
            QMessageBox::warning(this, "创建图表", "时间数据列和数值数据列不能是同一列！");
            return;
        }

        // 创建压力曲线
        CurveData pressureCurve = createCurveFromTableData(
            pressureTimeIndex, pressureDataIndex, pressureNameEdit->text().trimmed(),
            QColor(255, 152, 0), pressureTimeAxisType, pressureAxisType,
            timeLabelEdit->text().trimmed(), "压力",
            timeUnitEdit->text().trimmed(), pressureUnitEdit->text().trimmed(),
            pressureLineWidthSpin->value(), pressurePointSizeSpin->value()
            );

        // 创建产量曲线
        CurveData productionCurve;
        if (curveTypeStr == "时间段vs产量") {
            productionCurve = createStepProductionCurve(
                productionTimeIndex, productionDataIndex, productionNameEdit->text().trimmed(),
                QColor(76, 175, 80), pressureTimeAxisType, productionAxisType,
                timeLabelEdit->text().trimmed(), "产量",
                timeUnitEdit->text().trimmed(), productionUnitEdit->text().trimmed(),
                productionLineWidthSpin->value(), productionPointSizeSpin->value()
                );
        } else {
            productionCurve = createProductionCurve(
                productionTimeIndex, productionDataIndex, productionNameEdit->text().trimmed(),
                QColor(76, 175, 80), pressureTimeAxisType, productionAxisType,
                timeLabelEdit->text().trimmed(), "产量",
                timeUnitEdit->text().trimmed(), productionUnitEdit->text().trimmed(),
                curveTypeStr, productionLineWidthSpin->value(), productionPointSizeSpin->value()
                );
        }

        // 创建双图窗口或添加到当前窗口
        if (newWindowCheck->isChecked()) {
            QString windowTitle = "压力产量联合分析";
            DualPlotWindow *dualWindow = createDualPlotWindow(windowTitle);
            dualWindow->addPressureCurve(pressureCurve);
            dualWindow->addProductionCurve(productionCurve);
            dualWindow->setAxisSettings(timeLabelEdit->text().trimmed() + " (" + timeUnitEdit->text().trimmed() + ")",
                                        "压力 (" + pressureUnitEdit->text().trimmed() + ")",
                                        "产量 (" + productionUnitEdit->text().trimmed() + ")");
            dualWindow->updatePlots();
        } else {
            // 添加到当前窗口（如果需要的话）
            addCurve(pressureCurve);
            addCurve(productionCurve);
            calculateDataBounds();
        }

        QMessageBox::information(this, "创建成功", "压力产量联合曲线已成功创建！");
    }
}

void PlottingWidget::showDataSelectionDialog(const QString &plotType)
{
    QDialog dialog(this);
    applyDialogStyle(&dialog);  // 应用美化样式
    dialog.setWindowTitle(QString("🎨 创建%1图").arg(plotType));
    dialog.setModal(true);
    dialog.resize(650, 550);

    QVBoxLayout *layout = new QVBoxLayout(&dialog);

    QHBoxLayout *nameLayout = new QHBoxLayout();
    nameLayout->addWidget(new QLabel("曲线名称:"));
    QLineEdit *nameEdit = new QLineEdit();
    nameEdit->setText(plotType);
    nameLayout->addWidget(nameEdit);
    layout->addLayout(nameLayout);

    QGroupBox *axisGroup = new QGroupBox("坐标轴设置");
    QGridLayout *axisLayout = new QGridLayout(axisGroup);

    axisLayout->addWidget(new QLabel("X轴数据:"), 0, 0);
    QComboBox *xCombo = new QComboBox();
    for (int i = 0; i < m_tableData.headers.size(); ++i) {
        xCombo->addItem(QString("%1 (列%2)").arg(m_tableData.headers[i]).arg(i + 1));
    }
    axisLayout->addWidget(xCombo, 0, 1);

    axisLayout->addWidget(new QLabel("X轴类型:"), 0, 2);
    QComboBox *xAxisTypeCombo = new QComboBox();
    xAxisTypeCombo->addItems(QStringList() << "常规坐标系" << "对数坐标系");
    axisLayout->addWidget(xAxisTypeCombo, 0, 3);

    axisLayout->addWidget(new QLabel("Y轴数据:"), 1, 0);
    QComboBox *yCombo = new QComboBox();
    for (int i = 0; i < m_tableData.headers.size(); ++i) {
        yCombo->addItem(QString("%1 (列%2)").arg(m_tableData.headers[i]).arg(i + 1));
    }
    if (yCombo->count() > 1) {
        yCombo->setCurrentIndex(1);
    }
    axisLayout->addWidget(yCombo, 1, 1);

    axisLayout->addWidget(new QLabel("Y轴类型:"), 1, 2);
    QComboBox *yAxisTypeCombo = new QComboBox();
    yAxisTypeCombo->addItems(QStringList() << "常规坐标系" << "对数坐标系");
    axisLayout->addWidget(yAxisTypeCombo, 1, 3);

    axisLayout->addWidget(new QLabel("X轴标签:"), 2, 0);
    QLineEdit *xLabelEdit = new QLineEdit();
    xLabelEdit->setPlaceholderText("如：时间");
    axisLayout->addWidget(xLabelEdit, 2, 1);

    axisLayout->addWidget(new QLabel("X轴单位:"), 2, 2);
    QLineEdit *xUnitEdit = new QLineEdit();
    xUnitEdit->setPlaceholderText("如：小时");
    axisLayout->addWidget(xUnitEdit, 2, 3);

    axisLayout->addWidget(new QLabel("Y轴标签:"), 3, 0);
    QLineEdit *yLabelEdit = new QLineEdit();
    yLabelEdit->setPlaceholderText("如：产量");
    axisLayout->addWidget(yLabelEdit, 3, 1);

    axisLayout->addWidget(new QLabel("Y轴单位:"), 3, 2);
    QLineEdit *yUnitEdit = new QLineEdit();
    yUnitEdit->setPlaceholderText("如：m³/d");
    axisLayout->addWidget(yUnitEdit, 3, 3);

    layout->addWidget(axisGroup);

    QGroupBox *styleGroup = new QGroupBox("曲线样式");
    QGridLayout *styleLayout = new QGridLayout(styleGroup);

    styleLayout->addWidget(new QLabel("颜色:"), 0, 0);
    QPushButton *colorButton = new QPushButton();
    QColor selectedColor = QColor::fromHsv(m_curves.size() * 45 % 360, 200, 200);
    colorButton->setStyleSheet(QString("background-color: %1; min-width: 60px; min-height: 25px;").arg(selectedColor.name()));
    connect(colorButton, &QPushButton::clicked, [&selectedColor, colorButton]() {
        QColor newColor = QColorDialog::getColor(selectedColor);
        if (newColor.isValid()) {
            selectedColor = newColor;
            colorButton->setStyleSheet(QString("background-color: %1; min-width: 60px; min-height: 25px;").arg(newColor.name()));
        }
    });
    styleLayout->addWidget(colorButton, 0, 1);

    styleLayout->addWidget(new QLabel("线宽:"), 0, 2);
    QSpinBox *lineWidthSpin = new QSpinBox();
    lineWidthSpin->setRange(1, 10);
    lineWidthSpin->setValue(2);
    styleLayout->addWidget(lineWidthSpin, 0, 3);

    styleLayout->addWidget(new QLabel("点大小:"), 1, 0);
    QSpinBox *pointSizeSpin = new QSpinBox();
    pointSizeSpin->setRange(1, 20);
    pointSizeSpin->setValue(4);
    styleLayout->addWidget(pointSizeSpin, 1, 1);

    layout->addWidget(styleGroup);

    QCheckBox *newWindowCheck = new QCheckBox("在新窗口中显示");
    newWindowCheck->setChecked(true);
    layout->addWidget(newWindowCheck);

    QHBoxLayout *buttonLayout = new QHBoxLayout();
    QPushButton *okButton = new QPushButton("确定");
    QPushButton *cancelButton = new QPushButton("取消");
    buttonLayout->addStretch();
    buttonLayout->addWidget(okButton);
    buttonLayout->addWidget(cancelButton);
    layout->addLayout(buttonLayout);

    connect(okButton, &QPushButton::clicked, &dialog, &QDialog::accept);
    connect(cancelButton, &QPushButton::clicked, &dialog, &QDialog::reject);

    if (dialog.exec() == QDialog::Accepted) {
        QString curveName = nameEdit->text().trimmed();
        if (curveName.isEmpty()) {
            curveName = plotType;
        }

        int xIndex = xCombo->currentIndex();
        int yIndex = yCombo->currentIndex();
        AxisType xAxisType = (xAxisTypeCombo->currentIndex() == 0) ? AxisType::Linear : AxisType::Logarithmic;
        AxisType yAxisType = (yAxisTypeCombo->currentIndex() == 0) ? AxisType::Linear : AxisType::Logarithmic;

        if (xIndex == yIndex) {
            QMessageBox::warning(this, "创建图表", "X轴和Y轴不能选择相同的数据列！");
            return;
        }

        if (xIndex < 0 || xIndex >= m_tableData.columns.size() ||
            yIndex < 0 || yIndex >= m_tableData.columns.size()) {
            QMessageBox::warning(this, "创建图表", "请选择有效的数据列！");
            return;
        }

        CurveData newCurve = createCurveFromTableData(xIndex, yIndex, curveName, selectedColor,
                                                      xAxisType, yAxisType,
                                                      xLabelEdit->text().trimmed(), yLabelEdit->text().trimmed(),
                                                      xUnitEdit->text().trimmed(), yUnitEdit->text().trimmed(),
                                                      lineWidthSpin->value(), pointSizeSpin->value());

        if (newWindowCheck->isChecked()) {
            QString windowTitle = QString("%1 - %2").arg(curveName).arg(plotType);
            PlotWindow *plotWindow = createPlotWindow(windowTitle, plotType);
            plotWindow->addCurve(newCurve);
        } else {
            addCurve(newCurve);
            calculateDataBounds();
        }

        QMessageBox::information(this, "创建成功", QString("%1 '%2' 已成功创建！").arg(plotType).arg(curveName));
    }
}

CurveData PlottingWidget::createProductionCurve(int timeIndex, int productionIndex, const QString &curveName,
                                                const QColor &color, AxisType timeAxisType, AxisType productionAxisType,
                                                const QString &timeLabel, const QString &productionLabel,
                                                const QString &timeUnit, const QString &productionUnit,
                                                const QString &, int lineWidth, int pointSize)
{
    CurveData curve;
    curve.name = curveName;
    curve.color = color;
    curve.visible = true;
    curve.xAxisType = timeAxisType;
    curve.yAxisType = productionAxisType;
    curve.lineWidth = lineWidth;
    curve.pointSize = pointSize;

    const QVector<double> &timeData = m_tableData.columns[timeIndex];
    const QVector<double> &productionData = m_tableData.columns[productionIndex];

    int dataSize = qMin(timeData.size(), productionData.size());

    for (int i = 0; i < dataSize; ++i) {
        if (isValidDataPoint(timeData[i], productionData[i])) {
            curve.xData.append(timeData[i]);
            curve.yData.append(productionData[i]);
        }
    }

    curve.xLabel = timeLabel.isEmpty() ? "时间" : timeLabel;
    curve.yLabel = productionLabel.isEmpty() ? "产量" : productionLabel;
    curve.xUnit = timeUnit;
    curve.yUnit = productionUnit;

    return curve;
}

CurveData PlottingWidget::createStepProductionCurve(int timeIndex, int productionIndex, const QString &curveName,
                                                    const QColor &color, AxisType timeAxisType, AxisType productionAxisType,
                                                    const QString &timeLabel, const QString &productionLabel,
                                                    const QString &timeUnit, const QString &productionUnit,
                                                    int lineWidth, int pointSize)
{
    CurveData curve;
    curve.name = curveName;
    curve.color = color;
    curve.visible = true;
    curve.xAxisType = timeAxisType;
    curve.yAxisType = productionAxisType;
    curve.lineWidth = lineWidth;
    curve.pointSize = pointSize;
    curve.drawType = CurveType::Step;

    const QVector<double> &timeData = m_tableData.columns[timeIndex];
    const QVector<double> &productionData = m_tableData.columns[productionIndex];

    int dataSize = qMin(timeData.size(), productionData.size());

    double currentTime = 0.0;

    for (int i = 0; i < dataSize; ++i) {
        if (!isValidDataPoint(timeData[i], productionData[i])) {
            continue;
        }

        double duration = timeData[i];
        double production = productionData[i];
        double nextTime = currentTime + duration;

        curve.xData.append(currentTime);
        curve.yData.append(production);
        curve.xData.append(nextTime);
        curve.yData.append(production);

        currentTime = nextTime;
    }

    curve.xLabel = timeLabel.isEmpty() ? "时间" : timeLabel;
    curve.yLabel = productionLabel.isEmpty() ? "产量" : productionLabel;
    curve.xUnit = timeUnit;
    curve.yUnit = productionUnit;

    return curve;
}

PlotWindow* PlottingWidget::createPlotWindow(const QString &title, const QString &dataType)
{
    PlotWindow *plotWindow = new PlotWindow(title, this);
    m_plotWindows.append(plotWindow);

    if (dataType == "产量数据" || dataType == "生产数据") {
        plotWindow->setAxisSettings(false, false, "时间 (小时)", "产量");
        plotWindow->setPlotTitle("产量数据分析");
    } else if (dataType == "压力数据") {
        plotWindow->setAxisSettings(false, false, "时间 (小时)", "压力 (MPa)");
        plotWindow->setPlotTitle("压力数据分析");
    } else if (dataType == "压力导数") {
        plotWindow->setAxisSettings(false, false, "时间 (小时)", "压力导数 (MPa)");
        plotWindow->setPlotTitle("压力导数分析");
    } else {
        plotWindow->setAxisSettings(false, false, "X轴", "Y轴");
        plotWindow->setPlotTitle("自定义曲线分析");
    }

    plotWindow->show();
    return plotWindow;
}

DualPlotWindow* PlottingWidget::createDualPlotWindow(const QString &title)
{
    DualPlotWindow *dualWindow = new DualPlotWindow(title, this);
    m_dualPlotWindows.append(dualWindow);
    dualWindow->show();
    return dualWindow;
}

CurveData PlottingWidget::createCurveFromTableData(int xColumn, int yColumn, const QString &curveName,
                                                   const QColor &color, AxisType xAxisType, AxisType yAxisType,
                                                   const QString &xLabel, const QString &yLabel,
                                                   const QString &xUnit, const QString &yUnit,
                                                   int lineWidth, int pointSize)
{
    CurveData curve;
    curve.name = curveName;
    curve.color = color;
    curve.xData = m_tableData.columns[xColumn];
    curve.yData = m_tableData.columns[yColumn];
    curve.xLabel = xLabel.isEmpty() ? m_tableData.headers[xColumn] : xLabel;
    curve.yLabel = yLabel.isEmpty() ? m_tableData.headers[yColumn] : yLabel;
    curve.xUnit = xUnit;
    curve.yUnit = yUnit;
    curve.visible = true;
    curve.xAxisType = xAxisType;
    curve.yAxisType = yAxisType;
    curve.lineWidth = lineWidth;
    curve.pointSize = pointSize;
    return curve;
}

void PlottingWidget::addCurve(const CurveData &curve)
{
    m_curves.append(curve);

    if (m_curves.size() == 1) {
        m_plotSettings.logScaleX = (curve.xAxisType == AxisType::Logarithmic);
        m_plotSettings.logScaleY = (curve.yAxisType == AxisType::Logarithmic);
        m_plotSettings.xAxisType = curve.xAxisType;
        m_plotSettings.yAxisType = curve.yAxisType;
    }

    updateCurvesList();
    updatePlot();
    emit curveAdded(curve.name);
}

void PlottingWidget::removeCurve(int index)
{
    if (index >= 0 && index < m_curves.size()) {
        QString name = m_curves[index].name;
        m_curves.removeAt(index);
        updateCurvesList();
        updatePlot();
        emit curveRemoved(name);
    }
}

// 绘图相关函数实现
void PlottingWidget::paintPlotArea(QPaintEvent *)
{
    QPainter painter(ui->widget_plot);
    painter.setRenderHint(QPainter::Antialiasing, true);
    painter.setRenderHint(QPainter::TextAntialiasing, true);

    QRect widgetRect = ui->widget_plot->rect();
    m_plotArea = QRect(80, 50, widgetRect.width() - 160, widgetRect.height() - 100);

    drawBackground(painter);

    if (m_plotSettings.showGrid) {
        drawGrid(painter);
    }

    drawAxes(painter);

    if (!m_curves.isEmpty()) {
        drawAllCurves(painter);
    } else {
        drawNoDataMessage(painter);
    }

    drawMarkers(painter);
    drawAnnotations(painter);

    if (m_isSelecting) {
        drawSelection(painter);
    }

    drawCoordinates(painter);

    if (m_plotSettings.showLegend && !m_curves.isEmpty()) {
        drawLegend(painter);
    }
}

void PlottingWidget::drawBackground(QPainter &painter)
{
    painter.fillRect(m_plotArea, m_plotSettings.backgroundColor);
    // 移除外围蓝色边框，只保留绘图区域的边框
    painter.setPen(QPen(Qt::black, 1));
    painter.drawRect(m_plotArea);
}

// 修改后的drawGrid函数 - 正确处理对数坐标系
void PlottingWidget::drawGrid(QPainter &painter)
{
    painter.setPen(QPen(m_plotSettings.gridColor, 1, Qt::DotLine));

    QVector<double> xLabels = generateOptimizedAxisLabels(m_plotSettings.xMin, m_plotSettings.xMax, m_plotSettings.xAxisType);
    QVector<double> yLabels = generateOptimizedAxisLabels(m_plotSettings.yMin, m_plotSettings.yMax, m_plotSettings.yAxisType);

    // 绘制X轴网格线 - 使用正确的对数转换
    for (double value : xLabels) {
        double x;
        if (m_plotSettings.xAxisType == AxisType::Logarithmic && m_plotSettings.xMin > 0) {
            double normalizedX = (log10(value) - log10(m_plotSettings.xMin)) /
                                 (log10(m_plotSettings.xMax) - log10(m_plotSettings.xMin));
            x = m_plotArea.left() + normalizedX * m_plotArea.width();
        } else {
            x = m_plotArea.left() + (value - m_plotSettings.xMin) /
                                        (m_plotSettings.xMax - m_plotSettings.xMin) * m_plotArea.width();
        }

        if (x >= m_plotArea.left() && x <= m_plotArea.right()) {
            painter.drawLine(x, m_plotArea.top(), x, m_plotArea.bottom());
        }
    }

    // 绘制Y轴网格线 - 使用正确的对数转换
    for (double value : yLabels) {
        double y;
        if (m_plotSettings.yAxisType == AxisType::Logarithmic && m_plotSettings.yMin > 0) {
            double normalizedY = (log10(value) - log10(m_plotSettings.yMin)) /
                                 (log10(m_plotSettings.yMax) - log10(m_plotSettings.yMin));
            y = m_plotArea.bottom() - normalizedY * m_plotArea.height();
        } else {
            y = m_plotArea.bottom() - (value - m_plotSettings.yMin) /
                                          (m_plotSettings.yMax - m_plotSettings.yMin) * m_plotArea.height();
        }

        if (y >= m_plotArea.top() && y <= m_plotArea.bottom()) {
            painter.drawLine(m_plotArea.left(), y, m_plotArea.right(), y);
        }
    }
}

// 格式化坐标轴标签，参考Saphir软件风格
QString PlottingWidget::formatAxisLabel(double value, bool isLog)
{
    if (isLog) {
        // 对数坐标：对于10的整数次幂显示简洁格式
        double logValue = log10(qAbs(value));
        if (qAbs(logValue - qRound(logValue)) < 0.01) {
            int power = qRound(logValue);
            if (power == 0) return "1";
            if (power == 1) return "10";
            if (power == 2) return "100";
            if (power == 3) return "1000";
            if (power == -1) return "0.1";
            if (power == -2) return "0.01";
            return QString("10^%1").arg(power);
        } else {
            // 非10的整数次幂，显示实际数值
            if (value >= 1000) {
                return QString::number(value, 'g', 2);
            } else if (value >= 1) {
                return QString::number(value, 'f', 0);
            } else {
                return QString::number(value, 'g', 2);
            }
        }
    } else {
        // 线性坐标：智能格式化
        if (qAbs(value) >= 100000) {
            return QString::number(value, 'e', 1);
        } else if (qAbs(value) >= 1000) {
            return QString::number(value, 'f', 0);
        } else if (qAbs(value) >= 1) {
            return QString::number(value, 'f', 1);
        } else if (qAbs(value) >= 0.01) {
            return QString::number(value, 'f', 2);
        } else if (value == 0) {
            return "0";
        } else {
            return QString::number(value, 'g', 2);
        }
    }
}

// 修改后的drawAxes函数 - 正确处理对数坐标系
void PlottingWidget::drawAxes(QPainter &painter)
{
    painter.setPen(QPen(Qt::black, 2));
    painter.setFont(QFont("Arial", 9));

    QVector<double> xLabels = generateOptimizedAxisLabels(m_plotSettings.xMin, m_plotSettings.xMax, m_plotSettings.xAxisType);
    QVector<double> yLabels = generateOptimizedAxisLabels(m_plotSettings.yMin, m_plotSettings.yMax, m_plotSettings.yAxisType);

    // 绘制X轴标签 - 参考Saphir软件风格
    for (double value : xLabels) {
        double x;
        if (m_plotSettings.xAxisType == AxisType::Logarithmic && m_plotSettings.xMin > 0) {
            double normalizedX = (log10(value) - log10(m_plotSettings.xMin)) /
                                 (log10(m_plotSettings.xMax) - log10(m_plotSettings.xMin));
            x = m_plotArea.left() + normalizedX * m_plotArea.width();
        } else {
            x = m_plotArea.left() + (value - m_plotSettings.xMin) /
                                        (m_plotSettings.xMax - m_plotSettings.xMin) * m_plotArea.width();
        }

        if (x >= m_plotArea.left() && x <= m_plotArea.right()) {
            painter.drawLine(x, m_plotArea.bottom(), x, m_plotArea.bottom() - 8);
            QString label = formatAxisLabel(value, m_plotSettings.xAxisType == AxisType::Logarithmic);
            QRect textRect(x - 30, m_plotArea.bottom() + 5, 60, 15);
            painter.drawText(textRect, Qt::AlignCenter, label);
        }
    }

    // 绘制Y轴标签 - 参考Saphir软件风格
    for (double value : yLabels) {
        double y;
        if (m_plotSettings.yAxisType == AxisType::Logarithmic && m_plotSettings.yMin > 0) {
            double normalizedY = (log10(value) - log10(m_plotSettings.yMin)) /
                                 (log10(m_plotSettings.yMax) - log10(m_plotSettings.yMin));
            y = m_plotArea.bottom() - normalizedY * m_plotArea.height();
        } else {
            y = m_plotArea.bottom() - (value - m_plotSettings.yMin) /
                                          (m_plotSettings.yMax - m_plotSettings.yMin) * m_plotArea.height();
        }

        if (y >= m_plotArea.top() && y <= m_plotArea.bottom()) {
            painter.drawLine(m_plotArea.left(), y, m_plotArea.left() + 8, y);
            QString label = formatAxisLabel(value, m_plotSettings.yAxisType == AxisType::Logarithmic);
            QRect textRect(m_plotArea.left() - 75, y - 8, 70, 16);
            painter.drawText(textRect, Qt::AlignRight | Qt::AlignVCenter, label);
        }
    }

    painter.setFont(QFont("Arial", 10, QFont::Bold));

    QRect xTitleRect(m_plotArea.left(), m_plotArea.bottom() + 30, m_plotArea.width(), 15);
    painter.drawText(xTitleRect, Qt::AlignCenter, m_plotSettings.xAxisTitle);

    painter.save();
    painter.translate(15, m_plotArea.center().y());
    painter.rotate(-90);
    painter.drawText(-80, -3, 160, 15, Qt::AlignCenter, m_plotSettings.yAxisTitle);
    painter.restore();

    painter.setFont(QFont("Arial", 12, QFont::Bold));
    QRect titleRect(m_plotArea.left(), 5, m_plotArea.width(), 35);
    painter.drawText(titleRect, Qt::AlignCenter, m_plotSettings.plotTitle);
}

void PlottingWidget::drawAllCurves(QPainter &painter)
{
    for (const CurveData &curve : m_curves) {
        if (curve.visible) {
            drawCurve(painter, curve);
        }
    }
}

void PlottingWidget::drawCurve(QPainter &painter, const CurveData &curve)
{
    if (curve.drawType == CurveType::Step) {
        drawStepCurve(painter, curve);
        return;
    }

    if (curve.xData.isEmpty() || curve.yData.isEmpty()) {
        return;
    }

    Qt::PenStyle penStyle = lineStyleToQt(curve.lineStyle);
    painter.setPen(QPen(curve.color, curve.lineWidth, penStyle));

    QVector<QPointF> points;
    int maxPoints = qMin(curve.xData.size(), curve.yData.size());

    for (int i = 0; i < maxPoints; ++i) {
        if (!isValidDataPoint(curve.xData[i], curve.yData[i])) {
            continue;
        }

        // 对数坐标需要正值
        if (m_plotSettings.xAxisType == AxisType::Logarithmic && curve.xData[i] <= 0) continue;
        if (m_plotSettings.yAxisType == AxisType::Logarithmic && curve.yData[i] <= 0) continue;

        QPointF dataPoint(curve.xData[i], curve.yData[i]);
        QPointF pixelPoint = dataToPixel(dataPoint);

        if (pixelPoint.x() >= m_plotArea.left() - 50 && pixelPoint.x() <= m_plotArea.right() + 50 &&
            pixelPoint.y() >= m_plotArea.top() - 50 && pixelPoint.y() <= m_plotArea.bottom() + 50) {
            points.append(pixelPoint);
        }
    }

    if (points.size() > 1) {
        painter.setClipRect(m_plotArea.adjusted(-5, -5, 5, 5));
        for (int i = 0; i < points.size() - 1; ++i) {
            painter.drawLine(points[i], points[i + 1]);
        }
        painter.setClipping(false);
    }

    painter.setBrush(curve.color);
    for (const QPointF &point : points) {
        if (m_plotArea.contains(point.toPoint())) {
            painter.drawEllipse(point, curve.pointSize/2, curve.pointSize/2);
        }
    }
}

void PlottingWidget::drawStepCurve(QPainter &painter, const CurveData &curve)
{
    if (curve.xData.isEmpty() || curve.yData.isEmpty()) {
        return;
    }

    Qt::PenStyle penStyle = lineStyleToQt(curve.lineStyle);
    painter.setPen(QPen(curve.color, curve.lineWidth, penStyle));

    painter.setClipRect(m_plotArea.adjusted(-5, -5, 5, 5));

    int dataSize = qMin(curve.xData.size(), curve.yData.size());

    // 绘制阶梯状曲线
    for (int i = 0; i < dataSize - 1; i += 2) {
        if (i + 1 >= dataSize) break;

        if (!isValidDataPoint(curve.xData[i], curve.yData[i]) ||
            !isValidDataPoint(curve.xData[i+1], curve.yData[i+1])) {
            continue;
        }

        // 对数坐标检查
        if (m_plotSettings.xAxisType == AxisType::Logarithmic &&
            (curve.xData[i] <= 0 || curve.xData[i+1] <= 0)) continue;
        if (m_plotSettings.yAxisType == AxisType::Logarithmic &&
            (curve.yData[i] <= 0 || curve.yData[i+1] <= 0)) continue;

        QPointF segmentStart = dataToPixel(QPointF(curve.xData[i], curve.yData[i]));
        QPointF segmentEnd = dataToPixel(QPointF(curve.xData[i+1], curve.yData[i+1]));

        painter.drawLine(segmentStart, segmentEnd);

        if (i + 2 < dataSize && i + 3 < dataSize) {
            if (isValidDataPoint(curve.xData[i+2], curve.yData[i+2])) {
                // 对数坐标检查
                if (m_plotSettings.xAxisType == AxisType::Logarithmic && curve.xData[i+2] <= 0) continue;
                if (m_plotSettings.yAxisType == AxisType::Logarithmic && curve.yData[i+2] <= 0) continue;

                QPointF nextSegmentStart = dataToPixel(QPointF(curve.xData[i+2], curve.yData[i+2]));
                QPointF verticalStart = segmentEnd;
                QPointF verticalEnd = QPointF(segmentEnd.x(), nextSegmentStart.y());
                painter.drawLine(verticalStart, verticalEnd);
            }
        }
    }

    painter.setClipping(false);

    painter.setBrush(curve.color);
    for (int i = 0; i < dataSize; i += 2) {
        if (i < dataSize && isValidDataPoint(curve.xData[i], curve.yData[i])) {
            // 对数坐标检查
            if (m_plotSettings.xAxisType == AxisType::Logarithmic && curve.xData[i] <= 0) continue;
            if (m_plotSettings.yAxisType == AxisType::Logarithmic && curve.yData[i] <= 0) continue;

            QPointF pixelPoint = dataToPixel(QPointF(curve.xData[i], curve.yData[i]));
            if (m_plotArea.contains(pixelPoint.toPoint())) {
                painter.drawEllipse(pixelPoint, curve.pointSize/2, curve.pointSize/2);
            }
        }
    }
}

void PlottingWidget::drawLegend(QPainter &painter)
{
    if (m_curves.isEmpty()) return;

    int visibleCurveCount = 0;
    for (const CurveData &curve : m_curves) {
        if (curve.visible) {
            visibleCurveCount++;
        }
    }

    if (visibleCurveCount == 0) return;

    painter.setFont(QFont("Arial", 9));
    QFontMetrics fm(painter.font());

    int lineHeight = fm.height() + 4;
    int legendWidth = 0;
    int legendHeight = visibleCurveCount * lineHeight + 20;

    for (const CurveData &curve : m_curves) {
        if (curve.visible) {
            int textWidth = fm.horizontalAdvance(curve.name) + 50;
            legendWidth = qMax(legendWidth, textWidth);
        }
    }
    legendWidth += 20;

    int legendX = m_plotArea.right() - legendWidth - 10 + m_legendOffset.x();
    int legendY = m_plotArea.top() + 10 + m_legendOffset.y();

    legendX = qMax(m_plotArea.left() + 10, qMin(legendX, m_plotArea.right() - legendWidth - 10));
    legendY = qMax(m_plotArea.top() + 10, qMin(legendY, m_plotArea.bottom() - legendHeight - 10));

    m_legendArea = QRect(legendX, legendY, legendWidth, legendHeight);

    // 绘制图例背景（使用渐变和阴影效果）
    painter.save();

    // 绘制阴影
    QRect shadowRect = m_legendArea.adjusted(3, 3, 3, 3);
    painter.setPen(Qt::NoPen);
    painter.setBrush(QColor(0, 0, 0, 50));
    painter.drawRoundedRect(shadowRect, 8, 8);

    // 绘制渐变背景
    QLinearGradient gradient(legendX, legendY, legendX, legendY + legendHeight);
    gradient.setColorAt(0, QColor(255, 255, 255, 245));
    gradient.setColorAt(1, QColor(227, 242, 253, 245));
    painter.setBrush(gradient);
    painter.setPen(QPen(QColor(33, 150, 243), 2));
    painter.drawRoundedRect(m_legendArea, 8, 8);

    painter.restore();

    // 绘制图例标题
    painter.setPen(QPen(QColor(25, 118, 210), 1));
    painter.setFont(QFont("Arial", 10, QFont::Bold));
    painter.drawText(legendX + 10, legendY + 18, "📋 图例");

    // 绘制分隔线
    painter.setPen(QPen(QColor(187, 222, 251), 2));
    painter.drawLine(legendX + 10, legendY + 25, legendX + legendWidth - 10, legendY + 25);

    // 绘制曲线列表
    painter.setFont(QFont("Arial", 8));
    int currentY = legendY + 32;

    for (const CurveData &curve : m_curves) {
        if (!curve.visible) continue;

        Qt::PenStyle penStyle = lineStyleToQt(curve.lineStyle);
        painter.setPen(QPen(curve.color, qMax(2, curve.lineWidth), penStyle));
        painter.drawLine(legendX + 10, currentY + lineHeight/2 - 2,
                         legendX + 35, currentY + lineHeight/2 - 2);

        painter.setBrush(curve.color);
        painter.setPen(QPen(curve.color, 1));
        painter.drawEllipse(legendX + 22 - 3, currentY + lineHeight/2 - 5, 6, 6);

        painter.setPen(QPen(QColor(33, 33, 33), 1));
        painter.drawText(legendX + 40, currentY + lineHeight/2 + 4, curve.name);

        currentY += lineHeight;
    }
}

void PlottingWidget::applyDialogStyle(QDialog *dialog)
{
    dialog->setStyleSheet(
        "QDialog { "
        "   background-color: #FAFAFA; "
        "}"
        "QLabel { "
        "   color: #424242; "
        "   font-weight: bold; "
        "}"
        "QPushButton { "
        "   background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1, "
        "                              stop: 0 #42A5F5, stop: 1 #2196F3); "
        "   border: none; "
        "   border-radius: 6px; "
        "   padding: 8px 16px; "
        "   font-weight: bold; "
        "   color: white; "
        "   min-width: 80px; "
        "   min-height: 32px; "
        "}"
        "QPushButton:hover { "
        "   background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1, "
        "                              stop: 0 #64B5F6, stop: 1 #42A5F5); "
        "}"
        "QGroupBox { "
        "   border: 2px solid #2196F3; "
        "   border-radius: 8px; "
        "   margin-top: 16px; "
        "   font-weight: bold; "
        "   color: #1976D2; "
        "   background-color: white; "
        "   padding-top: 10px; "
        "}"
        "QGroupBox::title { "
        "   subcontrol-origin: margin; "
        "   subcontrol-position: top center; "
        "   padding: 0 15px; "
        "   background-color: #2196F3; "
        "   color: white; "
        "   border-radius: 4px; "
        "}"
        );
}

void PlottingWidget::drawMarkers(QPainter &painter)
{
    painter.setPen(QPen(Qt::red, 2));
    painter.setBrush(Qt::red);

    for (const QPointF &marker : m_markers) {
        QPointF pixelPoint = dataToPixel(marker);
        painter.drawEllipse(pixelPoint, 6, 6);
        painter.drawLine(pixelPoint.x() - 10, pixelPoint.y(), pixelPoint.x() + 10, pixelPoint.y());
        painter.drawLine(pixelPoint.x(), pixelPoint.y() - 10, pixelPoint.x(), pixelPoint.y() + 10);
    }
}

void PlottingWidget::drawAnnotations(QPainter &painter)
{
    painter.setPen(QPen(Qt::darkBlue, 1));
    painter.setFont(QFont("Arial", 8));

    for (const auto &annotation : m_annotations) {
        QPointF pixelPoint = dataToPixel(annotation.first);
        QString text = annotation.second;

        QFontMetrics fm(painter.font());
        QRect textRect = fm.boundingRect(text);
        textRect.moveCenter(pixelPoint.toPoint());
        textRect.adjust(-3, -1, 3, 1);

        painter.fillRect(textRect, QColor(255, 255, 255, 200));
        painter.drawRect(textRect);
        painter.drawText(textRect, Qt::AlignCenter, text);
    }
}

void PlottingWidget::drawSelection(QPainter &painter)
{
    painter.setPen(QPen(Qt::blue, 1, Qt::DashLine));
    painter.setBrush(QColor(0, 0, 255, 30));
    painter.drawRect(m_selectionRect);
}

void PlottingWidget::drawNoDataMessage(QPainter &painter)
{
    painter.setPen(QPen(Qt::gray, 1));
    painter.setFont(QFont("Arial", 12));

    QString message = "暂无曲线数据\n请先在数据页面加载数据文件\n然后点击相应按钮添加曲线";

    QRect textRect = m_plotArea;
    painter.drawText(textRect, Qt::AlignCenter, message);
}

void PlottingWidget::drawCoordinates(QPainter &)
{
    if (!m_plotArea.contains(m_lastMousePos)) {
        m_coordinateLabel->hide();
        return;
    }

    QPointF dataPos = pixelToData(m_lastMousePos);
    QString coordText = QString("X: %1, Y: %2")
                            .arg(formatScientific(dataPos.x(), 3))
                            .arg(formatScientific(dataPos.y(), 3));

    m_coordinateLabel->setText(coordText);
    m_coordinateLabel->adjustSize();

    QPoint labelPos = mapToGlobal(m_lastMousePos);
    labelPos.setX(labelPos.x() + 15);
    labelPos.setY(labelPos.y() - m_coordinateLabel->height() - 5);

    m_coordinateLabel->move(mapFromGlobal(labelPos));
    m_coordinateLabel->show();
}

// 其他必要的函数实现...
// 包括所有的坐标转换、数据处理、UI更新等函数

void PlottingWidget::updatePlot()
{
    m_plotSettings.showGrid = ui->checkBox_showGrid->isChecked();
    m_plotSettings.showLegend = ui->checkBox_showLegend->isChecked();

    if (ui->widget_plot) {
        ui->widget_plot->update();
    }
}

void PlottingWidget::updateCurvesList()
{
    if (!m_curvesListWidget) return;

    m_curvesListWidget->clear();

    for (int i = 0; i < m_curves.size(); ++i) {
        const CurveData &curve = m_curves[i];
        QListWidgetItem *item = new QListWidgetItem();

        QString displayText = QString("[%1] %2 %3")
                                  .arg(i + 1, 2, 10, QChar('0'))
                                  .arg(curve.visible ? "●" : "○")
                                  .arg(curve.name);
        item->setText(displayText);
        item->setData(Qt::UserRole, i);

        QPixmap colorPixmap(20, 16);
        colorPixmap.fill(curve.color);
        QPainter painter(&colorPixmap);
        painter.setPen(QPen(Qt::black, 1));
        painter.drawRect(0, 0, 19, 15);
        item->setIcon(QIcon(colorPixmap));

        QFont font = item->font();
        if (!curve.visible) {
            font.setItalic(true);
            item->setForeground(QColor(128, 128, 128));
        } else {
            font.setBold(true);
            item->setForeground(QColor(0, 0, 0));
        }
        item->setFont(font);

        m_curvesListWidget->addItem(item);
    }

    if (m_curves.isEmpty()) {
        QListWidgetItem *item = new QListWidgetItem("暂无曲线数据");
        item->setForeground(QColor(128, 128, 128));
        item->setFlags(Qt::NoItemFlags);
        m_curvesListWidget->addItem(item);
    }
}

// 数据管理函数
void PlottingWidget::setTableData(const TableData &data)
{
    m_tableData = data;
    m_hasTableData = true;

    QString dataInfo = QString("数据信息：\n文件：%1\n行数：%2\n列数：%3")
                           .arg(data.fileName.isEmpty() ? "未命名" : data.fileName)
                           .arg(data.rowCount)
                           .arg(data.headers.size());
    ui->label_dataInfo->setText(dataInfo);
}

void PlottingWidget::setTableDataFromModel(QStandardItemModel* model, const QString &fileName)
{
    if (!model) {
        return;
    }

    TableData data;
    data.fileName = fileName;
    data.rowCount = model->rowCount();

    for (int col = 0; col < model->columnCount(); ++col) {
        QString header = model->headerData(col, Qt::Horizontal).toString();
        if (header.isEmpty()) {
            header = QString("列%1").arg(col + 1);
        }
        data.headers.append(header);
    }

    data.columns.resize(model->columnCount());
    for (int col = 0; col < model->columnCount(); ++col) {
        for (int row = 0; row < model->rowCount(); ++row) {
            QStandardItem* item = model->item(row, col);
            if (item) {
                bool ok;
                double value = item->text().toDouble(&ok);
                if (ok) {
                    data.columns[col].append(value);
                } else {
                    data.columns[col].append(0.0);
                }
            } else {
                data.columns[col].append(0.0);
            }
        }
    }

    setTableData(data);
}

// 修改后的坐标转换函数 - 正确处理对数坐标系
QPointF PlottingWidget::dataToPixel(const QPointF &dataPoint)
{
    double x, y;

    // X轴坐标转换
    if (m_plotSettings.xAxisType == AxisType::Logarithmic && dataPoint.x() > 0 && m_plotSettings.xMin > 0) {
        double normalizedX = (log10(dataPoint.x()) - log10(m_plotSettings.xMin)) /
                             (log10(m_plotSettings.xMax) - log10(m_plotSettings.xMin));
        x = m_plotArea.left() + normalizedX * m_plotArea.width();
    } else {
        if (m_plotSettings.xMax > m_plotSettings.xMin) {
            x = m_plotArea.left() + (dataPoint.x() - m_plotSettings.xMin) /
                                        (m_plotSettings.xMax - m_plotSettings.xMin) * m_plotArea.width();
        } else {
            x = m_plotArea.left();
        }
    }

    // Y轴坐标转换
    if (m_plotSettings.yAxisType == AxisType::Logarithmic && dataPoint.y() > 0 && m_plotSettings.yMin > 0) {
        double normalizedY = (log10(dataPoint.y()) - log10(m_plotSettings.yMin)) /
                             (log10(m_plotSettings.yMax) - log10(m_plotSettings.yMin));
        y = m_plotArea.bottom() - normalizedY * m_plotArea.height();
    } else {
        if (m_plotSettings.yMax > m_plotSettings.yMin) {
            y = m_plotArea.bottom() - (dataPoint.y() - m_plotSettings.yMin) /
                                          (m_plotSettings.yMax - m_plotSettings.yMin) * m_plotArea.height();
        } else {
            y = m_plotArea.bottom();
        }
    }

    return QPointF(x, y);
}

// 修改后的像素到数据转换函数 - 正确处理对数坐标系
QPointF PlottingWidget::pixelToData(const QPointF &pixelPoint)
{
    double x, y;

    // X轴反向转换
    if (m_plotSettings.xAxisType == AxisType::Logarithmic && m_plotSettings.xMin > 0 && m_plotSettings.xMax > 0) {
        double normalizedX = (pixelPoint.x() - m_plotArea.left()) / m_plotArea.width();
        x = pow(10, log10(m_plotSettings.xMin) + normalizedX * (log10(m_plotSettings.xMax) - log10(m_plotSettings.xMin)));
    } else {
        x = m_plotSettings.xMin + (pixelPoint.x() - m_plotArea.left()) / m_plotArea.width() *
                                      (m_plotSettings.xMax - m_plotSettings.xMin);
    }

    // Y轴反向转换
    if (m_plotSettings.yAxisType == AxisType::Logarithmic && m_plotSettings.yMin > 0 && m_plotSettings.yMax > 0) {
        double normalizedY = (m_plotArea.bottom() - pixelPoint.y()) / m_plotArea.height();
        y = pow(10, log10(m_plotSettings.yMin) + normalizedY * (log10(m_plotSettings.yMax) - log10(m_plotSettings.yMin)));
    } else {
        y = m_plotSettings.yMin + (m_plotArea.bottom() - pixelPoint.y()) / m_plotArea.height() *
                                      (m_plotSettings.yMax - m_plotSettings.yMin);
    }

    return QPointF(x, y);
}

// 修改后的数据边界计算 - 正确处理对数坐标系
void PlottingWidget::calculateDataBounds()
{
    if (m_curves.isEmpty()) return;

    double minX = 1e10, maxX = -1e10;
    double minY = 1e10, maxY = -1e10;
    bool hasValidData = false;

    for (const CurveData &curve : m_curves) {
        if (!curve.visible || curve.xData.isEmpty() || curve.yData.isEmpty()) {
            continue;
        }

        for (int i = 0; i < qMin(curve.xData.size(), curve.yData.size()); ++i) {
            double x = curve.xData[i];
            double y = curve.yData[i];

            if (!qIsFinite(x) || !qIsFinite(y)) continue;

            // 对数坐标系需要正值
            if (m_plotSettings.xAxisType == AxisType::Logarithmic && x <= 0) continue;
            if (m_plotSettings.yAxisType == AxisType::Logarithmic && y <= 0) continue;

            minX = qMin(minX, x);
            maxX = qMax(maxX, x);
            minY = qMin(minY, y);
            maxY = qMax(maxY, y);
            hasValidData = true;
        }
    }

    if (hasValidData && minX < maxX && minY < maxY) {
        QPair<double, double> xRange = calculateOptimalRange(minX, maxX, m_plotSettings.xAxisType == AxisType::Logarithmic);
        QPair<double, double> yRange = calculateOptimalRange(minY, maxY, m_plotSettings.yAxisType == AxisType::Logarithmic);

        m_plotSettings.xMin = xRange.first;
        m_plotSettings.xMax = xRange.second;
        m_plotSettings.yMin = yRange.first;
        m_plotSettings.yMax = yRange.second;
    }
}

// 修改后的最优范围计算函数 - 正确处理对数坐标系
QPair<double, double> PlottingWidget::calculateOptimalRange(double min, double max, bool isLog)
{
    if (max <= min) {
        return QPair<double, double>(min, max);
    }

    if (isLog) {
        // 对数坐标系：扩展到最近的10的整数次幂，参考Saphir软件
        if (min <= 0) min = 1e-10;  // 确保正值
        if (max <= 0) max = 1;

        double logMin = log10(min);
        double logMax = log10(max);

        // 向下取整到最近的10的次幂
        double rangeMin = pow(10, floor(logMin));

        // 向上取整到最近的10的次幂
        double rangeMax = pow(10, ceil(logMax));

        // 如果数据范围太小，至少显示一个数量级
        if (rangeMax / rangeMin < 10) {
            rangeMin = rangeMin / 10;
            rangeMax = rangeMax * 10;
        }

        return QPair<double, double>(rangeMin, rangeMax);
    } else {
        // 线性坐标系：智能计算合适的范围
        double range = max - min;

        if (range == 0) {
            // 如果所有数据相同，扩展范围
            return QPair<double, double>(min - 1, max + 1);
        }

        // 计算合适的刻度间隔
        double orderOfMagnitude = pow(10, floor(log10(range)));
        double normalizedRange = range / orderOfMagnitude;

        double tickInterval;
        if (normalizedRange <= 1.5) {
            tickInterval = orderOfMagnitude * 0.2;
        } else if (normalizedRange <= 3) {
            tickInterval = orderOfMagnitude * 0.5;
        } else if (normalizedRange <= 7) {
            tickInterval = orderOfMagnitude;
        } else {
            tickInterval = orderOfMagnitude * 2;
        }

        // 计算范围，确保包含所有数据点并留有余量
        double rangeMin = floor(min / tickInterval) * tickInterval;
        double rangeMax = ceil(max / tickInterval) * tickInterval;

        // 添加5%的边距，使数据点不会贴边
        double margin = (rangeMax - rangeMin) * 0.05;
        rangeMin -= margin;
        rangeMax += margin;

        // 如果范围包含0，确保0点对齐到刻度
        if (rangeMin < 0 && rangeMax > 0) {
            rangeMin = floor(rangeMin / tickInterval) * tickInterval;
            rangeMax = ceil(rangeMax / tickInterval) * tickInterval;
        }

        // 如果原始最小值接近0，将范围最小值设为0
        if (min > 0 && rangeMin < 0 && min < range * 0.2) {
            rangeMin = 0;
        }

        return QPair<double, double>(rangeMin, rangeMax);
    }
}

// 优化的对数坐标标签生成函数，参考Saphir软件
QVector<double> PlottingWidget::generateOptimizedAxisLabels(double min, double max, AxisType axisType)
{
    QVector<double> labels;

    if (max <= min) {
        return labels;
    }

    if (axisType == AxisType::Logarithmic) {
        // 对数坐标系：生成Saphir软件风格的对数标签
        if (min <= 0) min = 1e-10;
        if (max <= 0) max = 1;

        double logMin = log10(min);
        double logMax = log10(max);

        int startPower = qFloor(logMin);
        int endPower = qCeil(logMax);

        // 添加主要的10的整数次幂
        for (int power = startPower; power <= endPower; ++power) {
            double value = qPow(10, power);
            if (value >= min * 0.999 && value <= max * 1.001) {
                labels.append(value);
            }
        }

        // 根据范围决定是否显示中间刻度
        bool showMinorTicks = (endPower - startPower) <= 3;

        if (showMinorTicks) {
            // 添加次要刻度（2, 3, 5）
            for (int power = startPower; power < endPower; ++power) {
                for (int mult : {2, 3, 5}) {
                    double value = mult * qPow(10, power);
                    if (value > min && value < max) {
                        labels.append(value);
                    }
                }
            }
        }

        // 排序去重
        std::sort(labels.begin(), labels.end());
        labels.erase(std::unique(labels.begin(), labels.end(), [](double a, double b) {
                         return qAbs(a - b) < 1e-10;
                     }), labels.end());

    } else {
        // 线性坐标系：智能生成合适的刻度
        double range = max - min;
        double logRange = log10(range);
        int magnitude = qFloor(logRange);
        double step = qPow(10, magnitude);

        double labelCount = range / step;
        if (labelCount < 4) {
            step = step / 5;
        } else if (labelCount < 6) {
            step = step / 2;
        } else if (labelCount > 10) {
            step = step * 2;
        }

        double startValue = qFloor(min / step) * step;
        if (startValue < min) {
            startValue += step;
        }

        // 确保包含0（如果在范围内）
        if (min <= 0 && max >= 0) {
            labels.append(0);
        }

        for (double value = startValue; value <= max + step * 0.001; value += step) {
            if (value >= min && value <= max && qAbs(value) > 1e-10) {
                labels.append(value);
            }
        }

        std::sort(labels.begin(), labels.end());
        labels.erase(std::unique(labels.begin(), labels.end(), [](double a, double b) {
                         return qAbs(a - b) < 1e-10;
                     }), labels.end());
    }

    return labels;
}

// 鼠标事件处理
void PlottingWidget::mousePressEvent(QMouseEvent *event)
{
    QPoint plotPos = ui->widget_plot->mapFromParent(event->pos());

    if (!ui->widget_plot->rect().contains(plotPos)) {
        return;
    }

    m_lastMousePos = plotPos;

    if (m_legendArea.contains(plotPos)) {
        if (event->button() == Qt::LeftButton) {
            m_isDraggingLegend = true;
            m_legendDragStart = plotPos;
            return;
        }
    }

    if (!m_plotArea.contains(plotPos)) {
        return;
    }

    if (event->button() == Qt::LeftButton) {
        if (event->modifiers() & Qt::ControlModifier) {
            m_isSelecting = true;
            m_selectionStart = plotPos;
            m_selectionRect = QRect(m_selectionStart, m_selectionStart);
        } else {
            m_isDragging = true;
            m_isPanning = true;
        }
    }
}

void PlottingWidget::mouseMoveEvent(QMouseEvent *event)
{
    QPoint plotPos = ui->widget_plot->mapFromParent(event->pos());

    if (m_isDraggingLegend) {
        QPoint delta = plotPos - m_legendDragStart;
        m_legendOffset += delta;
        m_legendDragStart = plotPos;
        updatePlot();
        return;
    }

    m_lastMousePos = plotPos;

    if (m_isSelecting) {
        m_selectionRect = QRect(m_selectionStart, plotPos).normalized();
        updatePlot();
    } else if (m_isDragging && m_isPanning) {
        QPointF delta = plotPos - m_lastMousePos;
        panView(delta);
        m_lastMousePos = plotPos;
    }

    if (ui->widget_plot->rect().contains(plotPos) && m_plotArea.contains(plotPos)) {
        QPointF dataPos = pixelToData(plotPos);
        emit dataPointClicked(dataPos.x(), dataPos.y());
    }

    updatePlot();
}

void PlottingWidget::mouseReleaseEvent(QMouseEvent *event)
{
    Q_UNUSED(event)

    if (m_isDraggingLegend) {
        m_isDraggingLegend = false;
        return;
    }

    if (m_isSelecting && !m_selectionRect.isEmpty()) {
        QPointF topLeft = pixelToData(m_selectionRect.topLeft());
        QPointF bottomRight = pixelToData(m_selectionRect.bottomRight());

        m_plotSettings.xMin = qMin(topLeft.x(), bottomRight.x());
        m_plotSettings.xMax = qMax(topLeft.x(), bottomRight.x());
        m_plotSettings.yMin = qMin(topLeft.y(), bottomRight.y());
        m_plotSettings.yMax = qMax(topLeft.y(), bottomRight.y());
    }

    m_isDragging = false;
    m_isSelecting = false;
    m_isPanning = false;
    updatePlot();
}

void PlottingWidget::mouseDoubleClickEvent(QMouseEvent *event)
{
    QPoint plotPos = ui->widget_plot->mapFromParent(event->pos());

    if (m_plotArea.contains(plotPos)) {
        resetZoom();
    }
}

void PlottingWidget::wheelEvent(QWheelEvent *event)
{
    QPoint plotPos = ui->widget_plot->mapFromParent(event->position().toPoint());

    if (ui->widget_plot->rect().contains(plotPos) && m_plotArea.contains(plotPos)) {
        double factor = 1.0 + event->angleDelta().y() / 1200.0;
        zoomAtPoint(plotPos, factor);
    }
}

void PlottingWidget::contextMenuEvent(QContextMenuEvent *event)
{
    QPoint plotPos = ui->widget_plot->mapFromParent(event->pos());

    if (ui->widget_plot->rect().contains(plotPos) && m_plotArea.contains(plotPos)) {
        m_lastMousePos = plotPos;
        m_contextMenu->exec(event->globalPos());
    }
}

// 缩放和平移函数
void PlottingWidget::resetZoom()
{
    m_zoomFactor = 1.0;
    m_zoomFactorX = 1.0;
    m_zoomFactorY = 1.0;
    m_viewCenter = QPointF(0, 0);
    m_panOffset = QPointF(0, 0);
    calculateDataBounds();
    updatePlot();
}

void PlottingWidget::zoomIn()
{
    double factor = 1.25;
    QPointF center = m_plotArea.center();
    zoomAtPoint(center, factor);
}

void PlottingWidget::zoomOut()
{
    double factor = 0.8;
    QPointF center = m_plotArea.center();
    zoomAtPoint(center, factor);
}

void PlottingWidget::zoomToFit()
{
    calculateDataBounds();
    updatePlot();
}

// 新增单独缩放功能
void PlottingWidget::zoomXIn()
{
    if (m_plotSettings.xAxisType == AxisType::Logarithmic && m_plotSettings.xMin > 0 && m_plotSettings.xMax > 0) {
        // 对数坐标系下的X轴缩放
        double logMin = log10(m_plotSettings.xMin);
        double logMax = log10(m_plotSettings.xMax);
        double logRange = logMax - logMin;
        double newLogRange = logRange / 1.25;
        double logCenter = (logMin + logMax) / 2;
        m_plotSettings.xMin = pow(10, logCenter - newLogRange / 2);
        m_plotSettings.xMax = pow(10, logCenter + newLogRange / 2);
    } else {
        double xRange = m_plotSettings.xMax - m_plotSettings.xMin;
        double newXRange = xRange / 1.25;
        double centerX = (m_plotSettings.xMin + m_plotSettings.xMax) / 2;
        m_plotSettings.xMin = centerX - newXRange / 2;
        m_plotSettings.xMax = centerX + newXRange / 2;
    }
    updatePlot();
}

void PlottingWidget::zoomXOut()
{
    if (m_plotSettings.xAxisType == AxisType::Logarithmic && m_plotSettings.xMin > 0 && m_plotSettings.xMax > 0) {
        // 对数坐标系下的X轴缩放
        double logMin = log10(m_plotSettings.xMin);
        double logMax = log10(m_plotSettings.xMax);
        double logRange = logMax - logMin;
        double newLogRange = logRange * 1.25;
        double logCenter = (logMin + logMax) / 2;
        m_plotSettings.xMin = pow(10, logCenter - newLogRange / 2);
        m_plotSettings.xMax = pow(10, logCenter + newLogRange / 2);
    } else {
        double xRange = m_plotSettings.xMax - m_plotSettings.xMin;
        double newXRange = xRange * 1.25;
        double centerX = (m_plotSettings.xMin + m_plotSettings.xMax) / 2;
        m_plotSettings.xMin = centerX - newXRange / 2;
        m_plotSettings.xMax = centerX + newXRange / 2;
    }
    updatePlot();
}

void PlottingWidget::zoomYIn()
{
    if (m_plotSettings.yAxisType == AxisType::Logarithmic && m_plotSettings.yMin > 0 && m_plotSettings.yMax > 0) {
        // 对数坐标系下的Y轴缩放
        double logMin = log10(m_plotSettings.yMin);
        double logMax = log10(m_plotSettings.yMax);
        double logRange = logMax - logMin;
        double newLogRange = logRange / 1.25;
        double logCenter = (logMin + logMax) / 2;
        m_plotSettings.yMin = pow(10, logCenter - newLogRange / 2);
        m_plotSettings.yMax = pow(10, logCenter + newLogRange / 2);
    } else {
        double yRange = m_plotSettings.yMax - m_plotSettings.yMin;
        double newYRange = yRange / 1.25;
        double centerY = (m_plotSettings.yMin + m_plotSettings.yMax) / 2;
        m_plotSettings.yMin = centerY - newYRange / 2;
        m_plotSettings.yMax = centerY + newYRange / 2;
    }
    updatePlot();
}

void PlottingWidget::zoomYOut()
{
    if (m_plotSettings.yAxisType == AxisType::Logarithmic && m_plotSettings.yMin > 0 && m_plotSettings.yMax > 0) {
        // 对数坐标系下的Y轴缩放
        double logMin = log10(m_plotSettings.yMin);
        double logMax = log10(m_plotSettings.yMax);
        double logRange = logMax - logMin;
        double newLogRange = logRange * 1.25;
        double logCenter = (logMin + logMax) / 2;
        m_plotSettings.yMin = pow(10, logCenter - newLogRange / 2);
        m_plotSettings.yMax = pow(10, logCenter + newLogRange / 2);
    } else {
        double yRange = m_plotSettings.yMax - m_plotSettings.yMin;
        double newYRange = yRange * 1.25;
        double centerY = (m_plotSettings.yMin + m_plotSettings.yMax) / 2;
        m_plotSettings.yMin = centerY - newYRange / 2;
        m_plotSettings.yMax = centerY + newYRange / 2;
    }
    updatePlot();
}

void PlottingWidget::zoomAtPoint(const QPointF &point, double factor)
{
    QPointF dataPoint = pixelToData(point);

    if (m_plotSettings.xAxisType == AxisType::Logarithmic && m_plotSettings.xMin > 0 && m_plotSettings.xMax > 0 && dataPoint.x() > 0) {
        // 对数坐标系下的缩放
        double logMin = log10(m_plotSettings.xMin);
        double logMax = log10(m_plotSettings.xMax);
        double logPoint = log10(dataPoint.x());
        double logRange = logMax - logMin;
        double newLogRange = logRange / factor;

        m_plotSettings.xMin = pow(10, logPoint - newLogRange * (logPoint - logMin) / logRange);
        m_plotSettings.xMax = pow(10, logPoint + newLogRange * (logMax - logPoint) / logRange);
    } else {
        double xRange = m_plotSettings.xMax - m_plotSettings.xMin;
        double newXRange = xRange / factor;
        m_plotSettings.xMin = dataPoint.x() - newXRange * (dataPoint.x() - m_plotSettings.xMin) / xRange;
        m_plotSettings.xMax = dataPoint.x() + newXRange * (m_plotSettings.xMax - dataPoint.x()) / xRange;
    }

    if (m_plotSettings.yAxisType == AxisType::Logarithmic && m_plotSettings.yMin > 0 && m_plotSettings.yMax > 0 && dataPoint.y() > 0) {
        // 对数坐标系下的缩放
        double logMin = log10(m_plotSettings.yMin);
        double logMax = log10(m_plotSettings.yMax);
        double logPoint = log10(dataPoint.y());
        double logRange = logMax - logMin;
        double newLogRange = logRange / factor;

        m_plotSettings.yMin = pow(10, logPoint - newLogRange * (logPoint - logMin) / logRange);
        m_plotSettings.yMax = pow(10, logPoint + newLogRange * (logMax - logPoint) / logRange);
    } else {
        double yRange = m_plotSettings.yMax - m_plotSettings.yMin;
        double newYRange = yRange / factor;
        m_plotSettings.yMin = dataPoint.y() - newYRange * (dataPoint.y() - m_plotSettings.yMin) / yRange;
        m_plotSettings.yMax = dataPoint.y() + newYRange * (m_plotSettings.yMax - dataPoint.y()) / yRange;
    }

    updatePlot();
}

void PlottingWidget::panView(const QPointF &delta)
{
    if (m_plotSettings.xAxisType == AxisType::Logarithmic && m_plotSettings.xMin > 0 && m_plotSettings.xMax > 0) {
        // 对数坐标系下的平移
        double logMin = log10(m_plotSettings.xMin);
        double logMax = log10(m_plotSettings.xMax);
        double logRange = logMax - logMin;
        double logDelta = -delta.x() / m_plotArea.width() * logRange;
        m_plotSettings.xMin = pow(10, logMin + logDelta);
        m_plotSettings.xMax = pow(10, logMax + logDelta);
    } else {
        QPointF dataDelta = pixelToData(delta) - pixelToData(QPointF(0, 0));
        m_plotSettings.xMin -= dataDelta.x();
        m_plotSettings.xMax -= dataDelta.x();
    }

    if (m_plotSettings.yAxisType == AxisType::Logarithmic && m_plotSettings.yMin > 0 && m_plotSettings.yMax > 0) {
        // 对数坐标系下的平移
        double logMin = log10(m_plotSettings.yMin);
        double logMax = log10(m_plotSettings.yMax);
        double logRange = logMax - logMin;
        double logDelta = delta.y() / m_plotArea.height() * logRange;
        m_plotSettings.yMin = pow(10, logMin + logDelta);
        m_plotSettings.yMax = pow(10, logMax + logDelta);
    } else {
        QPointF dataDelta = pixelToData(delta) - pixelToData(QPointF(0, 0));
        m_plotSettings.yMin -= dataDelta.y();
        m_plotSettings.yMax -= dataDelta.y();
    }

    updatePlot();
}

// 槽函数实现
void PlottingWidget::onEditCurve()
{
    QMessageBox::information(this, "编辑曲线", "编辑曲线功能正在开发中！");
}

void PlottingWidget::onRemoveCurve()
{
    if (!m_curvesListWidget) return;

    QListWidgetItem *currentItem = m_curvesListWidget->currentItem();
    if (!currentItem) {
        QMessageBox::information(this, "删除曲线", "请先选择要删除的曲线！");
        return;
    }

    int index = currentItem->data(Qt::UserRole).toInt();
    if (index >= 0 && index < m_curves.size()) {
        QString curveName = m_curves[index].name;

        int ret = QMessageBox::question(this, "删除曲线",
                                        QString("确定要删除曲线 '%1' 吗？").arg(curveName),
                                        QMessageBox::Yes | QMessageBox::No);

        if (ret == QMessageBox::Yes) {
            removeCurve(index);
            QMessageBox::information(this, "删除曲线", QString("曲线 '%1' 已删除！").arg(curveName));
        }
    }
}

void PlottingWidget::onExportPlot()
{
    QString fileName = QFileDialog::getSaveFileName(this,
                                                    "导出图像",
                                                    QString("数据曲线_%1.png").arg(QDateTime::currentDateTime().toString("yyyyMMdd_hhmmss")),
                                                    "PNG图像 (*.png);;JPEG图像 (*.jpg)");

    if (!fileName.isEmpty()) {
        QPixmap pixmap(ui->widget_plot->size());
        pixmap.fill(Qt::white);
        ui->widget_plot->render(&pixmap);

        if (pixmap.save(fileName)) {
            QMessageBox::information(this, "导出成功", "图像已成功导出到: " + fileName);
            emit plotExported(fileName);
        } else {
            QMessageBox::warning(this, "导出失败", "无法保存文件: " + fileName);
        }
    }
}

void PlottingWidget::onCurveSelectionChanged()
{
    updatePlot();
}

void PlottingWidget::onMarkerAdded()
{
    QPointF dataPos = pixelToData(m_lastMousePos);
    m_markers.append(dataPos);
    updatePlot();

    QString message = QString("标记已添加在 (%1, %2)")
                          .arg(formatScientific(dataPos.x(), 3))
                          .arg(formatScientific(dataPos.y(), 3));
    QMessageBox::information(this, "添加标记", message);
}

void PlottingWidget::onAnnotationAdded()
{
    bool ok;
    QString text = QInputDialog::getText(this, "添加注释", "请输入注释文本:", QLineEdit::Normal, "", &ok);

    if (ok && !text.isEmpty()) {
        QPointF dataPos = pixelToData(m_lastMousePos);
        m_annotations.append(QPair<QPointF, QString>(dataPos, text));
        updatePlot();

        QString message = QString("注释已添加在 (%1, %2)")
                              .arg(formatScientific(dataPos.x(), 3))
                              .arg(formatScientific(dataPos.y(), 3));
        QMessageBox::information(this, "添加注释", message);
    }
}

void PlottingWidget::onRemoveAllMarkers()
{
    if (m_markers.isEmpty()) {
        QMessageBox::information(this, "删除标记", "当前没有标记！");
        return;
    }

    int ret = QMessageBox::question(this, "删除所有标记",
                                    QString("确定要删除所有 %1 个标记吗？").arg(m_markers.size()),
                                    QMessageBox::Yes | QMessageBox::No);

    if (ret == QMessageBox::Yes) {
        m_markers.clear();
        updatePlot();
        QMessageBox::information(this, "删除标记", "所有标记已删除！");
    }
}

void PlottingWidget::onRemoveLastMarker()
{
    if (m_markers.isEmpty()) {
        QMessageBox::information(this, "删除标记", "当前没有标记！");
        return;
    }

    m_markers.removeLast();
    updatePlot();
    QMessageBox::information(this, "删除标记", "最后一个标记已删除！");
}

void PlottingWidget::onRemoveAllAnnotations()
{
    if (m_annotations.isEmpty()) {
        QMessageBox::information(this, "删除注释", "当前没有注释！");
        return;
    }

    int ret = QMessageBox::question(this, "删除所有注释",
                                    QString("确定要删除所有 %1 个注释吗？").arg(m_annotations.size()),
                                    QMessageBox::Yes | QMessageBox::No);

    if (ret == QMessageBox::Yes) {
        m_annotations.clear();
        updatePlot();
        QMessageBox::information(this, "删除注释", "所有注释已删除！");
    }
}

void PlottingWidget::onColorSettingsChanged()
{
    QPushButton *button = qobject_cast<QPushButton*>(sender());
    if (!button) return;

    QColor currentColor = m_plotSettings.gridColor;
    QString title = "选择网格颜色";

    QColor newColor = QColorDialog::getColor(currentColor, this, title);

    if (newColor.isValid()) {
        button->setStyleSheet(QString("background-color: %1;").arg(newColor.name()));
        m_plotSettings.gridColor = newColor;
        updatePlot();
    }
}

void PlottingWidget::updateControlsFromSettings()
{
    ui->checkBox_showGrid->setChecked(m_plotSettings.showGrid);
    ui->checkBox_showLegend->setChecked(m_plotSettings.showLegend);
}

void PlottingWidget::paintEvent(QPaintEvent *event)
{
    QWidget::paintEvent(event);
}

// 辅助函数
QString PlottingWidget::formatScientific(double value, int decimals)
{
    if (qIsNaN(value) || !qIsFinite(value)) {
        return "N/A";
    }

    if (qAbs(value) >= 1000 || (qAbs(value) < 0.01 && value != 0)) {
        return QString::number(value, 'e', decimals);
    } else {
        return QString::number(value, 'f', decimals);
    }
}

bool PlottingWidget::isValidDataPoint(double x, double y)
{
    return !qIsNaN(x) && !qIsNaN(y) && qIsFinite(x) && qIsFinite(y);
}

double PlottingWidget::transformToLogScale(double value, bool isLog)
{
    if (isLog && value > 0) {
        return log10(value);
    }
    return value;
}

double PlottingWidget::transformFromLogScale(double value, bool isLog)
{
    if (isLog) {
        return pow(10, value);
    }
    return value;
}

// 其他数据管理函数
void PlottingWidget::setWellTestData(const WellTestData &data)
{
    m_currentData = data;
    m_wellTestDataSets.clear();
    m_wellTestDataSets.append(data);
    updatePlot();
}

void PlottingWidget::addWellTestData(const WellTestData &data)
{
    m_wellTestDataSets.append(data);
    updatePlot();
}

void PlottingWidget::clearAllData()
{
    m_wellTestDataSets.clear();
    m_currentData = WellTestData();
    m_markers.clear();
    m_annotations.clear();
    m_hasTableData = false;
    m_tableData = TableData();
    m_curves.clear();

    ui->label_dataInfo->setText("数据信息：未加载数据");
    updateCurvesList();
    updatePlot();
}

void PlottingWidget::removeDataSet(int index)
{
    if (index >= 0 && index < m_wellTestDataSets.size()) {
        m_wellTestDataSets.removeAt(index);
        updatePlot();
    }
}

void PlottingWidget::removeCurve(const QString &name)
{
    for (int i = 0; i < m_curves.size(); ++i) {
        if (m_curves[i].name == name) {
            removeCurve(i);
            break;
        }
    }
}

void PlottingWidget::updateCurve(int index, const CurveData &curve)
{
    if (index >= 0 && index < m_curves.size()) {
        m_curves[index] = curve;
        updatePlot();
    }
}

void PlottingWidget::setCurveVisible(int index, bool visible)
{
    if (index >= 0 && index < m_curves.size()) {
        m_curves[index].visible = visible;
        updatePlot();
    }
}

void PlottingWidget::setCurveVisible(const QString &name, bool visible)
{
    for (int i = 0; i < m_curves.size(); ++i) {
        if (m_curves[i].name == name) {
            setCurveVisible(i, visible);
            break;
        }
    }
}

QVector<CurveData> PlottingWidget::getAllCurves() const
{
    return m_curves;
}

int PlottingWidget::getCurveCount() const
{
    return m_curves.size();
}

void PlottingWidget::setPlotSettings(const PlotSettings &settings)
{
    m_plotSettings = settings;
    updatePlot();
}

PlotSettings PlottingWidget::getPlotSettings() const
{
    return m_plotSettings;
}

void PlottingWidget::resetToDefaultSettings()
{
    setupDefaultSettings();
    updatePlot();
}

void PlottingWidget::exportPlot(const QString &fileName, const QString &format)
{
    Q_UNUSED(format)

    QString ext = QFileInfo(fileName).suffix().toLower();

    if (ext == "pdf") {
        QPrinter printer(QPrinter::HighResolution);
        printer.setOutputFormat(QPrinter::PdfFormat);
        printer.setOutputFileName(fileName);
        printer.setPageSize(QPageSize::A4);
        printer.setPageOrientation(QPageLayout::Landscape);

        QPainter painter(&printer);
        ui->widget_plot->render(&painter);
    } else if (ext == "svg") {
        QSvgGenerator generator;
        generator.setFileName(fileName);
        generator.setSize(ui->widget_plot->size());
        generator.setViewBox(ui->widget_plot->rect());
        generator.setTitle("数据曲线分析");
        generator.setDescription("Data Curve Analysis");

        QPainter painter(&generator);
        ui->widget_plot->render(&painter);
    } else {
        QPixmap pixmap(ui->widget_plot->size());
        pixmap.fill(Qt::white);
        ui->widget_plot->render(&pixmap);

        if (pixmap.save(fileName)) {
            QMessageBox::information(this, "导出成功", "图像已成功导出到: " + fileName);
            emit plotExported(fileName);
        } else {
            QMessageBox::warning(this, "导出失败", "无法保存文件: " + fileName);
        }
    }
}

// 分析函数实现（占位符）
void PlottingWidget::performLogLogAnalysis()
{
    QMap<QString, double> results;
    QString message = "双对数分析完成\n使用当前显示的数据进行分析";

    QMessageBox::information(this, "分析完成", "双对数分析功能正在开发中！");
    emit analysisCompleted("双对数分析", results);
}

void PlottingWidget::performSemiLogAnalysis()
{
    QMap<QString, double> results;
    QMessageBox::information(this, "分析完成", "半对数分析功能正在开发中！");
    emit analysisCompleted("半对数分析", results);
}

void PlottingWidget::performCartesianAnalysis()
{
    QMap<QString, double> results;
    QMessageBox::information(this, "分析完成", "直角坐标分析功能正在开发中！");
    emit analysisCompleted("直角坐标分析", results);
}

void PlottingWidget::performDerivativeAnalysis()
{
    QMap<QString, double> results;
    QMessageBox::information(this, "分析完成", "压力导数分析功能正在开发中！");
    emit analysisCompleted("压力导数分析", results);
}

void PlottingWidget::performModelMatching()
{
    QMap<QString, double> results;
    QMessageBox::information(this, "分析完成", "模型匹配功能正在开发中！");
    emit analysisCompleted("模型匹配", results);
}
