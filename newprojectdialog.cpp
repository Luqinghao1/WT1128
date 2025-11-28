#include "newprojectdialog.h"
#include "ui_newprojectdialog.h"
#include <QDebug>
#include <QRegularExpressionValidator>
#include <QGraphicsOpacityEffect>
#include <QEasingCurve>

NewProjectDialog::NewProjectDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::NewProjectDialog),
    m_currentPageIndex(0),
    m_isModified(false),
    m_isApplyingPreset(false),
    m_isInitializing(true),
    m_fadeAnimation(nullptr),
    m_validationTimer(new QTimer(this))
{
    ui->setupUi(this);
    setWindowFlags(windowFlags() & ~Qt::WindowContextHelpButtonHint);
    init();
}

NewProjectDialog::~NewProjectDialog()
{
    delete ui;
}

void NewProjectDialog::init()
{
    // 设置当前日期和时间
    ui->dateTimeEdit->setDateTime(QDateTime::currentDateTime());
    ui->testDateEdit->setDate(QDate::currentDate());

    // 创建页面映射
    m_pageMap["项目基本信息"] = 0;
    m_pageMap["井筒参数"] = 1;
    m_pageMap["压裂参数"] = 2;
    m_pageMap["储层参数"] = 3;
    m_pageMap["流体特性"] = 4;
    m_pageMap["分析设置"] = 5;

    // 选中第一个导航项
    ui->navigationList->setCurrentRow(0);
    ui->stackedWidget->setCurrentIndex(0);

    // 初始化页面验证状态
    for (int i = 0; i < 6; ++i) {
        m_pageValidationStatus[i] = false;
    }

    // 设置验证定时器
    m_validationTimer->setSingleShot(true);
    m_validationTimer->setInterval(500);
    connect(m_validationTimer, &QTimer::timeout, this, [this]() {
        if (!m_isInitializing) {  // 只有在初始化完成后才进行验证
            validateCurrentPage();
            updateProgressBar();
        }
    });

    // 设置默认项目路径
    QString documentsPath = QStandardPaths::writableLocation(QStandardPaths::DocumentsLocation);
    QString defaultProjectPath = documentsPath + "/试井解释项目";
    ui->projectPathEdit->setText(defaultProjectPath);

    // 先填充默认值，再设置连接，避免初始化时触发验证
    populateDefaultValues();

    // 设置连接
    setupConnections();

    // 设置动画
    setupAnimations();

    // 更新控件状态
    updateControlsVisibility();

    // 更新状态标签
    updateStatusLabel();

    // 初始化完成，允许验证和数据变化响应
    m_isInitializing = false;

    // 初始化完成后进行一次无错误提示的验证
    validateCurrentPageSilently();
    updateProgressBar();
}

void NewProjectDialog::setupConnections()
{
    // 导航项点击
    connect(ui->navigationList, &QListWidget::itemClicked, this, &NewProjectDialog::onNavigationItemClicked);

    // 预设按钮
    connect(ui->presetShaleOilButton, &QPushButton::clicked, this, &NewProjectDialog::applyShaleOilPreset);
    connect(ui->presetTightOilButton, &QPushButton::clicked, this, &NewProjectDialog::applyTightOilPreset);
    connect(ui->presetConventionalButton, &QPushButton::clicked, this, &NewProjectDialog::applyConventionalOilPreset);

    // 项目路径选择
    connect(ui->selectPathButton, &QPushButton::clicked, this, &NewProjectDialog::onSelectProjectPath);

    // 井型选择
    connect(ui->wellTypeCombo, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &NewProjectDialog::onWellTypeChanged);

    // 储层类型选择
    connect(ui->reservoirTypeCombo, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &NewProjectDialog::onReservoirTypeChanged);

    // 压裂类型选择
    connect(ui->fractureTypeCombo, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &NewProjectDialog::onFractureTypeChanged);

    // 流体相选择
    connect(ui->oilCheckBox, &QCheckBox::toggled, this, &NewProjectDialog::onFluidPhaseChanged);
    connect(ui->gasCheckBox, &QCheckBox::toggled, this, &NewProjectDialog::onFluidPhaseChanged);
    connect(ui->waterCheckBox, &QCheckBox::toggled, this, &NewProjectDialog::onFluidPhaseChanged);

    // 分析类型变化
    connect(ui->standardRadioButton, &QRadioButton::toggled, this, &NewProjectDialog::onAnalysisTypeChanged);
    connect(ui->nonlinearRadioButton, &QRadioButton::toggled, this, &NewProjectDialog::onAnalysisTypeChanged);
    connect(ui->multiLayerRadioButton, &QRadioButton::toggled, this, &NewProjectDialog::onAnalysisTypeChanged);
    connect(ui->fractureRadioButton, &QRadioButton::toggled, this, &NewProjectDialog::onAnalysisTypeChanged);

    // 总压缩系数计算
    connect(ui->compressibilitySpin, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &NewProjectDialog::calculateTotalCompressibility);
    connect(ui->compressibilityFluidSpin, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &NewProjectDialog::calculateTotalCompressibility);

    // 压裂参数计算
    connect(ui->fractureStagesSpin, QOverload<int>::of(&QSpinBox::valueChanged), this, &NewProjectDialog::calculateFractureParameters);
    connect(ui->horizontalLengthSpin, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &NewProjectDialog::calculateFractureParameters);

    // 导航按钮
    connect(ui->nextButton, &QPushButton::clicked, this, &NewProjectDialog::onNextButtonClicked);
    connect(ui->backButton, &QPushButton::clicked, this, &NewProjectDialog::onBackButtonClicked);
    connect(ui->cancelButton, &QPushButton::clicked, this, &NewProjectDialog::onCancelButtonClicked);
    connect(ui->helpButton, &QPushButton::clicked, this, &NewProjectDialog::onHelpButtonClicked);

    // 连接所有输入控件的变化信号到数据变化槽
    QList<QLineEdit*> lineEdits = this->findChildren<QLineEdit*>();
    for (QLineEdit* lineEdit : lineEdits) {
        connect(lineEdit, &QLineEdit::textChanged, this, &NewProjectDialog::onDataChanged);
    }

    QList<QComboBox*> comboBoxes = this->findChildren<QComboBox*>();
    for (QComboBox* comboBox : comboBoxes) {
        connect(comboBox, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &NewProjectDialog::onDataChanged);
    }

    QList<QDoubleSpinBox*> doubleSpinBoxes = this->findChildren<QDoubleSpinBox*>();
    for (QDoubleSpinBox* spinBox : doubleSpinBoxes) {
        connect(spinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged), this, &NewProjectDialog::onDataChanged);
    }

    QList<QSpinBox*> spinBoxes = this->findChildren<QSpinBox*>();
    for (QSpinBox* spinBox : spinBoxes) {
        connect(spinBox, QOverload<int>::of(&QSpinBox::valueChanged), this, &NewProjectDialog::onDataChanged);
    }
}

void NewProjectDialog::setupAnimations()
{
    // 设置页面切换动画
    m_fadeAnimation = new QPropertyAnimation(this);
    m_fadeAnimation->setDuration(300);
    m_fadeAnimation->setEasingCurve(QEasingCurve::InOutQuad);
    connect(m_fadeAnimation, &QPropertyAnimation::finished, this, &NewProjectDialog::onAnimationFinished);
}

void NewProjectDialog::populateDefaultValues()
{
    // 设置项目信息的默认值
    ui->projectNameEdit->setText(tr("页岩油多段压裂水平井试井项目"));
    ui->wellNameEdit->setText("Demo-001");  // 给一个默认井名
    ui->fieldNameEdit->setText("示例油田");  // 给一个默认油田名
    ui->engineerEdit->setText("工程师");    // 给一个默认工程师名

    // 设置井型和储层类型默认值
    ui->wellTypeCombo->setCurrentText(tr("水平井"));
    ui->reservoirTypeCombo->setCurrentText(tr("页岩油藏"));
    ui->fractureTypeCombo->setCurrentText(tr("多段压裂"));

    // 设置其他默认值
    ui->testTypeCombo->setCurrentIndex(0);
    ui->analysisTypeCombo->setCurrentIndex(0);
    ui->testDurationSpin->setValue(72.0);

    // 井筒参数默认值（统一使用国际标准单位制）
    ui->wellRadiusSpin->setValue(0.0912);         // m (约0.3ft)
    ui->horizontalLengthSpin->setValue(457.2);    // m (约1500ft)
    ui->perforationTopSpin->setValue(2438.4);     // m (约8000ft)
    ui->perforationBottomSpin->setValue(2468.9);  // m (约8100ft)

    // 压裂参数默认值
    ui->fractureStagesSpin->setValue(15);
    ui->fractureSpacingSpin->setValue(91.4);      // m (约300ft)
    ui->fractureHalfLengthSpin->setValue(61.0);   // m (约200ft)
    ui->fractureConductivitySpin->setValue(15.2); // mD·m (约50 md·ft)

    // 井下设备默认值
    ui->tuningSpin->setValue(114.3);              // mm (约4.5inch)
    ui->casingSpin->setValue(177.8);              // mm (约7inch)
    ui->completionCombo->setCurrentIndex(2);
    ui->bhpGaugeSpin->setValue(2453.6);           // m (约8050ft)

    // 储层参数默认值
    ui->reservoirModelCombo->setCurrentText(tr("双重介质模型"));
    ui->reservoirPressureSpin->setValue(34.5);    // MPa (约5000psi)
    ui->porositySpin->setValue(0.08);
    ui->matrixPorositySpin->setValue(0.06);
    ui->fracturePorositySpin->setValue(0.02);
    ui->reservoirTempSpin->setValue(93.3);        // °C (约200°F)
    ui->permeabilitySpin->setValue(0.001);        // mD
    ui->matrixPermeabilitySpin->setValue(0.0001);
    ui->fracturePermeabilitySpin->setValue(0.1);
    ui->compressibilitySpin->setValue(0.000145);  // 1/MPa (约0.000001/psi)

    // 边界条件默认值
    ui->boundaryTypeCombo->setCurrentIndex(2);
    ui->distanceToBoundarySpin->setValue(609.6);  // m (约2000ft)
    ui->anisotropySpin->setValue(0.1);
    ui->formationDipSpin->setValue(5.0);

    // 流体属性默认值
    ui->fluidSystemCombo->setCurrentIndex(0);
    ui->referencePressureSpin->setValue(34.5);    // MPa (约5000psi)
    ui->oilCheckBox->setChecked(true);
    ui->gasCheckBox->setChecked(true);
    ui->waterCheckBox->setChecked(false);
    ui->viscositySpin->setValue(0.8);             // mPa·s
    ui->fvfSpin->setValue(1.3);                   // m³/m³
    ui->compressibilityFluidSpin->setValue(0.00290); // 1/MPa (约0.00002/psi)

    // 分析设置默认值
    ui->fractureRadioButton->setChecked(true);
    ui->analysisModelCombo->setCurrentText(tr("压裂井模型"));
    ui->flowRegimeCombo->setCurrentText(tr("双线性流"));
    ui->methodCombo->setCurrentIndex(0);
    ui->optimizationMethodCombo->setCurrentIndex(0);
    ui->pressureDifferenceDefinitionCombo->setCurrentIndex(0);
    ui->timeFormatCombo->setCurrentIndex(0);
    ui->primaryPlotCombo->setCurrentIndex(0);
    ui->showGridCheckBox->setChecked(true);
    ui->adjustScaleCheckBox->setChecked(true);
    ui->generateReportCheckBox->setChecked(true);
    ui->saveParametersCheckBox->setChecked(true);

    // 高级设置默认值
    ui->enableTransientFlowCheckBox->setChecked(true);
    ui->enableNonDarcyFlowCheckBox->setChecked(false);

    // 计算相关参数
    calculateTotalCompressibility();
    calculateFractureParameters();
}

void NewProjectDialog::applyShaleOilPreset()
{
    m_isApplyingPreset = true;

    ui->wellTypeCombo->setCurrentText(tr("水平井"));
    ui->reservoirTypeCombo->setCurrentText(tr("页岩油藏"));
    ui->fractureTypeCombo->setCurrentText(tr("多段压裂"));

    ui->reservoirModelCombo->setCurrentText(tr("双重介质模型"));
    ui->porositySpin->setValue(0.08);
    ui->matrixPorositySpin->setValue(0.06);
    ui->fracturePorositySpin->setValue(0.02);
    ui->permeabilitySpin->setValue(0.001);
    ui->matrixPermeabilitySpin->setValue(0.0001);
    ui->fracturePermeabilitySpin->setValue(0.1);
    ui->boundaryTypeCombo->setCurrentText(tr("封闭边界"));
    ui->anisotropySpin->setValue(0.1);

    ui->fractureStagesSpin->setValue(15);
    ui->fractureSpacingSpin->setValue(300.0);
    ui->fractureHalfLengthSpin->setValue(200.0);
    ui->fractureConductivitySpin->setValue(50.0);

    ui->oilCheckBox->setChecked(true);
    ui->gasCheckBox->setChecked(true);
    ui->waterCheckBox->setChecked(false);
    ui->viscositySpin->setValue(0.8);
    ui->fvfSpin->setValue(1.3);

    ui->fractureRadioButton->setChecked(true);
    ui->analysisModelCombo->setCurrentText(tr("压裂井模型"));
    ui->flowRegimeCombo->setCurrentText(tr("双线性流"));

    calculateTotalCompressibility();
    calculateFractureParameters();

    QMessageBox::information(this, tr("预设配置"), tr("已应用页岩油多段压裂水平井预设配置"));
    m_isModified = true;
    m_isApplyingPreset = false;
}

void NewProjectDialog::applyTightOilPreset()
{
    m_isApplyingPreset = true;

    ui->wellTypeCombo->setCurrentText(tr("水平井"));
    ui->reservoirTypeCombo->setCurrentText(tr("致密油藏"));
    ui->fractureTypeCombo->setCurrentText(tr("多段压裂"));

    ui->reservoirModelCombo->setCurrentText(tr("双重介质模型"));
    ui->porositySpin->setValue(0.10);
    ui->permeabilitySpin->setValue(0.01);
    ui->boundaryTypeCombo->setCurrentText(tr("封闭边界"));
    ui->fractureStagesSpin->setValue(12);
    ui->fractureSpacingSpin->setValue(400.0);

    ui->oilCheckBox->setChecked(true);
    ui->gasCheckBox->setChecked(false);
    ui->waterCheckBox->setChecked(false);
    ui->viscositySpin->setValue(1.2);
    ui->fvfSpin->setValue(1.25);

    ui->fractureRadioButton->setChecked(true);
    ui->analysisModelCombo->setCurrentText(tr("压裂井模型"));
    ui->flowRegimeCombo->setCurrentText(tr("线性流"));

    calculateTotalCompressibility();
    calculateFractureParameters();

    QMessageBox::information(this, tr("预设配置"), tr("已应用致密油预设配置"));
    m_isModified = true;
    m_isApplyingPreset = false;
}

void NewProjectDialog::applyConventionalOilPreset()
{
    m_isApplyingPreset = true;

    ui->wellTypeCombo->setCurrentText(tr("垂直井"));
    ui->reservoirTypeCombo->setCurrentText(tr("常规油藏"));
    ui->fractureTypeCombo->setCurrentText(tr("未压裂"));

    ui->reservoirModelCombo->setCurrentText(tr("均质模型"));
    ui->porositySpin->setValue(0.18);
    ui->permeabilitySpin->setValue(50.0);
    ui->boundaryTypeCombo->setCurrentText(tr("无限大储层"));

    ui->oilCheckBox->setChecked(true);
    ui->gasCheckBox->setChecked(false);
    ui->waterCheckBox->setChecked(false);
    ui->viscositySpin->setValue(2.0);
    ui->fvfSpin->setValue(1.2);

    ui->standardRadioButton->setChecked(true);
    ui->analysisModelCombo->setCurrentText(tr("经典均质模型"));
    ui->flowRegimeCombo->setCurrentText(tr("径向流"));

    calculateTotalCompressibility();

    QMessageBox::information(this, tr("预设配置"), tr("已应用常规油藏预设配置"));
    m_isModified = true;
    m_isApplyingPreset = false;
}

void NewProjectDialog::onWellTypeChanged()
{
    if (m_isApplyingPreset || m_isInitializing) return;

    QString wellType = ui->wellTypeCombo->currentText();
    bool isHorizontal = (wellType == tr("水平井"));

    ui->horizontalLengthLabel->setVisible(isHorizontal);
    ui->horizontalLengthSpin->setVisible(isHorizontal);

    if (isHorizontal) {
        ui->fractureTypeCombo->setCurrentText(tr("多段压裂"));
    } else {
        ui->fractureTypeCombo->setCurrentText(tr("未压裂"));
    }

    updateControlsVisibility();
}

void NewProjectDialog::onReservoirTypeChanged()
{
    if (m_isApplyingPreset || m_isInitializing) return;

    QString reservoirType = ui->reservoirTypeCombo->currentText();

    if (reservoirType.contains(tr("页岩"))) {
        ui->reservoirModelCombo->setCurrentText(tr("分形介质模型"));
        ui->boundaryTypeCombo->setCurrentText(tr("封闭边界"));
        ui->anisotropySpin->setValue(0.1);
    } else if (reservoirType.contains(tr("致密"))) {
        ui->reservoirModelCombo->setCurrentText(tr("双重介质模型"));
        ui->boundaryTypeCombo->setCurrentText(tr("封闭边界"));
        ui->anisotropySpin->setValue(0.5);
    } else {
        ui->reservoirModelCombo->setCurrentText(tr("均质模型"));
        ui->boundaryTypeCombo->setCurrentText(tr("无限大储层"));
        ui->anisotropySpin->setValue(1.0);
    }
}

void NewProjectDialog::onFractureTypeChanged()
{
    if (m_isApplyingPreset || m_isInitializing) return;

    QString fractureType = ui->fractureTypeCombo->currentText();
    bool isFractured = !fractureType.contains(tr("未压裂"));

    ui->fractureParamsGroup->setVisible(isFractured);

    if (isFractured) {
        ui->fractureRadioButton->setChecked(true);
        if (fractureType.contains(tr("多段"))) {
            ui->fractureStagesSpin->setValue(15);
            ui->flowRegimeCombo->setCurrentText(tr("双线性流"));
        } else {
            ui->fractureStagesSpin->setValue(1);
            ui->flowRegimeCombo->setCurrentText(tr("线性流"));
        }
    } else {
        ui->standardRadioButton->setChecked(true);
        ui->flowRegimeCombo->setCurrentText(tr("径向流"));
    }

    updateControlsVisibility();
}

void NewProjectDialog::calculateTotalCompressibility()
{
    if (m_isApplyingPreset || m_isInitializing) return;

    double rockCompressibility = ui->compressibilitySpin->value();
    double fluidCompressibility = ui->compressibilityFluidSpin->value();
    double totalCompressibility = rockCompressibility + fluidCompressibility;

    ui->totalCompressibilitySpin->setValue(totalCompressibility);
}

void NewProjectDialog::calculateFractureParameters()
{
    if (m_isApplyingPreset || m_isInitializing) return;

    double horizontalLength = ui->horizontalLengthSpin->value();
    int stages = ui->fractureStagesSpin->value();

    if (stages > 1) {
        double spacing = horizontalLength / (stages - 1);
        ui->fractureSpacingSpin->setValue(spacing);
    }

    // 更新计算结果显示 - 使用国际标准单位制
    double totalVolume = horizontalLength * ui->fractureHalfLengthSpin->value() * 2;
    double effectiveVolume = totalVolume * ui->porositySpin->value();
    double totalFractureLength = stages * ui->fractureHalfLengthSpin->value() * 2;

    ui->totalFractureVolumeValue->setText(QString("%1 m³").arg(totalVolume, 0, 'e', 2));
    ui->effectiveReservoirVolumeValue->setText(QString("%1 m³").arg(effectiveVolume, 0, 'e', 2));
    ui->averageSpacingValue->setText(QString("%1 m").arg(ui->fractureSpacingSpin->value(), 0, 'f', 1));
    ui->totalFractureLengthValue->setText(QString("%1 m").arg(totalFractureLength, 0, 'f', 0));
}

void NewProjectDialog::onNavigationItemClicked(QListWidgetItem *item)
{
    QString pageName = item->text();
    if (m_pageMap.contains(pageName)) {
        if (m_currentPageIndex >= 0 && !validateCurrentPage()) {
            ui->navigationList->setCurrentRow(m_currentPageIndex);
            return;
        }

        int pageIndex = m_pageMap[pageName];
        animatePageTransition(m_currentPageIndex, pageIndex);
        m_currentPageIndex = pageIndex;

        ui->backButton->setEnabled(pageIndex > 0);

        if (pageIndex == ui->stackedWidget->count() - 1) {
            ui->nextButton->setText(tr("完成创建"));
            ui->nextButton->setIcon(QIcon(":/new/prefix1/Resource/complete_icon.png"));
        } else {
            ui->nextButton->setText(tr("下一步 >>"));
            ui->nextButton->setIcon(QIcon(":/new/prefix1/Resource/next_icon.png"));
        }

        updateStatusLabel();
        updateProgressBar();
    }
}

void NewProjectDialog::animatePageTransition(int fromPage, int toPage)
{
    if (m_fadeAnimation->state() == QPropertyAnimation::Running) {
        m_fadeAnimation->stop();
    }

    QWidget* currentWidget = ui->stackedWidget->widget(fromPage);
    QGraphicsOpacityEffect* effect = new QGraphicsOpacityEffect(currentWidget);
    currentWidget->setGraphicsEffect(effect);

    m_fadeAnimation->setTargetObject(effect);
    m_fadeAnimation->setPropertyName("opacity");
    m_fadeAnimation->setStartValue(1.0);
    m_fadeAnimation->setEndValue(0.0);
    m_fadeAnimation->start();

    m_fadeAnimation->setProperty("targetPage", toPage);
}

void NewProjectDialog::onAnimationFinished()
{
    int targetPage = m_fadeAnimation->property("targetPage").toInt();
    ui->stackedWidget->setCurrentIndex(targetPage);

    QWidget* currentWidget = ui->stackedWidget->currentWidget();
    currentWidget->setGraphicsEffect(nullptr);
}

void NewProjectDialog::onNextButtonClicked()
{
    if (!validateCurrentPage()) {
        return;
    }

    if (m_currentPageIndex == ui->stackedWidget->count() - 1) {
        // 完成创建
        collectAllData();

        if (createProjectFile()) {
            emit projectCreated(m_projectInfo);
            accept();
        }
    } else {
        int nextPage = m_currentPageIndex + 1;
        animatePageTransition(m_currentPageIndex, nextPage);
        ui->navigationList->setCurrentRow(nextPage);
        m_currentPageIndex = nextPage;

        ui->backButton->setEnabled(true);

        if (nextPage == ui->stackedWidget->count() - 1) {
            ui->nextButton->setText(tr("完成创建"));
            ui->nextButton->setIcon(QIcon(":/new/prefix1/Resource/complete_icon.png"));
        }

        updateStatusLabel();
        updateProgressBar();
    }
}

void NewProjectDialog::onBackButtonClicked()
{
    if (m_currentPageIndex > 0) {
        int prevPage = m_currentPageIndex - 1;
        animatePageTransition(m_currentPageIndex, prevPage);
        ui->navigationList->setCurrentRow(prevPage);
        m_currentPageIndex = prevPage;

        ui->backButton->setEnabled(prevPage > 0);
        ui->nextButton->setText(tr("下一步 >>"));
        ui->nextButton->setIcon(QIcon(":/new/prefix1/Resource/next_icon.png"));

        updateStatusLabel();
        updateProgressBar();
    }
}

void NewProjectDialog::onCancelButtonClicked()
{
    if (m_isModified) {
        QMessageBox::StandardButton reply = QMessageBox::question(this, tr("取消确认"),
                                                                  tr("您确定要取消创建项目吗？所有已输入的数据将丢失。"),
                                                                  QMessageBox::Yes | QMessageBox::No);
        if (reply == QMessageBox::Yes) {
            reject();
        }
    } else {
        reject();
    }
}

void NewProjectDialog::onHelpButtonClicked()
{
    QString helpTitle;
    QString helpContent;

    switch (m_currentPageIndex) {
    case 0:
        helpTitle = tr("项目基本信息帮助");
        helpContent = tr("在此页面中，您可以设置试井项目的基本信息和项目保存路径。\n\n"
                         "快速配置预设可以帮助您快速设置不同类型储层的典型参数。\n"
                         "项目文件将保存为 .wtproject 格式，包含所有项目配置信息。");
        break;
    case 1:
        helpTitle = tr("井筒参数帮助");
        helpContent = tr("设置井筒的几何参数和完井信息。\n\n"
                         "所有参数均使用国际标准单位制：\n"
                         "• 长度：米 (m)\n"
                         "• 管柱尺寸：毫米 (mm)\n"
                         "• 压力：兆帕 (MPa)");
        break;
    case 2:
        helpTitle = tr("压裂参数帮助");
        helpContent = tr("多段压裂水平井的关键压裂参数设置。\n\n"
                         "系统将自动计算压裂改造体积和相关参数。");
        break;
    case 3:
        helpTitle = tr("储层参数帮助");
        helpContent = tr("储层的孔渗特性和边界条件设置。\n\n"
                         "根据储层类型选择合适的介质模型。");
        break;
    case 4:
        helpTitle = tr("流体特性帮助");
        helpContent = tr("流体的物理化学性质参数。\n\n"
                         "系统将自动计算总压缩系数。");
        break;
    case 5:
        helpTitle = tr("分析设置帮助");
        helpContent = tr("试井解释的分析方法和输出选项。\n\n"
                         "建议使用智能自动匹配方法。");
        break;
    default:
        helpTitle = tr("试井解释软件");
        helpContent = tr("专业的试井解释分析软件");
        break;
    }

    QMessageBox msgBox;
    msgBox.setWindowTitle(helpTitle);
    msgBox.setText(helpContent);
    msgBox.setIcon(QMessageBox::Information);
    msgBox.setStyleSheet("QLabel{color: black; font-size: 12px;}");
    msgBox.exec();
}

bool NewProjectDialog::validateCurrentPage()
{
    // 在初始化阶段不进行验证，避免错误提示
    if (m_isInitializing) {
        return true;
    }

    bool isValid = false;

    switch (m_currentPageIndex) {
    case 0: isValid = validateBasicInfo(); break;
    case 1: isValid = validateWellParams(); break;
    case 2: isValid = validateFractureParams(); break;
    case 3: isValid = validateReservoirParams(); break;
    case 4: isValid = validateFluidProperties(); break;
    case 5: isValid = validateAnalysisSettings(); break;
    default: isValid = true; break;
    }

    m_pageValidationStatus[m_currentPageIndex] = isValid;
    return isValid;
}

bool NewProjectDialog::validateCurrentPageSilently()
{
    // 静默验证，不显示错误对话框
    bool isValid = false;

    switch (m_currentPageIndex) {
    case 0: isValid = validateBasicInfoSilently(); break;
    case 1: isValid = validateWellParamsSilently(); break;
    case 2: isValid = validateFractureParamsSilently(); break;
    case 3: isValid = validateReservoirParamsSilently(); break;
    case 4: isValid = validateFluidPropertiesSilently(); break;
    case 5: isValid = validateAnalysisSettingsSilently(); break;
    default: isValid = true; break;
    }

    m_pageValidationStatus[m_currentPageIndex] = isValid;
    return isValid;
}

bool NewProjectDialog::validateBasicInfo()
{
    if (ui->projectNameEdit->text().trimmed().isEmpty()) {
        QMessageBox::warning(this, tr("验证错误"), tr("请输入项目名称。"));
        ui->projectNameEdit->setFocus();
        return false;
    }

    if (ui->wellNameEdit->text().trimmed().isEmpty()) {
        QMessageBox::warning(this, tr("验证错误"), tr("请输入井名。"));
        ui->wellNameEdit->setFocus();
        return false;
    }

    if (ui->projectPathEdit->text().trimmed().isEmpty()) {
        QMessageBox::warning(this, tr("验证错误"), tr("请选择项目保存路径。"));
        return false;
    }

    if (ui->testDurationSpin->value() <= 0) {
        QMessageBox::warning(this, tr("验证错误"), tr("测试时长必须大于零。"));
        ui->testDurationSpin->setFocus();
        return false;
    }

    return true;
}

bool NewProjectDialog::validateBasicInfoSilently()
{
    return !ui->projectNameEdit->text().trimmed().isEmpty() &&
           !ui->wellNameEdit->text().trimmed().isEmpty() &&
           !ui->projectPathEdit->text().trimmed().isEmpty() &&
           ui->testDurationSpin->value() > 0;
}

bool NewProjectDialog::validateWellParams()
{
    if (ui->wellRadiusSpin->value() <= 0) {
        QMessageBox::warning(this, tr("验证错误"), tr("井半径必须大于零。"));
        ui->wellRadiusSpin->setFocus();
        return false;
    }

    if (ui->wellTypeCombo->currentText() == tr("水平井") && ui->horizontalLengthSpin->value() <= 0) {
        QMessageBox::warning(this, tr("验证错误"), tr("水平井的水平段长度必须大于零。"));
        ui->horizontalLengthSpin->setFocus();
        return false;
    }

    return true;
}

bool NewProjectDialog::validateWellParamsSilently()
{
    bool isValid = ui->wellRadiusSpin->value() > 0;

    if (ui->wellTypeCombo->currentText() == tr("水平井")) {
        isValid = isValid && ui->horizontalLengthSpin->value() > 0;
    }

    return isValid;
}

bool NewProjectDialog::validateFractureParams()
{
    if (ui->fractureTypeCombo->currentText().contains(tr("压裂"))) {
        if (ui->fractureStagesSpin->value() <= 0) {
            QMessageBox::warning(this, tr("验证错误"), tr("压裂段数必须大于零。"));
            ui->fractureStagesSpin->setFocus();
            return false;
        }

        if (ui->fractureHalfLengthSpin->value() <= 0) {
            QMessageBox::warning(this, tr("验证错误"), tr("裂缝半长必须大于零。"));
            ui->fractureHalfLengthSpin->setFocus();
            return false;
        }
    }
    return true;
}

bool NewProjectDialog::validateFractureParamsSilently()
{
    if (ui->fractureTypeCombo->currentText().contains(tr("压裂"))) {
        return ui->fractureStagesSpin->value() > 0 && ui->fractureHalfLengthSpin->value() > 0;
    }
    return true;
}

bool NewProjectDialog::validateReservoirParams()
{
    if (ui->porositySpin->value() <= 0 || ui->porositySpin->value() >= 1) {
        QMessageBox::warning(this, tr("验证错误"), tr("孔隙度必须在0到1之间。"));
        ui->porositySpin->setFocus();
        return false;
    }

    if (ui->permeabilitySpin->value() <= 0) {
        QMessageBox::warning(this, tr("验证错误"), tr("渗透率必须大于零。"));
        ui->permeabilitySpin->setFocus();
        return false;
    }

    if (ui->compressibilitySpin->value() <= 0) {
        QMessageBox::warning(this, tr("验证错误"), tr("岩石压缩率必须大于零。"));
        ui->compressibilitySpin->setFocus();
        return false;
    }

    return true;
}

bool NewProjectDialog::validateReservoirParamsSilently()
{
    return ui->porositySpin->value() > 0 && ui->porositySpin->value() < 1 &&
           ui->permeabilitySpin->value() > 0 &&
           ui->compressibilitySpin->value() > 0;
}

bool NewProjectDialog::validateFluidProperties()
{
    if (!ui->oilCheckBox->isChecked() && !ui->gasCheckBox->isChecked() && !ui->waterCheckBox->isChecked()) {
        QMessageBox::warning(this, tr("验证错误"), tr("必须至少选择一种流体相。"));
        return false;
    }

    if (ui->viscositySpin->value() <= 0) {
        QMessageBox::warning(this, tr("验证错误"), tr("粘度必须大于零。"));
        ui->viscositySpin->setFocus();
        return false;
    }

    if (ui->fvfSpin->value() <= 0) {
        QMessageBox::warning(this, tr("验证错误"), tr("体积系数必须大于零。"));
        ui->fvfSpin->setFocus();
        return false;
    }

    return true;
}

bool NewProjectDialog::validateFluidPropertiesSilently()
{
    return (ui->oilCheckBox->isChecked() || ui->gasCheckBox->isChecked() || ui->waterCheckBox->isChecked()) &&
           ui->viscositySpin->value() > 0 &&
           ui->fvfSpin->value() > 0;
}

bool NewProjectDialog::validateAnalysisSettings()
{
    return true;
}

bool NewProjectDialog::validateAnalysisSettingsSilently()
{
    return true;
}

void NewProjectDialog::onDataChanged()
{
    if (!m_isApplyingPreset && !m_isInitializing) {
        m_isModified = true;
        m_validationTimer->start();
    }
}

void NewProjectDialog::updateStatusLabel()
{
    QString statusText;

    switch (m_currentPageIndex) {
    case 0: statusText = tr("设置试井项目基本信息和保存路径"); break;
    case 1: statusText = tr("配置井筒几何参数和完井信息"); break;
    case 2: statusText = tr("设置多段压裂参数"); break;
    case 3: statusText = tr("配置储层物性参数"); break;
    case 4: statusText = tr("设置流体特性参数"); break;
    case 5: statusText = tr("选择试井分析方法和输出选项"); break;
    default: statusText = tr("完成项目创建"); break;
    }

    ui->statusLabel->setText(statusText);
}

void NewProjectDialog::updateProgressBar()
{
    int validPages = 0;
    for (auto it = m_pageValidationStatus.begin(); it != m_pageValidationStatus.end(); ++it) {
        if (it.value()) validPages++;
    }

    int progress = (validPages * 100) / m_pageValidationStatus.size();
    ui->progressBar->setValue(progress);

    QString progressText = tr("完成进度: %1% (%2/%3)").arg(progress).arg(validPages).arg(m_pageValidationStatus.size());
    ui->progressLabel->setText(progressText);
}

void NewProjectDialog::updateControlsVisibility()
{
    bool isHorizontal = (ui->wellTypeCombo->currentText() == tr("水平井"));
    ui->horizontalLengthLabel->setVisible(isHorizontal);
    ui->horizontalLengthSpin->setVisible(isHorizontal);

    bool isFractured = !ui->fractureTypeCombo->currentText().contains(tr("未压裂"));
    ui->fractureParamsGroup->setVisible(isFractured);

    QString model = ui->reservoirModelCombo->currentText();
    bool isDualPorosity = model.contains(tr("双重")) || model.contains(tr("分形"));

    ui->matrixPorosityLabel->setVisible(isDualPorosity);
    ui->matrixPorositySpin->setVisible(isDualPorosity);
    ui->fracturePorosityLabel->setVisible(isDualPorosity);
    ui->fracturePorositySpin->setVisible(isDualPorosity);

    ui->matrixPermeabilityLabel->setVisible(isDualPorosity);
    ui->matrixPermeabilitySpin->setVisible(isDualPorosity);
    ui->fracturePermeabilityLabel->setVisible(isDualPorosity);
    ui->fracturePermeabilitySpin->setVisible(isDualPorosity);
}

void NewProjectDialog::collectAllData()
{
    // 收集所有表单数据到ProjectInfo结构体
    m_projectInfo.projectName = ui->projectNameEdit->text();
    m_projectInfo.wellName = ui->wellNameEdit->text();
    m_projectInfo.fieldName = ui->fieldNameEdit->text();
    m_projectInfo.engineerName = ui->engineerEdit->text();
    m_projectInfo.creationDate = ui->dateTimeEdit->dateTime();
    m_projectInfo.projectPath = ui->projectPathEdit->text();

    m_projectInfo.wellType = ui->wellTypeCombo->currentText();
    m_projectInfo.reservoirType = ui->reservoirTypeCombo->currentText();
    m_projectInfo.fractureType = ui->fractureTypeCombo->currentText();

    m_projectInfo.testType = ui->testTypeCombo->currentText();
    m_projectInfo.analysisType = ui->analysisTypeCombo->currentText();
    m_projectInfo.testDate = ui->testDateEdit->date();
    m_projectInfo.testDuration = ui->testDurationSpin->value();
    m_projectInfo.testDurationUnit = "小时";

    m_projectInfo.wellRadius = ui->wellRadiusSpin->value();
    m_projectInfo.horizontalLength = ui->horizontalLengthSpin->value();
    m_projectInfo.perforationTopDepth = ui->perforationTopSpin->value();
    m_projectInfo.perforationBottomDepth = ui->perforationBottomSpin->value();

    m_projectInfo.fractureStages = ui->fractureStagesSpin->value();
    m_projectInfo.fractureSpacing = ui->fractureSpacingSpin->value();
    m_projectInfo.fractureHalfLength = ui->fractureHalfLengthSpin->value();
    m_projectInfo.fractureConductivity = ui->fractureConductivitySpin->value();

    m_projectInfo.tuningSize = ui->tuningSpin->value();
    m_projectInfo.casingSize = ui->casingSpin->value();
    m_projectInfo.completionType = ui->completionCombo->currentText();
    m_projectInfo.bhpGaugeDepth = ui->bhpGaugeSpin->value();

    m_projectInfo.reservoirModel = ui->reservoirModelCombo->currentText();
    m_projectInfo.initialPressure = ui->reservoirPressureSpin->value();
    m_projectInfo.porosity = ui->porositySpin->value();
    m_projectInfo.matrixPorosity = ui->matrixPorositySpin->value();
    m_projectInfo.fracturePorosity = ui->fracturePorositySpin->value();
    m_projectInfo.reservoirTemp = ui->reservoirTempSpin->value();
    m_projectInfo.permeability = ui->permeabilitySpin->value();
    m_projectInfo.matrixPermeability = ui->matrixPermeabilitySpin->value();
    m_projectInfo.fracturePermeability = ui->fracturePermeabilitySpin->value();
    m_projectInfo.rockCompressibility = ui->compressibilitySpin->value();

    m_projectInfo.boundaryType = ui->boundaryTypeCombo->currentText();
    m_projectInfo.distanceToBoundary = ui->distanceToBoundarySpin->value();
    m_projectInfo.anisotropyRatio = ui->anisotropySpin->value();
    m_projectInfo.formationDip = ui->formationDipSpin->value();

    m_projectInfo.fluidSystem = ui->fluidSystemCombo->currentText();
    m_projectInfo.referencePressure = ui->referencePressureSpin->value();
    m_projectInfo.hasOil = ui->oilCheckBox->isChecked();
    m_projectInfo.hasGas = ui->gasCheckBox->isChecked();
    m_projectInfo.hasWater = ui->waterCheckBox->isChecked();
    m_projectInfo.viscosity = ui->viscositySpin->value();
    m_projectInfo.fvf = ui->fvfSpin->value();
    m_projectInfo.fluidCompressibility = ui->compressibilityFluidSpin->value();
    m_projectInfo.totalCompressibility = ui->totalCompressibilitySpin->value();

    m_projectInfo.analysisModel = ui->analysisModelCombo->currentText();
    m_projectInfo.flowRegime = ui->flowRegimeCombo->currentText();
    m_projectInfo.solutionMethod = ui->methodCombo->currentText();
    m_projectInfo.optimizationMethod = ui->optimizationMethodCombo->currentText();
    m_projectInfo.pressureDifferenceDefinition = ui->pressureDifferenceDefinitionCombo->currentText();
    m_projectInfo.timeFormat = ui->timeFormatCombo->currentText();
    m_projectInfo.primaryPlotType = ui->primaryPlotCombo->currentText();
    m_projectInfo.showGrid = ui->showGridCheckBox->isChecked();
    m_projectInfo.autoAdjustScale = ui->adjustScaleCheckBox->isChecked();
    m_projectInfo.generateReport = ui->generateReportCheckBox->isChecked();
    m_projectInfo.saveParameters = ui->saveParametersCheckBox->isChecked();

    m_projectInfo.enableTransientFlow = ui->enableTransientFlowCheckBox->isChecked();
    m_projectInfo.enableNonDarcyFlow = ui->enableNonDarcyFlowCheckBox->isChecked();

    qDebug() << "项目数据收集完成：" << m_projectInfo.projectName;
}

bool NewProjectDialog::createProjectFile()
{
    QString fileName = generateProjectFileName();
    QString fullPath = QDir(m_projectInfo.projectPath).absoluteFilePath(fileName);

    // 确保目录存在
    QDir dir(m_projectInfo.projectPath);
    if (!dir.exists()) {
        if (!dir.mkpath(".")) {
            QMessageBox::critical(this, tr("错误"), tr("无法创建项目目录：%1").arg(m_projectInfo.projectPath));
            return false;
        }
    }

    return saveProjectToFile(fullPath);
}

QString NewProjectDialog::generateProjectFileName()
{
    QString baseName = m_projectInfo.projectName;
    // 移除文件名中不允许的字符
    baseName = baseName.replace(QRegularExpression("[<>:\"/\\|?*]"), "_");

    QString fileName = baseName + ".wtproject";

    // 如果文件已存在，添加序号
    QString fullPath = QDir(m_projectInfo.projectPath).absoluteFilePath(fileName);
    int counter = 1;
    while (QFile::exists(fullPath)) {
        fileName = QString("%1_%2.wtproject").arg(baseName).arg(counter);
        fullPath = QDir(m_projectInfo.projectPath).absoluteFilePath(fileName);
        counter++;
    }

    return fileName;
}

bool NewProjectDialog::saveProjectToFile(const QString &filePath)
{
    QJsonObject projectJson;

    // 基本信息
    projectJson["projectName"] = m_projectInfo.projectName;
    projectJson["wellName"] = m_projectInfo.wellName;
    projectJson["fieldName"] = m_projectInfo.fieldName;
    projectJson["engineerName"] = m_projectInfo.engineerName;
    projectJson["creationDate"] = m_projectInfo.creationDate.toString(Qt::ISODate);
    projectJson["projectPath"] = m_projectInfo.projectPath;

    // 井型信息
    projectJson["wellType"] = m_projectInfo.wellType;
    projectJson["reservoirType"] = m_projectInfo.reservoirType;
    projectJson["fractureType"] = m_projectInfo.fractureType;

    // 测试信息
    projectJson["testType"] = m_projectInfo.testType;
    projectJson["analysisType"] = m_projectInfo.analysisType;
    projectJson["testDate"] = m_projectInfo.testDate.toString(Qt::ISODate);
    projectJson["testDuration"] = m_projectInfo.testDuration;
    projectJson["testDurationUnit"] = m_projectInfo.testDurationUnit;

    // 井筒参数
    projectJson["wellRadius"] = m_projectInfo.wellRadius;
    projectJson["horizontalLength"] = m_projectInfo.horizontalLength;
    projectJson["perforationTopDepth"] = m_projectInfo.perforationTopDepth;
    projectJson["perforationBottomDepth"] = m_projectInfo.perforationBottomDepth;

    // 压裂参数
    projectJson["fractureStages"] = m_projectInfo.fractureStages;
    projectJson["fractureSpacing"] = m_projectInfo.fractureSpacing;
    projectJson["fractureHalfLength"] = m_projectInfo.fractureHalfLength;
    projectJson["fractureConductivity"] = m_projectInfo.fractureConductivity;

    // 储层参数
    projectJson["reservoirModel"] = m_projectInfo.reservoirModel;
    projectJson["initialPressure"] = m_projectInfo.initialPressure;
    projectJson["porosity"] = m_projectInfo.porosity;
    projectJson["matrixPorosity"] = m_projectInfo.matrixPorosity;
    projectJson["fracturePorosity"] = m_projectInfo.fracturePorosity;
    projectJson["reservoirTemp"] = m_projectInfo.reservoirTemp;
    projectJson["permeability"] = m_projectInfo.permeability;
    projectJson["matrixPermeability"] = m_projectInfo.matrixPermeability;
    projectJson["fracturePermeability"] = m_projectInfo.fracturePermeability;
    projectJson["rockCompressibility"] = m_projectInfo.rockCompressibility;

    // 流体参数
    projectJson["fluidSystem"] = m_projectInfo.fluidSystem;
    projectJson["referencePressure"] = m_projectInfo.referencePressure;
    projectJson["hasOil"] = m_projectInfo.hasOil;
    projectJson["hasGas"] = m_projectInfo.hasGas;
    projectJson["hasWater"] = m_projectInfo.hasWater;
    projectJson["viscosity"] = m_projectInfo.viscosity;
    projectJson["fvf"] = m_projectInfo.fvf;
    projectJson["fluidCompressibility"] = m_projectInfo.fluidCompressibility;
    projectJson["totalCompressibility"] = m_projectInfo.totalCompressibility;

    // 分析设置
    projectJson["analysisModel"] = m_projectInfo.analysisModel;
    projectJson["flowRegime"] = m_projectInfo.flowRegime;
    projectJson["solutionMethod"] = m_projectInfo.solutionMethod;
    projectJson["optimizationMethod"] = m_projectInfo.optimizationMethod;
    projectJson["generateReport"] = m_projectInfo.generateReport;
    projectJson["saveParameters"] = m_projectInfo.saveParameters;

    // 高级设置
    projectJson["enableTransientFlow"] = m_projectInfo.enableTransientFlow;
    projectJson["enableNonDarcyFlow"] = m_projectInfo.enableNonDarcyFlow;

    // 写入文件
    QJsonDocument document(projectJson);
    QFile file(filePath);
    if (!file.open(QIODevice::WriteOnly)) {
        QMessageBox::critical(this, tr("错误"), tr("无法创建项目文件：%1\n%2").arg(filePath).arg(file.errorString()));
        return false;
    }

    file.write(document.toJson());
    file.close();

    // 更新项目信息中的完整路径
    m_projectInfo.projectPath = filePath;

    QMessageBox::information(this, tr("成功"), tr("项目文件已创建：\n%1").arg(filePath));
    return true;
}

void NewProjectDialog::onSelectProjectPath()
{
    QString defaultPath = QStandardPaths::writableLocation(QStandardPaths::DocumentsLocation);
    QString selectedPath = QFileDialog::getExistingDirectory(this, tr("选择项目保存路径"), defaultPath);

    if (!selectedPath.isEmpty()) {
        ui->projectPathEdit->setText(selectedPath);
    }
}

void NewProjectDialog::onFluidPhaseChanged()
{
    if (m_isApplyingPreset || m_isInitializing) return;

    if (!ui->oilCheckBox->isChecked() && !ui->gasCheckBox->isChecked() && !ui->waterCheckBox->isChecked()) {
        ui->oilCheckBox->setChecked(true);
        QMessageBox::warning(this, tr("流体相选择"), tr("必须至少选择一种流体相。"));
    }
}

void NewProjectDialog::onAnalysisTypeChanged()
{
    if (m_isApplyingPreset || m_isInitializing) return;
    // 根据分析类型调整相关设置
}

NewProjectDialog::ProjectInfo NewProjectDialog::getProjectInfo() const
{
    return m_projectInfo;
}

void NewProjectDialog::closeEvent(QCloseEvent *event)
{
    if (m_isModified) {
        QMessageBox::StandardButton reply = QMessageBox::question(this, tr("关闭确认"),
                                                                  tr("您确定要关闭此对话框吗？所有已输入的数据将丢失。"),
                                                                  QMessageBox::Yes | QMessageBox::No);
        if (reply == QMessageBox::No) {
            event->ignore();
            return;
        }
    }
    event->accept();
}

void NewProjectDialog::showEvent(QShowEvent *event)
{
    QDialog::showEvent(event);
    updateProgressBar();
}
