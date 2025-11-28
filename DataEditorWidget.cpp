#include "dataeditorwidget.h"
#include "ui_dataeditorwidget.h"
#include "PressureDerivativeCalculator.h"
#include <QDebug>
#include <QFileDialog>
#include <QMessageBox>
#include <QFile>
#include <QTextStream>
#include <QStandardItemModel>
#include <QHeaderView>
#include <QStyledItemDelegate>
#include <QPainter>
#include <QInputDialog>
#include <QTimer>
#include <QFileInfo>
#include <QRegularExpression>
#include <QPropertyAnimation>
#include <QGraphicsDropShadowEffect>
#include <QSortFilterProxyModel>
#include <QUndoStack>
#include <QProgressDialog>
#include <QApplication>
#include <QClipboard>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QPrintDialog>
#include <QPrinter>
#include <QTextDocument>
#include <QTextDocumentWriter>
#include <QMimeData>
#include <QScrollBar>
#include <QSet>
#include <QList>
#include <QFormLayout>
#include <QGroupBox>
#include <QCheckBox>
#include <QSpinBox>
#include <QDateTimeEdit>
#include <QMovie>
#include <QGraphicsOpacityEffect>
#include <QScrollArea>
#include <QCloseEvent>
#include <QLineEdit>
#include <QSplitter>
#include <QTabWidget>
#include <QButtonGroup>
#include <QRadioButton>
#include <QSlider>
#include <QDial>
#include <QTextEdit>
#include <QPlainTextEdit>
#include <cmath>
#include <algorithm>

// Qt6兼容性处理
#if QT_VERSION >= QT_VERSION_CHECK(6, 0, 0)
#include <QStringConverter>
#else
#include <QTextCodec>
#endif

// Windows Excel COM组件支持
#ifdef Q_OS_WIN
#include <QAxObject>
#endif

// ============================================================================
// 撤销重做命令实现
// ============================================================================

DataEditCommand::DataEditCommand(QStandardItemModel* model, QUndoCommand* parent)
    : QUndoCommand(parent), m_model(model)
{
}

CellEditCommand::CellEditCommand(QStandardItemModel* model, int row, int column,
                                 const QString& oldValue, const QString& newValue,
                                 QUndoCommand* parent)
    : DataEditCommand(model, parent), m_row(row), m_column(column),
    m_oldValue(oldValue), m_newValue(newValue)
{
    setText(QString("编辑单元格 (%1, %2)").arg(row + 1).arg(column + 1));
}

void CellEditCommand::undo()
{
    if (m_model && m_row < m_model->rowCount() && m_column < m_model->columnCount()) {
        QStandardItem* item = m_model->item(m_row, m_column);
        if (!item) {
            item = new QStandardItem();
            m_model->setItem(m_row, m_column, item);
        }
        item->setText(m_oldValue);
    }
}

void DataEditorWidget::loadDataWithConfig(const QString& filePath, const QString& fileType, const DataLoadConfigDialog::LoadConfig& config)
{
    qDebug() << "开始加载文件:" << filePath << "类型:" << fileType << "起始行:" << config.startRow;

    QFileInfo fileInfo(filePath);
    if (!fileInfo.exists() || !fileInfo.isReadable()) {
        showStyledMessageBox("文件加载失败",
                             QString("文件不存在或无法读取: %1").arg(filePath),
                             QMessageBox::Warning);
        return;
    }

    // 显示进度对话框
    showAnimatedProgress("加载数据文件", "正在读取文件数据，请稍候...");

    clearData();

    m_currentFilePath = filePath;
    m_currentFileType = fileType;
    ui->filePathLineEdit->setText(filePath);

    bool loadSuccess = false;
    QString errorMessage;

    updateProgress(20, "正在分析文件格式...");

    QString lowerType = fileType.toLower();
    if (lowerType == "txt" || lowerType == "csv") {
        loadSuccess = loadCsvFileWithConfig(filePath, config, errorMessage);
    } else {
        // 其他类型使用默认方法
        if (lowerType == "excel") {
            loadSuccess = loadExcelFileOptimized(filePath, errorMessage);
        } else if (lowerType == "json") {
            loadSuccess = loadJsonFile(filePath, errorMessage);
        } else {
            errorMessage = QString("不支持的文件类型: %1").arg(fileType);
        }
    }

    hideAnimatedProgress();

    if (loadSuccess) {
        updateStatus(QString("文件加载成功 - %1行 × %2列")
                         .arg(m_dataModel->rowCount())
                         .arg(m_dataModel->columnCount()), "success");

        setButtonsEnabled(true);
        m_dataModified = false;

        applyColumnStyles();
        optimizeColumnWidths();
        optimizeTableDisplay();

        // 弹出列定义对话框
        QTimer::singleShot(500, this, &DataEditorWidget::onDefineColumns);

        emitDataChanged();

        qDebug() << "文件加载成功，数据行数:" << m_dataModel->rowCount()
                 << "列数:" << m_dataModel->columnCount();
    } else {
        updateStatus("文件加载失败", "error");
        showStyledMessageBox("文件加载失败",
                             QString("无法加载文件: %1").arg(filePath),
                             QMessageBox::Critical,
                             errorMessage);
        qDebug() << "文件加载失败:" << errorMessage;
    }
}

void CellEditCommand::redo()
{
    if (m_model && m_row < m_model->rowCount() && m_column < m_model->columnCount()) {
        QStandardItem* item = m_model->item(m_row, m_column);
        if (!item) {
            item = new QStandardItem();
            m_model->setItem(m_row, m_column, item);
        }
        item->setText(m_newValue);
    }
}

RowEditCommand::RowEditCommand(QStandardItemModel* model, Operation op, int row,
                               const QStringList& rowData, QUndoCommand* parent)
    : DataEditCommand(model, parent), m_operation(op), m_row(row), m_rowData(rowData)
{
    QString operationText = (op == Insert) ? "插入" : "删除";
    setText(QString("%1行 %2").arg(operationText).arg(row + 1));
}

void RowEditCommand::undo()
{
    if (!m_model) return;

    if (m_operation == Insert) {
        if (m_row < m_model->rowCount()) {
            m_model->removeRow(m_row);
        }
    } else {
        m_model->insertRow(m_row);
        for (int col = 0; col < m_rowData.size(); ++col) {
            QStandardItem* item = new QStandardItem(m_rowData[col]);
            m_model->setItem(m_row, col, item);
        }
    }
}

void RowEditCommand::redo()
{
    if (!m_model) return;

    if (m_operation == Insert) {
        m_model->insertRow(m_row);
        for (int col = 0; col < m_model->columnCount(); ++col) {
            QStandardItem* item = new QStandardItem("");
            m_model->setItem(m_row, col, item);
        }
    } else {
        if (m_row < m_model->rowCount()) {
            m_rowData.clear();
            for (int col = 0; col < m_model->columnCount(); ++col) {
                QStandardItem* item = m_model->item(m_row, col);
                m_rowData.append(item ? item->text() : "");
            }
            m_model->removeRow(m_row);
        }
    }
}

ColumnEditCommand::ColumnEditCommand(QStandardItemModel* model, Operation op, int column,
                                     const QString& headerName, const QStringList& columnData,
                                     QUndoCommand* parent)
    : DataEditCommand(model, parent), m_operation(op), m_column(column),
    m_headerName(headerName), m_columnData(columnData)
{
    QString operationText = (op == Insert) ? "插入" : "删除";
    setText(QString("%1列 %2").arg(operationText).arg(column + 1));
}

void ColumnEditCommand::undo()
{
    if (!m_model) return;

    if (m_operation == Insert) {
        if (m_column < m_model->columnCount()) {
            m_model->removeColumn(m_column);
        }
    } else {
        m_model->insertColumn(m_column);
        QStandardItem* headerItem = new QStandardItem(m_headerName);
        m_model->setHorizontalHeaderItem(m_column, headerItem);

        for (int row = 0; row < m_columnData.size(); ++row) {
            QStandardItem* item = new QStandardItem(m_columnData[row]);
            m_model->setItem(row, m_column, item);
        }
    }
}

void ColumnEditCommand::redo()
{
    if (!m_model) return;

    if (m_operation == Insert) {
        m_model->insertColumn(m_column);
        QStandardItem* headerItem = new QStandardItem(m_headerName);
        m_model->setHorizontalHeaderItem(m_column, headerItem);

        for (int row = 0; row < m_model->rowCount(); ++row) {
            QStandardItem* item = new QStandardItem("");
            m_model->setItem(row, m_column, item);
        }
    } else {
        if (m_column < m_model->columnCount()) {
            QStandardItem* headerItem = m_model->horizontalHeaderItem(m_column);
            m_headerName = headerItem ? headerItem->text() : QString("列%1").arg(m_column + 1);

            m_columnData.clear();
            for (int row = 0; row < m_model->rowCount(); ++row) {
                QStandardItem* item = m_model->item(row, m_column);
                m_columnData.append(item ? item->text() : "");
            }

            m_model->removeColumn(m_column);
        }
    }
}

// ============================================================================
// 列定义对话框实现 - 优化版本
// ============================================================================

ColumnDefinitionDialog::ColumnDefinitionDialog(const QStringList& columnNames,
                                               const QList<ColumnDefinition>& definitions,
                                               QWidget* parent)
    : QDialog(parent), m_columnNames(columnNames), m_definitions(definitions)
{
    setupUI();
    if (m_definitions.isEmpty()) {
        // 创建默认定义
        for (const QString& name : columnNames) {
            ColumnDefinition def;
            def.name = name;
            def.type = WellTestColumnType::Custom;
            def.unit = "";
            def.description = "";
            m_definitions.append(def);
        }
    }

    // 设置初始值
    for (int i = 0; i < m_columnNames.size() && i < m_definitions.size(); ++i) {
        if (i < m_typeComboBoxes.size()) {
            m_typeComboBoxes[i]->setCurrentIndex(static_cast<int>(m_definitions[i].type));
            WellTestColumnType type = static_cast<WellTestColumnType>(m_typeComboBoxes[i]->currentIndex());
            updateUnitsForType(type, m_unitComboBoxes[i]);
        }
        if (i < m_unitComboBoxes.size()) {
            int unitIndex = m_unitComboBoxes[i]->findText(m_definitions[i].unit);
            if (unitIndex >= 0) {
                m_unitComboBoxes[i]->setCurrentIndex(unitIndex);
            }
        }
        if (i < m_requiredChecks.size()) {
            m_requiredChecks[i]->setChecked(m_definitions[i].isRequired);
        }
        // 更新预览
        updatePreviewLabel(i);
    }
}

// ============================================================================
// 列定义对话框实现 - 优化版本（修改部分）
// ============================================================================

void ColumnDefinitionDialog::setupUI()
{
    setWindowTitle("数据列定义");
    setModal(true);
    resize(750, 500);

    m_mainLayout = new QVBoxLayout(this);

    // 说明标签
    QLabel* infoLabel = new QLabel("为每列数据定义物理意义和单位，将替换原列名：");
    infoLabel->setStyleSheet("font-size: 14px; font-weight: bold; color: #2c3e50; margin: 10px;");
    m_mainLayout->addWidget(infoLabel);

    // 滚动区域
    QScrollArea* scrollArea = new QScrollArea;
    QWidget* scrollWidget = new QWidget;
    m_formLayout = new QFormLayout(scrollWidget);

    // 中文类型名称 - 添加日期、时刻选项，调整顺序
    QStringList typeNames = {
        "序号", "日期", "时刻", "时间", "压力", "温度", "流量", "深度", "粘度", "密度",
        "渗透率", "孔隙度", "井半径", "表皮系数", "距离", "体积", "压降", "自定义"
    };

    for (int i = 0; i < m_columnNames.size(); ++i) {
        // 创建行容器
        QWidget* rowWidget = new QWidget;
        QHBoxLayout* rowLayout = new QHBoxLayout(rowWidget);
        rowLayout->setContentsMargins(0, 0, 0, 0);

        // 原列名标签
        QLabel* originalNameLabel = new QLabel(QString("原列名: %1").arg(m_columnNames[i]));
        originalNameLabel->setFixedWidth(150);
        originalNameLabel->setStyleSheet("font-weight: bold; color: #6c757d; font-size: 11px;");
        rowLayout->addWidget(originalNameLabel);

        // 类型下拉框
        QComboBox* typeCombo = new QComboBox;
        typeCombo->addItems(typeNames);
        typeCombo->setFixedWidth(120);
        typeCombo->setCurrentIndex(17); // 默认选择"自定义"
        typeCombo->setEditable(false); // 初始状态不可编辑
        rowLayout->addWidget(typeCombo);

        connect(typeCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
                this, &ColumnDefinitionDialog::onTypeChanged);

        // 连接文本编辑信号（当下拉框变为可编辑时使用）
        connect(typeCombo, &QComboBox::editTextChanged,
                this, &ColumnDefinitionDialog::onCustomTypeTextChanged);

        m_typeComboBoxes.append(typeCombo);

        // 创建自定义类型输入框（保留但不再使用）
        QLineEdit* customTypeEdit = new QLineEdit;
        customTypeEdit->setFixedWidth(120);
        customTypeEdit->setPlaceholderText("输入自定义类型");
        customTypeEdit->setVisible(false);
        m_customTypeEdits.append(customTypeEdit);

        // 单位下拉框和自定义输入框容器
        QWidget* unitWidget = new QWidget;
        QVBoxLayout* unitLayout = new QVBoxLayout(unitWidget);
        unitLayout->setContentsMargins(0, 0, 0, 0);
        unitLayout->setSpacing(2);

        QComboBox* unitCombo = new QComboBox;
        unitCombo->setFixedWidth(100);
        unitCombo->setEditable(false);
        unitLayout->addWidget(unitCombo);

        QLineEdit* customUnitEdit = new QLineEdit;
        customUnitEdit->setFixedWidth(100);
        customUnitEdit->setPlaceholderText("输入单位");
        customUnitEdit->setVisible(false);
        unitLayout->addWidget(customUnitEdit);

        unitWidget->setFixedWidth(100);
        rowLayout->addWidget(unitWidget);

        connect(unitCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
                this, &ColumnDefinitionDialog::onUnitChanged);
        connect(customUnitEdit, &QLineEdit::textChanged,
                this, &ColumnDefinitionDialog::onCustomUnitChanged);

        m_unitComboBoxes.append(unitCombo);
        m_customUnitEdits.append(customUnitEdit);

        // 必需复选框
        QCheckBox* requiredCheck = new QCheckBox("必需");
        m_requiredChecks.append(requiredCheck);
        rowLayout->addWidget(requiredCheck);

        // 预览标签
        QLabel* previewLabel = new QLabel("自定义\\-");
        previewLabel->setFixedWidth(120);
        previewLabel->setStyleSheet("color: #28a745; font-weight: bold; font-size: 11px;");
        m_previewLabels.append(previewLabel);
        rowLayout->addWidget(previewLabel);

        m_formLayout->addRow(rowWidget);

        // 初始化单位选项
        updateUnitsForType(WellTestColumnType::Custom, unitCombo);
        if (unitCombo->count() > 0) {
            unitCombo->setCurrentIndex(0);
        }

        // 设置预览标签的初始值
        updatePreviewLabel(i);
    }

    scrollArea->setWidget(scrollWidget);
    scrollArea->setWidgetResizable(true);
    m_mainLayout->addWidget(scrollArea);

    // 按钮区域
    QHBoxLayout* buttonLayout = new QHBoxLayout;

    QPushButton* presetBtn = new QPushButton("自动识别");
    presetBtn->setStyleSheet("QPushButton { background-color: #4a90e2; color: white; border: none; border-radius: 4px; padding: 8px 16px; }");
    connect(presetBtn, &QPushButton::clicked, this, &ColumnDefinitionDialog::loadPresetDefinitions);
    buttonLayout->addWidget(presetBtn);

    QPushButton* resetBtn = new QPushButton("重置");
    resetBtn->setStyleSheet("QPushButton { background-color: #fd7e14; color: white; border: none; border-radius: 4px; padding: 8px 16px; }");
    connect(resetBtn, &QPushButton::clicked, this, &ColumnDefinitionDialog::onResetClicked);
    buttonLayout->addWidget(resetBtn);

    buttonLayout->addStretch();

    QPushButton* okBtn = new QPushButton("确定");
    okBtn->setStyleSheet("QPushButton { background-color: #28a745; color: white; border: none; border-radius: 4px; padding: 8px 16px; }");
    connect(okBtn, &QPushButton::clicked, this, &QDialog::accept);
    buttonLayout->addWidget(okBtn);

    QPushButton* cancelBtn = new QPushButton("取消");
    cancelBtn->setStyleSheet("QPushButton { background-color: #fd7e14; color: white; border: none; border-radius: 4px; padding: 8px 16px; }");
    connect(cancelBtn, &QPushButton::clicked, this, &QDialog::reject);
    buttonLayout->addWidget(cancelBtn);

    m_mainLayout->addLayout(buttonLayout);
}


void ColumnDefinitionDialog::onTypeChanged(int index)
{
    // 找到触发信号的组合框
    QComboBox* senderCombo = qobject_cast<QComboBox*>(sender());
    if (!senderCombo) return;

    int comboIndex = -1;
    for (int i = 0; i < m_typeComboBoxes.size(); ++i) {
        if (m_typeComboBoxes[i] == senderCombo) {
            comboIndex = i;
            break;
        }
    }

    if (comboIndex >= 0 && comboIndex < m_unitComboBoxes.size()) {
        WellTestColumnType type = static_cast<WellTestColumnType>(index);

        // 处理自定义类型的显示
        bool isCustomType = (index == 15); // "自定义"是最后一个选项

        if (isCustomType) {
            // 让下拉框立即变为可编辑状态
            senderCombo->setEditable(true);
            senderCombo->setCurrentText("自定义");

            // 立即聚焦并选择文本，让用户可以直接编辑
            QTimer::singleShot(10, [senderCombo]() {
                if (senderCombo->lineEdit()) {
                    senderCombo->lineEdit()->selectAll();
                    senderCombo->lineEdit()->setFocus();
                }
            });
        } else {
            // 恢复为不可编辑状态
            senderCombo->setEditable(false);
        }

        // 隐藏自定义输入框（不再需要）
        if (comboIndex < m_customTypeEdits.size()) {
            m_customTypeEdits[comboIndex]->setVisible(false);
        }

        updateUnitsForType(type, m_unitComboBoxes[comboIndex]);
        updatePreviewLabel(comboIndex);
    }
}

void ColumnDefinitionDialog::onCustomTypeTextChanged(const QString& text)
{
    // 找到触发信号的组合框
    QComboBox* senderCombo = qobject_cast<QComboBox*>(sender());
    if (!senderCombo) return;

    int comboIndex = -1;
    for (int i = 0; i < m_typeComboBoxes.size(); ++i) {
        if (m_typeComboBoxes[i] == senderCombo) {
            comboIndex = i;
            break;
        }
    }

    if (comboIndex >= 0) {
        updatePreviewLabel(comboIndex);
    }
}

void ColumnDefinitionDialog::onUnitChanged(int index)
{
    // 找到触发信号的组合框
    QComboBox* senderCombo = qobject_cast<QComboBox*>(sender());
    if (!senderCombo) return;

    int comboIndex = -1;
    for (int i = 0; i < m_unitComboBoxes.size(); ++i) {
        if (m_unitComboBoxes[i] == senderCombo) {
            comboIndex = i;
            break;
        }
    }

    if (comboIndex >= 0) {
        // 处理自定义单位的显示
        QString unitText = senderCombo->itemText(index);
        bool isCustomUnit = (unitText == "自定义");
        if (comboIndex < m_customUnitEdits.size()) {
            m_customUnitEdits[comboIndex]->setVisible(isCustomUnit);
            if (isCustomUnit) {
                senderCombo->setEditable(true);
                senderCombo->setCurrentText("");
            } else {
                senderCombo->setEditable(false);
            }
        }

        updatePreviewLabel(comboIndex);
    }
}

void ColumnDefinitionDialog::onCustomTypeChanged(const QString& text)
{
    // 这个方法现在不再需要，因为我们直接使用ComboBox的editTextChanged信号
}

void ColumnDefinitionDialog::onCustomUnitChanged(const QString& text)
{
    QLineEdit* senderEdit = qobject_cast<QLineEdit*>(sender());
    if (!senderEdit) return;

    int editIndex = -1;
    for (int i = 0; i < m_customUnitEdits.size(); ++i) {
        if (m_customUnitEdits[i] == senderEdit) {
            editIndex = i;
            break;
        }
    }

    if (editIndex >= 0) {
        updatePreviewLabel(editIndex);
    }
}

void ColumnDefinitionDialog::updateUnitsForType(WellTestColumnType type, QComboBox* unitCombo)
{
    unitCombo->clear();

    switch (type) {
    case WellTestColumnType::SerialNumber:
        unitCombo->addItems({"-", "个", "项", "自定义"});
        break;
    case WellTestColumnType::Date:  // 新增日期类型
        unitCombo->addItems({"-", "yyyy-MM-dd", "yyyy/MM/dd", "dd/MM/yyyy", "自定义"});
        break;
    case WellTestColumnType::TimeOfDay:  // 新增时刻类型
        unitCombo->addItems({"-", "hh:mm:ss", "hh:mm:ss.zzz", "hh:mm", "自定义"});
        break;
    case WellTestColumnType::Time:
        unitCombo->addItems({"h", "min", "s", "day", "-", "自定义"});
        break;
    case WellTestColumnType::Pressure:
        unitCombo->addItems({"MPa", "kPa", "Pa", "psi", "bar", "atm", "-", "自定义"});
        break;
    case WellTestColumnType::PressureDrop:
        unitCombo->addItems({"MPa", "kPa", "Pa", "psi", "bar", "atm", "-", "自定义"});
        break;
    case WellTestColumnType::Temperature:
        unitCombo->addItems({"°C", "°F", "K", "-", "自定义"});
        break;
    case WellTestColumnType::FlowRate:
        unitCombo->addItems({"m³/d", "m³/h", "L/s", "bbl/d", "ft³/d", "-", "自定义"});
        break;
    case WellTestColumnType::Depth:
        unitCombo->addItems({"m", "ft", "km", "mm", "-", "自定义"});
        break;
    case WellTestColumnType::Viscosity:
        unitCombo->addItems({"mPa·s", "cP", "Pa·s", "-", "自定义"});
        break;
    case WellTestColumnType::Density:
        unitCombo->addItems({"kg/m³", "g/cm³", "lb/ft³", "-", "自定义"});
        break;
    case WellTestColumnType::Permeability:
        unitCombo->addItems({"mD", "D", "μm²", "-", "自定义"});
        break;
    case WellTestColumnType::Porosity:
        unitCombo->addItems({"%", "fraction", "-", "自定义"});
        break;
    case WellTestColumnType::WellRadius:
        unitCombo->addItems({"m", "ft", "cm", "in", "-", "自定义"});
        break;
    case WellTestColumnType::SkinFactor:
        unitCombo->addItems({"dimensionless", "-", "自定义"});
        break;
    case WellTestColumnType::Distance:
        unitCombo->addItems({"m", "ft", "km", "mm", "-", "自定义"});
        break;
    case WellTestColumnType::Volume:
        unitCombo->addItems({"m³", "L", "bbl", "ft³", "-", "自定义"});
        break;
    default: // Custom
        unitCombo->addItems({"-", "个", "项", "次", "自定义"});
        break;
    }

    // 连接单位变化信号到预览更新
    connect(unitCombo, QOverload<int>::of(&QComboBox::currentIndexChanged),
            this, &ColumnDefinitionDialog::onUnitChanged, Qt::UniqueConnection);
}

void ColumnDefinitionDialog::updatePreviewLabel(int index)
{
    if (index < 0 || index >= m_previewLabels.size()) {
        return;
    }

    QString typeName;
    QString unitName;

    // 获取类型名称
    if (index < m_typeComboBoxes.size()) {
        QComboBox* typeCombo = m_typeComboBoxes[index];
        if (typeCombo->isEditable()) {
            // 如果是可编辑状态，使用当前编辑的文本
            typeName = typeCombo->currentText();
            if (typeName.isEmpty()) {
                typeName = "自定义";
            }
        } else {
            // 否则使用选中的项目文本
            typeName = typeCombo->currentText();
        }
    }

    // 获取单位名称
    if (index < m_unitComboBoxes.size() && index < m_customUnitEdits.size()) {
        if (m_customUnitEdits[index]->isVisible()) {
            unitName = m_customUnitEdits[index]->text();
        } else {
            QString unit = m_unitComboBoxes[index]->currentText();
            if (unit == "自定义") {
                unitName = m_unitComboBoxes[index]->currentText();
            } else {
                unitName = (unit == "-") ? "" : unit;
            }
        }
    }

    // 生成预览文本
    QString preview;
    if (unitName.isEmpty() || unitName == "-") {
        preview = typeName;
    } else {
        preview = QString("%1\\%2").arg(typeName).arg(unitName);
    }

    m_previewLabels[index]->setText(preview);
}

void ColumnDefinitionDialog::loadPresetDefinitions()
{
    // 智能识别列名并应用预设 - 增加日期、时刻识别
    for (int i = 0; i < m_columnNames.size(); ++i) {
        QString name = m_columnNames[i].toLower();
        int typeIndex = 17; // 默认"自定义"
        QString suggestedUnit = "-";

        if (name.contains("序号") || name.contains("编号") || name.contains("number") || name == "no" || name == "id") {
            typeIndex = 0; // 序号
            suggestedUnit = "-";
        } else if (name.contains("日期") || name.contains("date") || name.contains("年月日")) {
            typeIndex = 1; // 日期
            suggestedUnit = "yyyy-MM-dd";
        } else if (name.contains("时刻") || name.contains("时分秒") || name.contains("timeofday") || name.contains("clock")) {
            typeIndex = 2; // 时刻
            suggestedUnit = "hh:mm:ss";
        } else if (name.contains("time") || name.contains("时间") || name == "t") {
            typeIndex = 3; // 时间（调整索引）
            suggestedUnit = "h";
        } else if (name.contains("pressure") || name.contains("压力") || name == "p") {
            typeIndex = 4; // 压力
            suggestedUnit = "MPa";
        } else if (name.contains("temp") || name.contains("温度")) {
            typeIndex = 5; // 温度
            suggestedUnit = "°C";
        } else if (name.contains("flow") || name.contains("流量") || name == "q") {
            typeIndex = 6; // 流量
            suggestedUnit = "m³/d";
        } else if (name.contains("depth") || name.contains("深度")) {
            typeIndex = 7; // 深度
            suggestedUnit = "m";
        } else if (name.contains("viscosity") || name.contains("粘度")) {
            typeIndex = 8; // 粘度
            suggestedUnit = "mPa·s";
        } else if (name.contains("density") || name.contains("密度")) {
            typeIndex = 9; // 密度
            suggestedUnit = "kg/m³";
        } else if (name.contains("perm") || name.contains("渗透")) {
            typeIndex = 10; // 渗透率
            suggestedUnit = "mD";
        } else if (name.contains("porosity") || name.contains("孔隙")) {
            typeIndex = 11; // 孔隙度
            suggestedUnit = "%";
        } else if (name.contains("radius") || name.contains("半径")) {
            typeIndex = 12; // 井半径
            suggestedUnit = "m";
        } else if (name.contains("skin") || name.contains("表皮")) {
            typeIndex = 13; // 表皮系数
            suggestedUnit = "dimensionless";
        } else if (name.contains("distance") || name.contains("距离")) {
            typeIndex = 14; // 距离
            suggestedUnit = "m";
        } else if (name.contains("volume") || name.contains("体积")) {
            typeIndex = 15; // 体积
            suggestedUnit = "m³";
        } else if (name.contains("drop") || name.contains("降") || name.contains("差")) {
            typeIndex = 16; // 压降
            suggestedUnit = "MPa";
        }

        if (i < m_typeComboBoxes.size()) {
            // 恢复下拉框为不可编辑状态
            m_typeComboBoxes[i]->setEditable(false);

            m_typeComboBoxes[i]->setCurrentIndex(typeIndex);

            // 更新单位选项
            WellTestColumnType type = static_cast<WellTestColumnType>(typeIndex);
            updateUnitsForType(type, m_unitComboBoxes[i]);

            // 设置建议的单位
            if (i < m_unitComboBoxes.size()) {
                // 确保自定义单位输入框隐藏，显示下拉框
                if (i < m_customUnitEdits.size()) {
                    m_customUnitEdits[i]->setVisible(false);
                    m_unitComboBoxes[i]->setVisible(true);
                    m_unitComboBoxes[i]->setEditable(false);
                }

                int unitIndex = m_unitComboBoxes[i]->findText(suggestedUnit);
                if (unitIndex >= 0) {
                    m_unitComboBoxes[i]->setCurrentIndex(unitIndex);
                }
            }

            // 更新预览
            updatePreviewLabel(i);
        }
    }
}

void ColumnDefinitionDialog::onResetClicked()
{
    for (int i = 0; i < m_typeComboBoxes.size(); ++i) {
        // 重置类型为自定义
        m_typeComboBoxes[i]->setEditable(false);
        m_typeComboBoxes[i]->setCurrentIndex(17); // 设置为"自定义"

        // 重置单位
        updateUnitsForType(WellTestColumnType::Custom, m_unitComboBoxes[i]);
        if (i < m_customUnitEdits.size()) {
            m_customUnitEdits[i]->setVisible(false);
            m_customUnitEdits[i]->clear();
            m_unitComboBoxes[i]->setVisible(true);
            m_unitComboBoxes[i]->setEditable(false);
        }
        if (m_unitComboBoxes[i]->count() > 0) {
            m_unitComboBoxes[i]->setCurrentIndex(0); // 选择第一个单位
        }

        updatePreviewLabel(i);
    }

    for (auto check : m_requiredChecks) {
        check->setChecked(false);
    }
}

QList<ColumnDefinition> ColumnDefinitionDialog::getColumnDefinitions() const
{
    QList<ColumnDefinition> definitions;

    for (int i = 0; i < m_columnNames.size(); ++i) {
        ColumnDefinition def;

        // 获取类型名称
        QString typeName;
        if (i < m_typeComboBoxes.size()) {
            QComboBox* typeCombo = m_typeComboBoxes[i];
            if (typeCombo->isEditable()) {
                // 如果是可编辑状态，使用当前编辑的文本
                typeName = typeCombo->currentText();
                if (typeName.isEmpty()) {
                    typeName = "自定义";
                }
            } else {
                // 否则使用选中的项目文本
                typeName = typeCombo->currentText();
            }
        }

        // 获取单位名称
        QString unitName;
        if (i < m_unitComboBoxes.size() && i < m_customUnitEdits.size()) {
            if (m_customUnitEdits[i]->isVisible()) {
                unitName = m_customUnitEdits[i]->text();
            } else {
                QString unit = m_unitComboBoxes[i]->currentText();
                if (unit == "自定义") {
                    unitName = m_unitComboBoxes[i]->currentText();
                } else {
                    unitName = (unit == "-") ? "" : unit;
                }
            }
        }

        // 组合类型和单位作为新的列名
        if (unitName.isEmpty()) {
            def.name = typeName;
        } else {
            def.name = QString("%1\\%2").arg(typeName).arg(unitName);
        }

        if (i < m_typeComboBoxes.size()) {
            if (m_typeComboBoxes[i]->isEditable()) {
                // 如果是自定义类型，设置为Custom
                def.type = WellTestColumnType::Custom;
            } else {
                def.type = static_cast<WellTestColumnType>(m_typeComboBoxes[i]->currentIndex());
            }
        }

        def.unit = unitName;

        if (i < m_requiredChecks.size()) {
            def.isRequired = m_requiredChecks[i]->isChecked();
        }

        definitions.append(def);
    }

    return definitions;
}

// ============================================================================
// 修改后的时间转换对话框实现
// ============================================================================

TimeConversionDialog::TimeConversionDialog(const QStringList& columnNames, QWidget* parent)
    : QDialog(parent), m_columnNames(columnNames)
{
    setupUI();
    m_config.dateColumnIndex = -1;
    m_config.timeColumnIndex = -1;
    m_config.sourceTimeColumnIndex = 0;
    m_config.outputUnit = "s";
    m_config.newColumnName = "时间";
    m_config.useDateAndTime = false;
}

void TimeConversionDialog::setupUI()
{
    setWindowTitle("时间转换设置");
    setModal(true);
    resize(500, 400);

    QVBoxLayout* mainLayout = new QVBoxLayout(this);

    // 转换模式选择
    QGroupBox* modeGroup = new QGroupBox("转换模式");
    QVBoxLayout* modeLayout = new QVBoxLayout(modeGroup);

    m_dateTimeRadio = new QRadioButton("日期+时刻模式");
    m_timeOnlyRadio = new QRadioButton("仅时间模式");
    m_timeOnlyRadio->setChecked(true); // 默认选择仅时间模式

    modeLayout->addWidget(m_dateTimeRadio);
    modeLayout->addWidget(m_timeOnlyRadio);

    connect(m_dateTimeRadio, &QRadioButton::toggled, this, &TimeConversionDialog::onConversionModeChanged);
    connect(m_timeOnlyRadio, &QRadioButton::toggled, this, &TimeConversionDialog::onConversionModeChanged);

    mainLayout->addWidget(modeGroup);

    // 配置区域
    QGroupBox* configGroup = new QGroupBox("配置参数");
    QFormLayout* formLayout = new QFormLayout(configGroup);

    // 日期列选择（日期+时刻模式）
    m_dateColumnCombo = new QComboBox;
    m_dateColumnCombo->addItems(m_columnNames);
    formLayout->addRow("日期列:", m_dateColumnCombo);

    // 时刻列选择（日期+时刻模式）
    m_timeColumnCombo = new QComboBox;
    m_timeColumnCombo->addItems(m_columnNames);
    formLayout->addRow("时刻列:", m_timeColumnCombo);

    // 源时间列选择（仅时间模式）
    m_sourceColumnCombo = new QComboBox;
    m_sourceColumnCombo->addItems(m_columnNames);
    formLayout->addRow("源时间列:", m_sourceColumnCombo);

    // 新列名输入
    m_newColumnNameEdit = new QLineEdit("时间");
    formLayout->addRow("新列名:", m_newColumnNameEdit);

    // 输出单位选择
    m_outputUnitCombo = new QComboBox;
    m_outputUnitCombo->addItems({"s", "m", "h"});
    m_outputUnitCombo->setCurrentText("s");
    formLayout->addRow("输出单位:", m_outputUnitCombo);

    mainLayout->addWidget(configGroup);

    // 预览区域
    QGroupBox* previewGroup = new QGroupBox("预览");
    QVBoxLayout* previewLayout = new QVBoxLayout(previewGroup);

    m_previewButton = new QPushButton("生成预览");
    connect(m_previewButton, &QPushButton::clicked, this, &TimeConversionDialog::onPreviewClicked);
    previewLayout->addWidget(m_previewButton);

    m_previewLabel = new QLabel("点击'生成预览'查看转换效果");
    m_previewLabel->setStyleSheet("color: #6c757d; font-size: 11px; padding: 8px; border: 1px solid #e1e8ed; border-radius: 4px;");
    m_previewLabel->setWordWrap(true);
    previewLayout->addWidget(m_previewLabel);

    mainLayout->addWidget(previewGroup);

    // 按钮区域
    QHBoxLayout* buttonLayout = new QHBoxLayout;
    buttonLayout->addStretch();

    QPushButton* okBtn = new QPushButton("确定");
    connect(okBtn, &QPushButton::clicked, this, &QDialog::accept);
    buttonLayout->addWidget(okBtn);

    QPushButton* cancelBtn = new QPushButton("取消");
    connect(cancelBtn, &QPushButton::clicked, this, &QDialog::reject);
    buttonLayout->addWidget(cancelBtn);

    mainLayout->addLayout(buttonLayout);

    // 初始化UI状态
    updateUIForMode();
}

void TimeConversionDialog::onConversionModeChanged()
{
    updateUIForMode();
}

void TimeConversionDialog::updateUIForMode()
{
    bool useDateAndTime = m_dateTimeRadio->isChecked();

    // 显示/隐藏相应的控件
    m_dateColumnCombo->setEnabled(useDateAndTime);
    m_timeColumnCombo->setEnabled(useDateAndTime);
    m_sourceColumnCombo->setEnabled(!useDateAndTime);
}

void TimeConversionDialog::onPreviewClicked()
{
    // 生成预览文本
    QString unitText;
    QString unit = m_outputUnitCombo->currentText();
    if (unit == "s") {
        unitText = "秒";
    } else if (unit == "m") {
        unitText = "分钟";
    } else if (unit == "h") {
        unitText = "小时";
    }

    QString preview;

    if (m_dateTimeRadio->isChecked()) {
        // 日期+时刻模式
        preview = QString("将基于日期列 '%1' 和时刻列 '%2' 创建新列 '%3'，单位为%4。\n\n")
                      .arg(m_dateColumnCombo->currentText())
                      .arg(m_timeColumnCombo->currentText())
                      .arg(m_newColumnNameEdit->text())
                      .arg(unitText);

        preview += "转换规则：\n";
        preview += "• 第1行时间 = 0（基准时间）\n";
        preview += "• 第n行时间 = (第n行日期-第1行日期)*24 + (第n行时刻-第1行时刻)\n";
        preview += "• 如果是同一天数据，日期差为0，主要计算时刻差\n\n";

        preview += "示例（假设日期格式为 yyyy-MM-dd，时刻格式为 HH:mm:ss）：\n";
        preview += previewConversion("2006-07-18", "10:25:10", unit) + "\n";
        preview += previewConversion("2006-07-18", "10:25:15", unit) + "\n";
        preview += previewConversion("2006-07-18", "10:25:20", unit) + "\n";
    } else {
        // 仅时间模式
        preview = QString("将基于源列 '%1' 创建新列 '%2'，单位为%3。\n\n")
                      .arg(m_sourceColumnCombo->currentText())
                      .arg(m_newColumnNameEdit->text())
                      .arg(unitText);

        preview += "转换规则：\n";
        preview += "• 第1行时间 = 0（基准时间）\n";
        preview += "• 第2行时间 = 第2行原始时间 - 第1行原始时间\n";
        preview += "• 第3行时间 = 第3行原始时间 - 第1行原始时间\n";
        preview += "• 以此类推...\n\n";

        preview += "示例（假设原始时间格式为 HH:MM:SS）：\n";
        preview += previewConversion("", "10:25:10", unit) + "\n";
        preview += previewConversion("", "10:25:15", unit) + "\n";
        preview += previewConversion("", "10:25:20", unit) + "\n";
    }

    m_previewLabel->setText(preview);
}

QString TimeConversionDialog::previewConversion(const QString& sampleDateInput, const QString& sampleTimeInput, const QString& unit)
{
    static QString baseDate = "2006-07-18";
    static QString baseTime = "10:25:10";

    if (m_dateTimeRadio->isChecked()) {
        // 日期+时刻模式
        if (sampleDateInput == baseDate && sampleTimeInput == baseTime) {
            return QString("日期: %1, 时刻: %2 => 转换时间: 0 %3").arg(sampleDateInput).arg(sampleTimeInput).arg(unit);
        }

        // 简单计算示例时间差
        QDate date1 = QDate::fromString(baseDate, "yyyy-MM-dd");
        QDate date2 = QDate::fromString(sampleDateInput, "yyyy-MM-dd");
        QTime time1 = QTime::fromString(baseTime, "hh:mm:ss");
        QTime time2 = QTime::fromString(sampleTimeInput, "hh:mm:ss");

        if (date1.isValid() && date2.isValid() && time1.isValid() && time2.isValid()) {
            int daysDiff = date1.daysTo(date2);
            int timeDiffSeconds = time1.secsTo(time2);
            int totalSeconds = daysDiff * 24 * 3600 + timeDiffSeconds;

            double value = 0;
            if (unit == "s") {
                value = totalSeconds;
            } else if (unit == "m") {
                value = totalSeconds / 60.0;
            } else if (unit == "h") {
                value = totalSeconds / 3600.0;
            }

            return QString("日期: %1, 时刻: %2 => 转换时间: %3 %4")
                .arg(sampleDateInput)
                .arg(sampleTimeInput)
                .arg(QString::number(value, 'f', 3))
                .arg(unit);
        }
    } else {
        // 仅时间模式
        if (sampleTimeInput == baseTime) {
            return QString("原始时间: %1 => 转换时间: 0 %2").arg(sampleTimeInput).arg(unit);
        }

        // 简单计算示例时间差
        QTime base = QTime::fromString(baseTime, "hh:mm:ss");
        QTime current = QTime::fromString(sampleTimeInput, "hh:mm:ss");

        if (base.isValid() && current.isValid()) {
            int diffSeconds = base.secsTo(current);
            double value = 0;

            if (unit == "s") {
                value = diffSeconds;
            } else if (unit == "m") {
                value = diffSeconds / 60.0;
            } else if (unit == "h") {
                value = diffSeconds / 3600.0;
            }

            return QString("原始时间: %1 => 转换时间: %2 %3")
                .arg(sampleTimeInput)
                .arg(QString::number(value, 'f', 3))
                .arg(unit);
        }
    }

    return QString("示例数据格式错误");
}

TimeConversionConfig TimeConversionDialog::getConversionConfig() const
{
    TimeConversionConfig config;

    config.useDateAndTime = m_dateTimeRadio->isChecked();

    if (config.useDateAndTime) {
        config.dateColumnIndex = m_dateColumnCombo->currentIndex();
        config.timeColumnIndex = m_timeColumnCombo->currentIndex();
        config.sourceTimeColumnIndex = -1;
    } else {
        config.dateColumnIndex = -1;
        config.timeColumnIndex = -1;
        config.sourceTimeColumnIndex = m_sourceColumnCombo->currentIndex();
    }

    config.outputUnit = m_outputUnitCombo->currentText();
    config.newColumnName = m_newColumnNameEdit->text().trimmed();

    if (config.newColumnName.isEmpty()) {
        config.newColumnName = "时间";
    }

    return config;
}
// ============================================================================
// 数据读取配置对话框实现
// ============================================================================

DataLoadConfigDialog::DataLoadConfigDialog(const QString& filePath, QWidget* parent)
    : QDialog(parent), m_filePath(filePath)
{
    setupUI();

    // 设置默认配置
    m_config.startRow = 1;
    m_config.hasHeader = true;
    m_config.encoding = detectEncoding(filePath);
    m_config.separator = detectSeparator(filePath);

    // 设置UI初始值
    m_startRowSpin->setValue(m_config.startRow);
    m_hasHeaderCheck->setChecked(m_config.hasHeader);

    int encodingIndex = m_encodingCombo->findText(m_config.encoding);
    if (encodingIndex >= 0) {
        m_encodingCombo->setCurrentIndex(encodingIndex);
    }

    int separatorIndex = m_separatorCombo->findData(m_config.separator);
    if (separatorIndex >= 0) {
        m_separatorCombo->setCurrentIndex(separatorIndex);
    }

    // 加载文件预览
    loadFilePreview();
}

void DataLoadConfigDialog::setupUI()
{
    setWindowTitle("数据读取配置");
    setModal(true);
    resize(600, 500);

    QVBoxLayout* mainLayout = new QVBoxLayout(this);

    // 说明标签
    QLabel* infoLabel = new QLabel("请配置数据读取参数，预览文件内容以确认设置：");
    infoLabel->setStyleSheet("font-size: 14px; font-weight: bold; color: #2c3e50; margin: 10px;");
    mainLayout->addWidget(infoLabel);

    // 配置区域
    QGroupBox* configGroup = new QGroupBox("读取配置");
    QFormLayout* configLayout = new QFormLayout(configGroup);

    // 起始行
    m_startRowSpin = new QSpinBox;
    m_startRowSpin->setRange(1, 1000);
    m_startRowSpin->setValue(1);
    m_startRowSpin->setToolTip("指定从第几行开始读取数据（包含该行）");
    configLayout->addRow("起始行:", m_startRowSpin);

    // 是否有表头
    m_hasHeaderCheck = new QCheckBox("文件包含表头");
    m_hasHeaderCheck->setChecked(true);
    m_hasHeaderCheck->setToolTip("勾选此项将使用第一行作为列标题");
    connect(m_hasHeaderCheck, &QCheckBox::toggled, this, &DataLoadConfigDialog::onHasHeaderChanged);
    configLayout->addRow("", m_hasHeaderCheck);

    // 编码格式
    m_encodingCombo = new QComboBox;
    m_encodingCombo->addItems({"UTF-8", "GBK", "GB2312", "ASCII"});
    configLayout->addRow("编码格式:", m_encodingCombo);

    // 分隔符
    m_separatorCombo = new QComboBox;
    m_separatorCombo->addItem("逗号 (,)", ",");
    m_separatorCombo->addItem("制表符 (Tab)", "\t");
    m_separatorCombo->addItem("分号 (;)", ";");
    m_separatorCombo->addItem("竖线 (|)", "|");
    m_separatorCombo->addItem("空格", " ");
    configLayout->addRow("分隔符:", m_separatorCombo);

    mainLayout->addWidget(configGroup);

    // 预览区域
    QGroupBox* previewGroup = new QGroupBox("文件预览");
    QVBoxLayout* previewLayout = new QVBoxLayout(previewGroup);

    m_previewButton = new QPushButton("刷新预览");
    connect(m_previewButton, &QPushButton::clicked, this, &DataLoadConfigDialog::onPreviewClicked);
    previewLayout->addWidget(m_previewButton);

    m_previewText = new QTextEdit;
    m_previewText->setReadOnly(true);
    m_previewText->setMaximumHeight(200);
    m_previewText->setStyleSheet("font-family: 'Consolas', 'Monaco', monospace; font-size: 10px;");
    previewLayout->addWidget(m_previewText);

    mainLayout->addWidget(previewGroup);

    // 按钮区域
    QHBoxLayout* buttonLayout = new QHBoxLayout;
    buttonLayout->addStretch();

    QPushButton* okBtn = new QPushButton("确定");
    okBtn->setStyleSheet("QPushButton { background-color: #28a745; color: white; border: none; border-radius: 4px; padding: 8px 16px; }");
    connect(okBtn, &QPushButton::clicked, this, &QDialog::accept);
    buttonLayout->addWidget(okBtn);

    QPushButton* cancelBtn = new QPushButton("取消");
    cancelBtn->setStyleSheet("QPushButton { background-color: #fd7e14; color: white; border: none; border-radius: 4px; padding: 8px 16px; }");
    connect(cancelBtn, &QPushButton::clicked, this, &QDialog::reject);
    buttonLayout->addWidget(cancelBtn);

    mainLayout->addLayout(buttonLayout);
}

void DataLoadConfigDialog::onPreviewClicked()
{
    loadFilePreview();
}

void DataLoadConfigDialog::onHasHeaderChanged(bool hasHeader)
{
    // 当表头选项改变时，自动调整起始行
    if (hasHeader && m_startRowSpin->value() == 1) {
        m_startRowSpin->setValue(2);
    } else if (!hasHeader && m_startRowSpin->value() == 2) {
        m_startRowSpin->setValue(1);
    }
    loadFilePreview();
}

void DataLoadConfigDialog::loadFilePreview()
{
    QFile file(m_filePath);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        m_previewText->setText("无法读取文件");
        return;
    }

    QTextStream in(&file);
    QString encoding = m_encodingCombo->currentText();

#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    if (encoding != "UTF-8") {
        QTextCodec* codec = QTextCodec::codecForName(encoding.toLocal8Bit());
        if (codec) {
            in.setCodec(codec);
        }
    }
#else
    if (encoding == "GBK" || encoding == "GB2312") {
        in.setEncoding(QStringConverter::System);
    }
#endif

    QString previewText = QString("文件: %1\n").arg(QFileInfo(m_filePath).fileName());
    previewText += QString("编码: %1\n").arg(encoding);
    previewText += QString("分隔符: %1\n").arg(m_separatorCombo->currentText());
    previewText += QString("起始行: %1\n").arg(m_startRowSpin->value());
    previewText += QString("包含表头: %1\n").arg(m_hasHeaderCheck->isChecked() ? "是" : "否");
    previewText += QString("-").repeated(50) + "\n";

    // 读取前20行进行预览
    QStringList lines;
    int lineCount = 0;
    while (!in.atEnd() && lineCount < 20) {
        lines.append(in.readLine());
        lineCount++;
    }

    int startRow = m_startRowSpin->value() - 1; // 转换为0基索引
    QString separator = m_separatorCombo->currentData().toString();

    for (int i = 0; i < lines.size(); ++i) {
        QString prefix;
        if (i < startRow) {
            prefix = QString("[跳过] 第%1行: ").arg(i + 1);
        } else if (i == startRow && m_hasHeaderCheck->isChecked()) {
            prefix = QString("[表头] 第%1行: ").arg(i + 1);
        } else {
            prefix = QString("[数据] 第%1行: ").arg(i + 1);
        }

        QString line = lines[i];
        if (separator != " ") {
            // 将分隔符替换为 | 以便更好地显示列分隔
            line = line.replace(separator, " | ");
        }

        previewText += prefix + line + "\n";
    }

    file.close();
    m_previewText->setText(previewText);
}

QString DataLoadConfigDialog::detectEncoding(const QString& filePath)
{
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly)) {
        return "UTF-8";
    }

    QByteArray data = file.read(1024); // 读取前1KB进行检测
    file.close();

    // 简单的编码检测
    if (data.contains('\0')) {
        return "UTF-8"; // 可能包含二进制数据
    }

    // 检测是否包含中文
    QString testUtf8 = QString::fromUtf8(data);
    QString testGbk;

#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    QTextCodec* gbkCodec = QTextCodec::codecForName("GBK");
    if (gbkCodec) {
        testGbk = gbkCodec->toUnicode(data);
    }
#endif

    // 如果包含明显的中文字符，尝试GBK
    if (data.size() > testUtf8.toUtf8().size()) {
        return "GBK";
    }

    return "UTF-8";
}

QString DataLoadConfigDialog::detectSeparator(const QString& filePath)
{
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        return ",";
    }

    // 读取前几行检测分隔符
    QTextStream in(&file);
    QString sampleText;
    int lineCount = 0;
    while (!in.atEnd() && lineCount < 5) {
        sampleText += in.readLine() + "\n";
        lineCount++;
    }
    file.close();

    // 统计各种分隔符的出现频率
    QMap<QString, int> separatorCounts;
    separatorCounts[","] = sampleText.count(',');
    separatorCounts["\t"] = sampleText.count('\t');
    separatorCounts[";"] = sampleText.count(';');
    separatorCounts["|"] = sampleText.count('|');
    separatorCounts[" "] = sampleText.count(' ');

    // 找到出现频率最高的分隔符
    QString bestSeparator = ",";
    int maxCount = 0;
    for (auto it = separatorCounts.begin(); it != separatorCounts.end(); ++it) {
        if (it.value() > maxCount) {
            maxCount = it.value();
            bestSeparator = it.key();
        }
    }

    return bestSeparator;
}

DataLoadConfigDialog::LoadConfig DataLoadConfigDialog::getLoadConfig() const
{
    LoadConfig config;
    config.startRow = m_startRowSpin->value();
    config.hasHeader = m_hasHeaderCheck->isChecked();
    config.encoding = m_encodingCombo->currentText();
    config.separator = m_separatorCombo->currentData().toString();
    return config;
}

// ============================================================================
// 其他对话框实现（简化版）
// ============================================================================

DataCleaningDialog::DataCleaningDialog(QWidget* parent) : QDialog(parent)
{
    setupUI();
}

void DataCleaningDialog::setupUI()
{
    setWindowTitle("数据清理选项");
    setModal(true);
    resize(380, 280);

    QVBoxLayout* mainLayout = new QVBoxLayout(this);

    m_removeEmptyRowsCheck = new QCheckBox("删除空行");
    m_removeEmptyRowsCheck->setChecked(true);
    mainLayout->addWidget(m_removeEmptyRowsCheck);

    m_removeEmptyColumnsCheck = new QCheckBox("删除空列");
    mainLayout->addWidget(m_removeEmptyColumnsCheck);

    m_removeDuplicatesCheck = new QCheckBox("删除重复行");
    m_removeDuplicatesCheck->setChecked(true);
    mainLayout->addWidget(m_removeDuplicatesCheck);

    m_fillMissingValuesCheck = new QCheckBox("填充缺失值");
    mainLayout->addWidget(m_fillMissingValuesCheck);

    QHBoxLayout* fillLayout = new QHBoxLayout;
    fillLayout->addWidget(new QLabel("填充方法:"));
    m_fillMethodCombo = new QComboBox;
    m_fillMethodCombo->addItems({"零值", "线性插值", "平均值", "前值填充"});
    m_fillMethodCombo->setCurrentIndex(1);
    fillLayout->addWidget(m_fillMethodCombo);
    mainLayout->addLayout(fillLayout);

    m_removeOutliersCheck = new QCheckBox("删除异常值");
    mainLayout->addWidget(m_removeOutliersCheck);

    QHBoxLayout* outlierLayout = new QHBoxLayout;
    outlierLayout->addWidget(new QLabel("异常值阈值:"));
    m_outlierThresholdSpin = new QSpinBox;
    m_outlierThresholdSpin->setRange(1, 5);
    m_outlierThresholdSpin->setValue(2);
    m_outlierThresholdSpin->setSuffix(" 倍标准差");
    outlierLayout->addWidget(m_outlierThresholdSpin);
    mainLayout->addLayout(outlierLayout);

    m_standardizeFormatCheck = new QCheckBox("标准化数据格式");
    mainLayout->addWidget(m_standardizeFormatCheck);

    mainLayout->addStretch();

    QHBoxLayout* buttonLayout = new QHBoxLayout;
    buttonLayout->addStretch();

    QPushButton* okBtn = new QPushButton("执行清理");
    connect(okBtn, &QPushButton::clicked, this, &QDialog::accept);
    buttonLayout->addWidget(okBtn);

    QPushButton* cancelBtn = new QPushButton("取消");
    connect(cancelBtn, &QPushButton::clicked, this, &QDialog::reject);
    buttonLayout->addWidget(cancelBtn);

    mainLayout->addLayout(buttonLayout);
}

DataCleaningDialog::CleaningOptions DataCleaningDialog::getCleaningOptions() const
{
    CleaningOptions options;
    options.removeEmptyRows = m_removeEmptyRowsCheck->isChecked();
    options.removeEmptyColumns = m_removeEmptyColumnsCheck->isChecked();
    options.removeDuplicates = m_removeDuplicatesCheck->isChecked();
    options.fillMissingValues = m_fillMissingValuesCheck->isChecked();
    options.removeOutliers = m_removeOutliersCheck->isChecked();
    options.standardizeFormat = m_standardizeFormatCheck->isChecked();

    QStringList fillMethods = {"zero", "interpolation", "average", "forward"};
    options.fillMethod = fillMethods[m_fillMethodCombo->currentIndex()];
    options.outlierThreshold = m_outlierThresholdSpin->value();

    return options;
}

// 动画进度对话框简化实现
AnimatedProgressDialog::AnimatedProgressDialog(const QString& title, const QString& message, QWidget* parent)
    : QDialog(parent)
{
    setWindowTitle(title);
    setModal(true);
    setFixedSize(320, 100);
    setWindowFlags(Qt::Dialog | Qt::CustomizeWindowHint | Qt::WindowTitleHint);

    setupUI();
    setMessage(message);
}

void AnimatedProgressDialog::setupUI()
{
    QVBoxLayout* mainLayout = new QVBoxLayout(this);

    m_messageLabel = new QLabel;
    m_messageLabel->setStyleSheet("font-size: 13px; color: #2c3e50;");
    m_messageLabel->setWordWrap(true);
    mainLayout->addWidget(m_messageLabel);

    m_progressBar = new QProgressBar;
    m_progressBar->setStyleSheet(R"(
        QProgressBar {
            border: 1px solid #e1e8ed;
            border-radius: 6px;
            text-align: center;
            background-color: #f8f9fa;
        }
        QProgressBar::chunk {
            background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                        stop:0 #4a90e2, stop:1 #357abd);
            border-radius: 5px;
        }
    )");
    m_progressBar->setMinimum(0);
    m_progressBar->setMaximum(100);
    mainLayout->addWidget(m_progressBar);
}

void AnimatedProgressDialog::setupAnimation() {}

void AnimatedProgressDialog::setProgress(int value)
{
    m_progressBar->setValue(value);
}

void AnimatedProgressDialog::setMessage(const QString& message)
{
    m_messageLabel->setText(message);
}

void AnimatedProgressDialog::setMaximum(int maximum)
{
    m_progressBar->setMaximum(maximum);
}

void AnimatedProgressDialog::closeEvent(QCloseEvent* event)
{
    event->ignore();
}

// ============================================================================
// DataEditorWidget 主类实现
// ============================================================================

DataEditorWidget::DataEditorWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::DataEditorWidget),
    m_dataModel(nullptr),
    m_proxyModel(nullptr),
    m_undoStack(nullptr),
    m_dataModified(false),
    m_searchTimer(nullptr),
    m_progressDialog(nullptr),
    m_largeFileMode(false),
    m_maxDisplayRows(10000),
    m_contextMenu(nullptr),
    m_addRowAboveAction(nullptr),
    m_addRowBelowAction(nullptr),
    m_deleteRowsAction(nullptr),
    m_addColumnLeftAction(nullptr),
    m_addColumnRightAction(nullptr),
    m_deleteColumnsAction(nullptr),
    m_pressureDerivativeCalculator(nullptr)
{
    ui->setupUi(this);
    init();
}

DataEditorWidget::~DataEditorWidget()
{
    delete ui;
    if (m_dataModel) {
        delete m_dataModel;
    }
    if (m_proxyModel) {
        delete m_proxyModel;
    }
    if (m_undoStack) {
        delete m_undoStack;
    }
    if (m_progressDialog) {
        delete m_progressDialog;
    }
    if (m_contextMenu) {
        delete m_contextMenu;
    }
    if (m_pressureDerivativeCalculator) {
        delete m_pressureDerivativeCalculator;
    }
}

void DataEditorWidget::init()
{
    setupModels();
    setupUI();
    setupConnections();
    setupContextMenu();
    setupPressureDerivativeCalculator();

    // 初始化搜索定时器
    m_searchTimer = new QTimer(this);
    m_searchTimer->setSingleShot(true);
    m_searchTimer->setInterval(300);
    connect(m_searchTimer, &QTimer::timeout, this, &DataEditorWidget::onSearchData);

    // 初始状态
    setButtonsEnabled(false);
    updateStatus("就绪", "success");
    updateDataInfo();
}

void DataEditorWidget::setupModels()
{
    // 创建数据模型
    m_dataModel = new QStandardItemModel(this);

    // 创建代理模型用于搜索和筛选
    m_proxyModel = new QSortFilterProxyModel(this);
    m_proxyModel->setSourceModel(m_dataModel);
    m_proxyModel->setFilterCaseSensitivity(Qt::CaseInsensitive);
    m_proxyModel->setFilterKeyColumn(-1);

    // 设置表格视图的模型
    ui->dataTableView->setModel(m_proxyModel);

    // 优化表格属性
    ui->dataTableView->setAlternatingRowColors(true);
    ui->dataTableView->horizontalHeader()->setStretchLastSection(false);
    ui->dataTableView->horizontalHeader()->setSectionResizeMode(QHeaderView::Interactive);
    ui->dataTableView->verticalHeader()->setSectionResizeMode(QHeaderView::Interactive);
    ui->dataTableView->setSortingEnabled(false);
    ui->dataTableView->verticalHeader()->setVisible(true);

    // 设置合适的行高和列宽
    ui->dataTableView->verticalHeader()->setDefaultSectionSize(24);
    ui->dataTableView->verticalHeader()->setMinimumSectionSize(20);
    ui->dataTableView->horizontalHeader()->setDefaultSectionSize(100);
    ui->dataTableView->horizontalHeader()->setMinimumSectionSize(60);

    // 设置表格选择模式
    ui->dataTableView->setSelectionBehavior(QAbstractItemView::SelectItems);
    ui->dataTableView->setSelectionMode(QAbstractItemView::ExtendedSelection);

    // 创建撤销栈
    m_undoStack = new QUndoStack(this);
}

void DataEditorWidget::setupUI()
{
    // 设置初始状态
    ui->filePathLineEdit->setReadOnly(true);
    ui->filePathLineEdit->setPlaceholderText("📁 未选择文件");

    // 优化表格显示设置
    ui->dataTableView->setShowGrid(true);
    ui->dataTableView->setGridStyle(Qt::SolidLine);
    ui->dataTableView->setWordWrap(false);

    // 设置表格字体 - 调整为更小的字体以确保序号显示完全
    QFont tableFont = ui->dataTableView->font();
    tableFont.setPointSize(10);  // 稍微增大字体
    ui->dataTableView->setFont(tableFont);

    // 设置行标题字体与表格数据一致
    QFont headerFont = ui->dataTableView->verticalHeader()->font();
    headerFont.setPointSize(10);  // 与数据字体大小一致
    ui->dataTableView->verticalHeader()->setFont(headerFont);
}

void DataEditorWidget::setupConnections()
{
    // 文件操作按钮
    connect(ui->btnOpenFile, &QPushButton::clicked, this, &DataEditorWidget::onOpenFile);
    connect(ui->btnSave, &QPushButton::clicked, this, &DataEditorWidget::onSave);
    connect(ui->btnExport, &QPushButton::clicked, this, &DataEditorWidget::onExport);

    // 数据处理按钮
    connect(ui->btnDefineColumns, &QPushButton::clicked, this, &DataEditorWidget::onDefineColumns);
    connect(ui->btnTimeConvert, &QPushButton::clicked, this, &DataEditorWidget::onTimeConvert);
    connect(ui->btnPressureDropCalc, &QPushButton::clicked, this, &DataEditorWidget::onPressureDropCalc);
    connect(ui->btnPressureDerivativeCalc, &QPushButton::clicked, this, &DataEditorWidget::onPressureDerivativeCalc);
    connect(ui->btnDataClean, &QPushButton::clicked, this, &DataEditorWidget::onDataClean);
    connect(ui->btnDataStatistics, &QPushButton::clicked, this, &DataEditorWidget::onDataStatistics);

    // 搜索功能
    connect(ui->searchLineEdit, &QLineEdit::textChanged, this, &DataEditorWidget::onSearchTextChanged);

    // 模型数据变化
    connect(m_dataModel, &QStandardItemModel::itemChanged, this, &DataEditorWidget::onCellDataChanged);
    connect(m_dataModel, &QStandardItemModel::dataChanged, this, &DataEditorWidget::onModelDataChanged);

    // 右键菜单连接
    connect(ui->dataTableView, &QTableView::customContextMenuRequested,
            this, &DataEditorWidget::onTableContextMenuRequested);
}

void DataEditorWidget::setupContextMenu()
{
    // 创建右键菜单
    m_contextMenu = new QMenu(this);

    // 创建动作
    m_addRowAboveAction = new QAction("在上方插入行", this);
    m_addRowBelowAction = new QAction("在下方插入行", this);
    m_deleteRowsAction = new QAction("删除选中行", this);

    QAction* separator1 = new QAction(this);
    separator1->setSeparator(true);

    m_addColumnLeftAction = new QAction("在左侧插入列", this);
    m_addColumnRightAction = new QAction("在右侧插入列", this);
    m_deleteColumnsAction = new QAction("删除选中列", this);

    // 添加动作到菜单
    m_contextMenu->addAction(m_addRowAboveAction);
    m_contextMenu->addAction(m_addRowBelowAction);
    m_contextMenu->addAction(m_deleteRowsAction);
    m_contextMenu->addAction(separator1);
    m_contextMenu->addAction(m_addColumnLeftAction);
    m_contextMenu->addAction(m_addColumnRightAction);
    m_contextMenu->addAction(m_deleteColumnsAction);

    // 连接动作槽函数
    connect(m_addRowAboveAction, &QAction::triggered, this, &DataEditorWidget::onAddRowAbove);
    connect(m_addRowBelowAction, &QAction::triggered, this, &DataEditorWidget::onAddRowBelow);
    connect(m_deleteRowsAction, &QAction::triggered, this, &DataEditorWidget::onDeleteSelectedRows);
    connect(m_addColumnLeftAction, &QAction::triggered, this, &DataEditorWidget::onAddColumnLeft);
    connect(m_addColumnRightAction, &QAction::triggered, this, &DataEditorWidget::onAddColumnRight);
    connect(m_deleteColumnsAction, &QAction::triggered, this, &DataEditorWidget::onDeleteSelectedColumns);

    // 设置动作样式
    QString actionStyle = R"(
        QMenu {
            background-color: white;
            border: 1px solid #e1e8ed;
            border-radius: 6px;
            padding: 4px;
        }
        QMenu::item {
            background-color: transparent;
            padding: 8px 16px;
            color: #2c3e50;
            border-radius: 3px;
        }
        QMenu::item:selected {
            background-color: #f0f8ff;
            color: #2c3e50;
        }
        QMenu::separator {
            height: 1px;
            background-color: #e1e8ed;
            margin: 4px 8px;
        }
    )";
    m_contextMenu->setStyleSheet(actionStyle);
}

// ============================================================================
// 修改后的时间转换功能实现
// ============================================================================

void DataEditorWidget::onTimeConvert()
{
    if (!hasData()) {
        showStyledMessageBox("时间转换", "请先加载数据文件", QMessageBox::Information);
        return;
    }

    QStringList columnNames;
    for (int i = 0; i < m_dataModel->columnCount(); ++i) {
        columnNames.append(m_dataModel->headerData(i, Qt::Horizontal).toString());
    }

    TimeConversionDialog dialog(columnNames, this);
    if (dialog.exec() == QDialog::Accepted) {
        TimeConversionConfig config = dialog.getConversionConfig();

        showAnimatedProgress("时间转换", "正在转换时间数据...");

        TimeConversionResult result;
        try {
            result = convertTimeColumn(config);
        } catch (const std::exception& e) {
            result.success = false;
            result.errorMessage = QString("转换过程中发生错误: %1").arg(e.what());
        } catch (...) {
            result.success = false;
            result.errorMessage = "转换过程中发生未知错误";
        }

        hideAnimatedProgress();

        if (result.success) {
            updateStatus(QString("时间转换完成 - 已添加列: %1").arg(result.columnName), "success");
            m_dataModified = true;

            // 安全地发射信号和更新界面
            try {
                emit timeConversionCompleted(result);

                // 延迟发射数据变化信号，避免立即执行可能导致的问题
                QTimer::singleShot(200, [this]() {
                    try {
                        emit dataChanged();
                    } catch (...) {
                        qDebug() << "发射dataChanged信号时出错";
                    }
                });

                showStyledMessageBox("时间转换完成",
                                     QString("时间转换成功完成！\n"
                                             "新增列：%1\n"
                                             "处理行数：%2")
                                         .arg(result.columnName)
                                         .arg(result.processedRows),
                                     QMessageBox::Information);
            } catch (...) {
                qDebug() << "时间转换完成后处理信号时出错";
            }
        } else {
            updateStatus("时间转换失败", "error");
            showStyledMessageBox("时间转换失败", result.errorMessage, QMessageBox::Warning);
        }
    }
}

TimeConversionResult DataEditorWidget::convertTimeColumn(const TimeConversionConfig& config)
{
    TimeConversionResult result;
    result.success = false;
    result.addedColumnIndex = -1;
    result.processedRows = 0;

    try {
        if (!m_dataModel) {
            result.errorMessage = "数据模型不存在";
            return result;
        }

        // 创建新列名，包含单位信息
        QString unitText;
        if (config.outputUnit == "s") {
            unitText = "s";
        } else if (config.outputUnit == "m") {
            unitText = "min";
        } else if (config.outputUnit == "h") {
            unitText = "h";
        }

        QString newColumnName = QString("%1\\%2").arg(config.newColumnName).arg(unitText);

        int newColumnIndex;

        if (config.useDateAndTime) {
            // 日期+时刻模式
            if (config.dateColumnIndex < 0 || config.dateColumnIndex >= m_dataModel->columnCount() ||
                config.timeColumnIndex < 0 || config.timeColumnIndex >= m_dataModel->columnCount()) {
                result.errorMessage = "日期或时刻列索引无效";
                return result;
            }

            // 在时刻列后面插入新列
            newColumnIndex = qMax(config.dateColumnIndex, config.timeColumnIndex) + 1;
            m_dataModel->insertColumn(newColumnIndex);

            // 设置列标题
            QStandardItem* headerItem = new QStandardItem(newColumnName);
            m_dataModel->setHorizontalHeaderItem(newColumnIndex, headerItem);

            // 获取基准日期和时刻（第一行的数据）
            QDate baseDate;
            QTime baseTime;
            bool baseSet = false;

            // 找到第一个有效的日期和时刻
            for (int row = 0; row < m_dataModel->rowCount(); ++row) {
                QStandardItem* dateItem = m_dataModel->item(row, config.dateColumnIndex);
                QStandardItem* timeItem = m_dataModel->item(row, config.timeColumnIndex);

                if (dateItem && timeItem) {
                    QString dateStr = dateItem->text().trimmed();
                    QString timeStr = timeItem->text().trimmed();

                    QDate parsedDate = parseDateString(dateStr);
                    QTime parsedTime = parseTimeString(timeStr);

                    if (parsedDate.isValid() && parsedTime.isValid()) {
                        baseDate = parsedDate;
                        baseTime = parsedTime;
                        baseSet = true;
                        break;
                    }
                }
            }

            if (!baseSet) {
                result.errorMessage = "未找到有效的日期和时刻数据";
                if (newColumnIndex < m_dataModel->columnCount()) {
                    m_dataModel->removeColumn(newColumnIndex);
                }
                return result;
            }

            // 计算每行的相对时间
            for (int row = 0; row < m_dataModel->rowCount(); ++row) {
                QString convertedValue;

                QStandardItem* dateItem = m_dataModel->item(row, config.dateColumnIndex);
                QStandardItem* timeItem = m_dataModel->item(row, config.timeColumnIndex);

                if (dateItem && timeItem) {
                    QString dateStr = dateItem->text().trimmed();
                    QString timeStr = timeItem->text().trimmed();

                    QDate currentDate = parseDateString(dateStr);
                    QTime currentTime = parseTimeString(timeStr);

                    if (currentDate.isValid() && currentTime.isValid()) {
                        if (row == 0) {
                            // 第一行时间为0
                            convertedValue = "0.000";
                        } else {
                            // 计算时间差：(当前日期-基准日期)*24 + (当前时刻-基准时刻)
                            QDateTime baseDateTime = combineDateAndTime(baseDate, baseTime);
                            QDateTime currentDateTime = combineDateAndTime(currentDate, currentTime);

                            double timeDiff = calculateDateTimeDifference(baseDateTime, currentDateTime, config.outputUnit);
                            convertedValue = QString::number(timeDiff, 'f', 3);
                        }
                        result.processedRows++;
                    } else {
                        // 无效数据设为空
                        convertedValue = "";
                    }
                } else {
                    convertedValue = "";
                }

                // 创建新的数据项
                QStandardItem* newItem = new QStandardItem(convertedValue);
                if (newItem) {
                    newItem->setForeground(QBrush(QColor("#2c3e50")));
                    m_dataModel->setItem(row, newColumnIndex, newItem);
                }
            }

        } else {
            // 仅时间模式（原有逻辑）
            if (config.sourceTimeColumnIndex < 0 || config.sourceTimeColumnIndex >= m_dataModel->columnCount()) {
                result.errorMessage = "源时间列索引无效";
                return result;
            }

            // 在源列后面插入新列
            newColumnIndex = config.sourceTimeColumnIndex + 1;
            m_dataModel->insertColumn(newColumnIndex);

            // 设置列标题
            QStandardItem* headerItem = new QStandardItem(newColumnName);
            m_dataModel->setHorizontalHeaderItem(newColumnIndex, headerItem);

            // 获取源列的所有时间数据
            QList<QTime> timeValues;
            QTime baseTime; // 基准时间（第一行的时间）
            bool baseTimeSet = false;

            // 首先解析所有时间数据
            for (int row = 0; row < m_dataModel->rowCount(); ++row) {
                QStandardItem* sourceItem = m_dataModel->item(row, config.sourceTimeColumnIndex);
                if (sourceItem) {
                    QString timeStr = sourceItem->text().trimmed();
                    QTime parsedTime = parseTimeString(timeStr);

                    if (parsedTime.isValid()) {
                        timeValues.append(parsedTime);

                        // 设置基准时间（第一个有效时间）
                        if (!baseTimeSet) {
                            baseTime = parsedTime;
                            baseTimeSet = true;
                        }
                    } else {
                        timeValues.append(QTime()); // 添加无效时间占位
                    }
                } else {
                    timeValues.append(QTime()); // 添加无效时间占位
                }
            }

            if (!baseTimeSet) {
                result.errorMessage = "未找到有效的时间数据";
                // 删除已创建的列
                if (newColumnIndex < m_dataModel->columnCount()) {
                    m_dataModel->removeColumn(newColumnIndex);
                }
                return result;
            }

            // 计算相对时间并填充新列
            for (int row = 0; row < m_dataModel->rowCount(); ++row) {
                QString convertedValue;

                if (row < timeValues.size() && timeValues[row].isValid()) {
                    if (row == 0) {
                        // 第一行时间为0
                        convertedValue = "0.000";
                    } else {
                        // 计算与基准时间的差值
                        double timeDiff = calculateTimeDifference(baseTime, timeValues[row], config.outputUnit);
                        convertedValue = QString::number(timeDiff, 'f', 3);
                    }
                    result.processedRows++;
                } else {
                    // 无效时间设为空
                    convertedValue = "";
                }

                // 创建新的数据项
                QStandardItem* newItem = new QStandardItem(convertedValue);
                if (newItem) {
                    newItem->setForeground(QBrush(QColor("#2c3e50")));
                    m_dataModel->setItem(row, newColumnIndex, newItem);
                }
            }
        }

        // 安全地添加列定义
        try {
            ColumnDefinition newColumnDef;
            newColumnDef.name = newColumnName;
            newColumnDef.type = WellTestColumnType::Time;
            newColumnDef.unit = unitText;
            newColumnDef.description = "相对时间";
            newColumnDef.isRequired = false;
            newColumnDef.minValue = 0;
            newColumnDef.maxValue = 999999;
            newColumnDef.decimalPlaces = 3;

            // 插入到列定义列表中
            if (newColumnIndex <= m_columnDefinitions.size()) {
                m_columnDefinitions.insert(newColumnIndex, newColumnDef);
            } else {
                m_columnDefinitions.append(newColumnDef);
            }
        } catch (...) {
            qDebug() << "添加列定义时出错，但转换仍然成功";
        }

        result.success = true;
        result.addedColumnIndex = newColumnIndex;
        result.columnName = newColumnName;

        // 安全地优化列宽
        try {
            QTimer::singleShot(100, [this]() {
                try {
                    optimizeColumnWidths();
                } catch (...) {
                    qDebug() << "优化列宽时出错";
                }
            });
        } catch (...) {
            qDebug() << "设置优化列宽定时器时出错";
        }

    } catch (const std::exception& e) {
        result.success = false;
        result.errorMessage = QString("转换过程中发生异常: %1").arg(e.what());
    } catch (...) {
        result.success = false;
        result.errorMessage = "转换过程中发生未知错误";
    }

    return result;
}

// ============================================================================
// 新增的日期和时间解析方法
// ============================================================================

QDate DataEditorWidget::parseDateString(const QString& dateStr) const
{
    if (dateStr.isEmpty()) {
        return QDate();
    }

    // 支持的日期格式
    QStringList dateFormats = {
        "yyyy-MM-dd",
        "yyyy/MM/dd",
        "yyyy-M-d",
        "yyyy/M/d",
        "dd/MM/yyyy",
        "dd-MM-yyyy",
        "MM/dd/yyyy",
        "MM-dd-yyyy",
        "d/M/yyyy",
        "d-M-yyyy"
    };

    for (const QString& format : dateFormats) {
        QDate date = QDate::fromString(dateStr, format);
        if (date.isValid()) {
            return date;
        }
    }

    return QDate(); // 返回无效日期
}

QDateTime DataEditorWidget::combineDateAndTime(const QDate& date, const QTime& time) const
{
    if (!date.isValid() || !time.isValid()) {
        return QDateTime();
    }
    return QDateTime(date, time);
}

double DataEditorWidget::calculateDateTimeDifference(const QDateTime& baseDateTime, const QDateTime& currentDateTime, const QString& unit) const
{
    if (!baseDateTime.isValid() || !currentDateTime.isValid()) {
        return 0.0;
    }

    // 计算秒差
    qint64 diffSeconds = baseDateTime.secsTo(currentDateTime);

    return convertTimeToUnit(static_cast<double>(diffSeconds), unit);
}

bool DataEditorWidget::isValidDateFormat(const QString& dateStr) const
{
    return parseDateString(dateStr).isValid();
}

// ============================================================================
// 更新默认列定义方法，增加日期和时刻类型
// ============================================================================

ColumnDefinition DataEditorWidget::getDefaultColumnDefinition(const QString& columnName)
{
    ColumnDefinition def;
    def.name = columnName;
    def.type = WellTestColumnType::Custom;
    def.unit = "";
    def.description = "";
    def.isRequired = false;
    def.minValue = -999999;
    def.maxValue = 999999;
    def.decimalPlaces = 3;

    // 智能识别列类型 - 增加日期、时刻识别
    QString lowerName = columnName.toLower();

    if (lowerName.contains("序号") || lowerName.contains("编号") || lowerName.contains("number") || lowerName == "no" || lowerName == "id") {
        def.type = WellTestColumnType::SerialNumber;
        def.unit = "";
        def.description = "序号";
        def.minValue = 1;
        def.maxValue = 99999;
        def.decimalPlaces = 0;
    } else if (lowerName.contains("日期") || lowerName.contains("date") || lowerName.contains("年月日")) {
        def.type = WellTestColumnType::Date;
        def.unit = "yyyy-MM-dd";
        def.description = "日期";
        def.minValue = 0;
        def.maxValue = 0;
        def.decimalPlaces = 0;
    } else if (lowerName.contains("时刻") || lowerName.contains("时分秒") || lowerName.contains("timeofday") || lowerName.contains("clock")) {
        def.type = WellTestColumnType::TimeOfDay;
        def.unit = "hh:mm:ss";
        def.description = "时刻";
        def.minValue = 0;
        def.maxValue = 0;
        def.decimalPlaces = 0;
    } else if (lowerName.contains("time") || lowerName.contains("时间") || lowerName == "t") {
        def.type = WellTestColumnType::Time;
        def.unit = "h";
        def.description = "测试时间";
        def.minValue = 0;
        def.maxValue = 10000;
    } else if (lowerName.contains("pressure") || lowerName.contains("压力") || lowerName == "p") {
        def.type = WellTestColumnType::Pressure;
        def.unit = "MPa";
        def.description = "压力数据";
        def.minValue = 0;
        def.maxValue = 100;
    } else if (lowerName.contains("temp") || lowerName.contains("温度")) {
        def.type = WellTestColumnType::Temperature;
        def.unit = "°C";
        def.description = "温度数据";
        def.minValue = -50;
        def.maxValue = 200;
    } else if (lowerName.contains("flow") || lowerName.contains("流量") || lowerName == "q") {
        def.type = WellTestColumnType::FlowRate;
        def.unit = "m³/d";
        def.description = "流量数据";
        def.minValue = 0;
        def.maxValue = 10000;
    }

    return def;
}

// ============================================================================
// 时间转换相关辅助方法实现
// ============================================================================

QTime DataEditorWidget::parseTimeString(const QString& timeStr) const
{
    if (timeStr.isEmpty()) {
        return QTime();
    }

    // 支持的时间格式
    QStringList formats = {
        "hh:mm:ss",
        "h:mm:ss",
        "hh:mm:ss.zzz",
        "h:mm:ss.zzz",
        "mm:ss",
        "m:ss"
    };

    for (const QString& format : formats) {
        QTime time = QTime::fromString(timeStr, format);
        if (time.isValid()) {
            return time;
        }
    }

    return QTime(); // 返回无效时间
}

double DataEditorWidget::calculateTimeDifference(const QTime& baseTime, const QTime& currentTime, const QString& unit) const
{
    if (!baseTime.isValid() || !currentTime.isValid()) {
        return 0.0;
    }

    // 计算秒差
    int diffSeconds = baseTime.secsTo(currentTime);

    // 处理跨日情况（当前时间小于基准时间时）
    if (diffSeconds < 0) {
        diffSeconds += 24 * 3600; // 加一天的秒数
    }

    return convertTimeToUnit(diffSeconds, unit);
}

double DataEditorWidget::convertTimeToUnit(double seconds, const QString& unit) const
{
    if (unit == "s") {
        return seconds;
    } else if (unit == "m") {
        return seconds / 60.0;
    } else if (unit == "h") {
        return seconds / 3600.0;
    }

    return seconds; // 默认返回秒
}

bool DataEditorWidget::isValidTimeFormat(const QString& timeStr) const
{
    return parseTimeString(timeStr).isValid();
}

// ============================================================================
// 压降计算功能实现 - 优化版本
// ============================================================================

void DataEditorWidget::onPressureDropCalc()
{
    if (!hasData()) {
        showStyledMessageBox("压降计算", "请先加载数据文件", QMessageBox::Information);
        return;
    }

    // 显示计算进度
    showAnimatedProgress("压降计算", "正在计算压力降...");

    // 执行压降计算
    PressureDropResult result = calculatePressureDrop();

    hideAnimatedProgress();

    if (result.success) {
        updateStatus(QString("压降计算完成 - 已添加列: %1").arg(result.columnName), "success");
        m_dataModified = true;
        emit pressureDropCalculated(result);
        emitDataChanged();

        showStyledMessageBox("压降计算完成",
                             QString("压降计算成功完成！\n"
                                     "新增列：%1\n"
                                     "处理行数：%2")
                                 .arg(result.columnName)
                                 .arg(result.processedRows),
                             QMessageBox::Information);
    } else {
        updateStatus("压降计算失败", "error");
        showStyledMessageBox("压降计算失败", result.errorMessage, QMessageBox::Warning);
    }
}

PressureDropResult DataEditorWidget::calculatePressureDrop()
{
    PressureDropResult result;
    result.success = false;
    result.addedColumnIndex = -1;
    result.processedRows = 0;

    if (!m_dataModel) {
        result.errorMessage = "数据模型不存在";
        return result;
    }

    // 查找压力列
    int pressureColumn = findPressureColumn();
    if (pressureColumn == -1) {
        result.errorMessage = "未找到压力列。";
        return result;
    }

    // 获取压力单位
    QString pressureUnit = getPressureUnit();

    // 创建压降列名
    QString dropColumnName = QString("压降\\%1").arg(pressureUnit.isEmpty() ? "MPa" : pressureUnit);

    // 在压力列后面插入新列
    int newColumnIndex = pressureColumn + 1;
    m_dataModel->insertColumn(newColumnIndex);

    // 设置列标题
    QStandardItem* headerItem = new QStandardItem(dropColumnName);
    m_dataModel->setHorizontalHeaderItem(newColumnIndex, headerItem);

    // 计算压降数据
    int rowCount = m_dataModel->rowCount();
    QList<double> pressureValues;

    // 收集所有压力数据
    for (int row = 0; row < rowCount; ++row) {
        QStandardItem* pressureItem = m_dataModel->item(row, pressureColumn);
        if (pressureItem) {
            QString pressureText = pressureItem->text().trimmed();
            if (isValidPressureData(pressureText)) {
                bool ok;
                double pressure = pressureText.toDouble(&ok);
                if (ok) {
                    pressureValues.append(pressure);
                } else {
                    pressureValues.append(0.0);
                }
            } else {
                pressureValues.append(0.0);
            }
        } else {
            pressureValues.append(0.0);
        }
    }

    // 计算压降值 - 修正的计算逻辑：每个时刻相对于初始时刻的压降
    double initialPressure = pressureValues.isEmpty() ? 0.0 : pressureValues[0]; // 获取初始压力

    for (int row = 0; row < rowCount; ++row) {
        double pressureDrop = 0.0;

        if (row == 0) {
            // 第1行（初始时刻）的压降为0
            pressureDrop = 0.0;
        } else {
            // 其他行的压降 = 初始时刻压力 - 当前时刻压力
            pressureDrop = initialPressure - pressureValues[row];
        }

        // 创建新的数据项
        QString dropValueText = QString::number(pressureDrop, 'f', 3);
        QStandardItem* dropItem = new QStandardItem(dropValueText);
        dropItem->setForeground(QBrush(QColor("#2c3e50")));
        m_dataModel->setItem(row, newColumnIndex, dropItem);

        result.processedRows++;
    }

    // 添加列定义
    ColumnDefinition newColumnDef;
    newColumnDef.name = dropColumnName;
    newColumnDef.type = WellTestColumnType::PressureDrop;
    newColumnDef.unit = pressureUnit.isEmpty() ? "MPa" : pressureUnit;
    newColumnDef.description = "压力降";
    newColumnDef.isRequired = false;
    newColumnDef.minValue = -999999;
    newColumnDef.maxValue = 999999;
    newColumnDef.decimalPlaces = 3;

    // 插入到列定义列表中
    if (newColumnIndex < m_columnDefinitions.size()) {
        m_columnDefinitions.insert(newColumnIndex, newColumnDef);
    } else {
        m_columnDefinitions.append(newColumnDef);
    }

    result.success = true;
    result.addedColumnIndex = newColumnIndex;
    result.columnName = dropColumnName;

    // 优化列宽
    optimizeColumnWidths();

    return result;
}

int DataEditorWidget::findPressureColumn() const
{
    if (!m_dataModel) {
        return -1;
    }

    // 优先查找已定义为压力的列
    for (int i = 0; i < m_columnDefinitions.size() && i < m_dataModel->columnCount(); ++i) {
        if (m_columnDefinitions[i].type == WellTestColumnType::Pressure) {
            return i;
        }
    }

    // 如果没有定义的压力列，尝试从列名推断
    for (int col = 0; col < m_dataModel->columnCount(); ++col) {
        QString headerText = m_dataModel->headerData(col, Qt::Horizontal).toString().toLower();
        if (headerText.contains("pressure") || headerText.contains("压力") ||
            headerText.contains("压强") || headerText == "p") {
            return col;
        }
    }

    return -1;
}

int DataEditorWidget::findTimeColumn() const
{
    if (!m_dataModel) {
        return -1;
    }

    // 优先查找已定义为时间的列
    for (int i = 0; i < m_columnDefinitions.size() && i < m_dataModel->columnCount(); ++i) {
        if (m_columnDefinitions[i].type == WellTestColumnType::Time) {
            return i;
        }
    }

    // 如果没有定义的时间列，尝试从列名推断
    for (int col = 0; col < m_dataModel->columnCount(); ++col) {
        QString headerText = m_dataModel->headerData(col, Qt::Horizontal).toString().toLower();
        if (headerText.contains("time") || headerText.contains("时间") || headerText == "t") {
            return col;
        }
    }

    return -1;
}

QString DataEditorWidget::getPressureUnit() const
{
    int pressureColumn = findPressureColumn();
    if (pressureColumn >= 0 && pressureColumn < m_columnDefinitions.size()) {
        return m_columnDefinitions[pressureColumn].unit;
    }
    return "MPa"; // 默认单位
}

bool DataEditorWidget::isValidPressureData(const QString& data) const
{
    if (data.isEmpty()) {
        return false;
    }

    bool ok;
    data.toDouble(&ok);
    return ok;
}

// ============================================================================
// 右键菜单槽函数实现
// ============================================================================

void DataEditorWidget::onTableContextMenuRequested(const QPoint& pos)
{
    if (!m_dataModel || !m_contextMenu) {
        return;
    }

    // 保存右键点击位置
    m_lastContextMenuPos = pos;

    // 获取点击位置的行列
    QModelIndex index = ui->dataTableView->indexAt(pos);

    // 启用/禁用相应的动作
    bool hasData = m_dataModel->rowCount() > 0 && m_dataModel->columnCount() > 0;
    bool hasSelection = ui->dataTableView->selectionModel()->hasSelection();

    m_addRowAboveAction->setEnabled(hasData || m_dataModel->columnCount() > 0);
    m_addRowBelowAction->setEnabled(hasData || m_dataModel->columnCount() > 0);
    m_deleteRowsAction->setEnabled(hasSelection && m_dataModel->rowCount() > 0);

    m_addColumnLeftAction->setEnabled(true);
    m_addColumnRightAction->setEnabled(true);
    m_deleteColumnsAction->setEnabled(hasSelection && m_dataModel->columnCount() > 0);

    // 显示菜单
    QPoint globalPos = ui->dataTableView->mapToGlobal(pos);
    m_contextMenu->exec(globalPos);
}

void DataEditorWidget::onAddRowAbove()
{
    if (!m_dataModel) return;

    int row = getRowFromPosition(m_lastContextMenuPos);
    if (row == -1) {
        row = 0; // 如果没有点击在有效行上，在顶部添加
    }

    RowEditCommand* command = new RowEditCommand(m_dataModel, RowEditCommand::Insert, row);
    m_undoStack->push(command);

    m_dataModified = true;
    updateStatus("已在上方添加一行", "success");
    ui->dataTableView->selectRow(row);
    updateDataInfo();
    emitDataChanged();
}

void DataEditorWidget::onAddRowBelow()
{
    if (!m_dataModel) return;

    int row = getRowFromPosition(m_lastContextMenuPos);
    if (row == -1) {
        row = m_dataModel->rowCount(); // 如果没有点击在有效行上，在底部添加
    } else {
        row++; // 在选中行的下方添加
    }

    RowEditCommand* command = new RowEditCommand(m_dataModel, RowEditCommand::Insert, row);
    m_undoStack->push(command);

    m_dataModified = true;
    updateStatus("已在下方添加一行", "success");
    ui->dataTableView->selectRow(row);
    updateDataInfo();
    emitDataChanged();
}

void DataEditorWidget::onDeleteSelectedRows()
{
    if (!m_dataModel || m_dataModel->rowCount() == 0) return;

    QList<int> selectedRows = getSelectedRows();
    if (selectedRows.isEmpty()) {
        return;
    }

    std::sort(selectedRows.begin(), selectedRows.end(), std::greater<int>());

    QMessageBox msgBox;
    msgBox.setWindowTitle("确认删除");
    QString deleteText = selectedRows.size() == 1 ?
                             QString("确定要删除第 %1 行吗？").arg(selectedRows.first() + 1) :
                             QString("确定要删除选中的 %1 行吗？").arg(selectedRows.size());
    msgBox.setText(deleteText);
    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msgBox.setDefaultButton(QMessageBox::No);

    if (msgBox.exec() == QMessageBox::Yes) {
        m_undoStack->beginMacro("删除多行");

        for (int row : selectedRows) {
            QStringList rowData;
            for (int col = 0; col < m_dataModel->columnCount(); ++col) {
                QStandardItem* item = m_dataModel->item(row, col);
                rowData.append(item ? item->text() : "");
            }

            RowEditCommand* command = new RowEditCommand(m_dataModel, RowEditCommand::Delete, row, rowData);
            m_undoStack->push(command);
        }

        m_undoStack->endMacro();

        m_dataModified = true;
        updateStatus(QString("已删除 %1 行").arg(selectedRows.size()), "success");
        updateDataInfo();
        emitDataChanged();
    }
}

void DataEditorWidget::onAddColumnLeft()
{
    if (!m_dataModel) return;

    int col = getColumnFromPosition(m_lastContextMenuPos);
    if (col == -1) {
        col = 0; // 如果没有点击在有效列上，在左侧添加
    }

    bool ok;
    QString headerText = QInputDialog::getText(this, "添加列",
                                               "请输入列标题:", QLineEdit::Normal,
                                               QString("列%1").arg(col + 1), &ok);
    if (!ok || headerText.isEmpty()) {
        return;
    }

    ColumnEditCommand* command = new ColumnEditCommand(m_dataModel, ColumnEditCommand::Insert, col, headerText);
    m_undoStack->push(command);

    m_dataModified = true;
    updateStatus("已在左侧添加一列", "success");
    ui->dataTableView->selectColumn(col);
    updateDataInfo();
    emitDataChanged();
}

void DataEditorWidget::onAddColumnRight()
{
    if (!m_dataModel) return;

    int col = getColumnFromPosition(m_lastContextMenuPos);
    if (col == -1) {
        col = m_dataModel->columnCount(); // 如果没有点击在有效列上，在右侧添加
    } else {
        col++; // 在选中列的右侧添加
    }

    bool ok;
    QString headerText = QInputDialog::getText(this, "添加列",
                                               "请输入列标题:", QLineEdit::Normal,
                                               QString("列%1").arg(col + 1), &ok);
    if (!ok || headerText.isEmpty()) {
        return;
    }

    ColumnEditCommand* command = new ColumnEditCommand(m_dataModel, ColumnEditCommand::Insert, col, headerText);
    m_undoStack->push(command);

    m_dataModified = true;
    updateStatus("已在右侧添加一列", "success");
    ui->dataTableView->selectColumn(col);
    updateDataInfo();
    emitDataChanged();
}

void DataEditorWidget::onDeleteSelectedColumns()
{
    if (!m_dataModel || m_dataModel->columnCount() == 0) return;

    QList<int> selectedColumns = getSelectedColumns();
    if (selectedColumns.isEmpty()) {
        return;
    }

    std::sort(selectedColumns.begin(), selectedColumns.end(), std::greater<int>());

    QString headerText = selectedColumns.size() == 1 ?
                             m_dataModel->headerData(selectedColumns.first(), Qt::Horizontal).toString() :
                             QString("%1个列").arg(selectedColumns.size());

    QMessageBox msgBox;
    msgBox.setWindowTitle("确认删除");
    msgBox.setText(QString("确定要删除列 \"%1\" 吗？").arg(headerText));
    msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
    msgBox.setDefaultButton(QMessageBox::No);

    if (msgBox.exec() == QMessageBox::Yes) {
        m_undoStack->beginMacro("删除多列");

        for (int col : selectedColumns) {
            QStandardItem* headerItem = m_dataModel->horizontalHeaderItem(col);
            QString headerName = headerItem ? headerItem->text() : QString("列%1").arg(col + 1);

            QStringList columnData;
            for (int row = 0; row < m_dataModel->rowCount(); ++row) {
                QStandardItem* item = m_dataModel->item(row, col);
                columnData.append(item ? item->text() : "");
            }

            ColumnEditCommand* command = new ColumnEditCommand(m_dataModel, ColumnEditCommand::Delete, col, headerName, columnData);
            m_undoStack->push(command);
        }

        m_undoStack->endMacro();

        m_dataModified = true;
        updateStatus(QString("已删除 %1 列").arg(selectedColumns.size()), "success");
        updateDataInfo();
        emitDataChanged();
    }
}

// ============================================================================
// 辅助方法实现
// ============================================================================

int DataEditorWidget::getRowFromPosition(const QPoint& pos) const
{
    QModelIndex index = ui->dataTableView->indexAt(pos);
    if (index.isValid()) {
        QModelIndex sourceIndex = m_proxyModel->mapToSource(index);
        return sourceIndex.row();
    }
    return -1;
}

int DataEditorWidget::getColumnFromPosition(const QPoint& pos) const
{
    QModelIndex index = ui->dataTableView->indexAt(pos);
    if (index.isValid()) {
        QModelIndex sourceIndex = m_proxyModel->mapToSource(index);
        return sourceIndex.column();
    }
    return -1;
}

// ============================================================================
// 文件操作槽函数实现 - 优化版本
// ============================================================================

void DataEditorWidget::onOpenFile()
{
    if (m_dataModified && !checkDataModifiedAndPrompt()) {
        return;
    }

    QString filter = "所有支持的文件 (*.xlsx *.xls *.csv *.txt *.json);;Excel Files (*.xlsx *.xls);;CSV Files (*.csv);;Text Files (*.txt);;JSON Files (*.json)";

    QString filePath = QFileDialog::getOpenFileName(this, "选择试井数据文件", QString(), filter);

    if (filePath.isEmpty()) {
        return;
    }

    QString fileType;
    QString extension = QFileInfo(filePath).suffix().toLower();

    if (extension == "xlsx" || extension == "xls") {
        fileType = "excel";
    } else if (extension == "csv") {
        fileType = "csv";
    } else if (extension == "txt") {
        fileType = "txt";
    } else if (extension == "json") {
        fileType = "json";
    } else {
        fileType = "txt";
    }

    // 对于CSV和TXT文件，显示读取配置对话框
    if (fileType == "csv" || fileType == "txt") {
        DataLoadConfigDialog configDialog(filePath, this);
        if (configDialog.exec() == QDialog::Accepted) {
            DataLoadConfigDialog::LoadConfig config = configDialog.getLoadConfig();
            loadDataWithConfig(filePath, fileType, config);
        }
    } else {
        // Excel和JSON文件直接加载
        loadData(filePath, fileType);
    }
}

void DataEditorWidget::loadData(const QString& filePath, const QString& fileType)
{
    qDebug() << "开始加载文件:" << filePath << "类型:" << fileType;

    QFileInfo fileInfo(filePath);
    if (!fileInfo.exists() || !fileInfo.isReadable()) {
        showStyledMessageBox("文件加载失败",
                             QString("文件不存在或无法读取: %1").arg(filePath),
                             QMessageBox::Warning);
        return;
    }

    // 显示进度对话框
    showAnimatedProgress("加载数据文件", "正在读取文件数据，请稍候...");

    clearData();

    m_currentFilePath = filePath;
    m_currentFileType = fileType;
    ui->filePathLineEdit->setText(filePath);

    bool loadSuccess = false;
    QString errorMessage;

    updateProgress(20, "正在分析文件格式...");

    QString lowerType = fileType.toLower();
    if (lowerType == "excel") {
        loadSuccess = loadExcelFileOptimized(filePath, errorMessage);
    } else if (lowerType == "txt" || lowerType == "csv") {
        loadSuccess = loadCsvFile(filePath, errorMessage);
    } else if (lowerType == "json") {
        loadSuccess = loadJsonFile(filePath, errorMessage);
    } else {
        errorMessage = QString("不支持的文件类型: %1").arg(fileType);
    }

    hideAnimatedProgress();

    if (loadSuccess) {
        updateStatus(QString("文件加载成功 - %1行 × %2列")
                         .arg(m_dataModel->rowCount())
                         .arg(m_dataModel->columnCount()), "success");

        setButtonsEnabled(true);
        m_dataModified = false;

        applyColumnStyles();
        optimizeColumnWidths();
        optimizeTableDisplay();

        // 弹出列定义对话框
        QTimer::singleShot(500, this, &DataEditorWidget::onDefineColumns);

        emitDataChanged();

        qDebug() << "文件加载成功，数据行数:" << m_dataModel->rowCount()
                 << "列数:" << m_dataModel->columnCount();
    } else {
        updateStatus("文件加载失败", "error");
        showStyledMessageBox("文件加载失败",
                             QString("无法加载文件: %1").arg(filePath),
                             QMessageBox::Critical,
                             errorMessage);
        qDebug() << "文件加载失败:" << errorMessage;
    }
}

// ============================================================================
// 优化的文件加载方法
// ============================================================================

bool DataEditorWidget::quickDetectFileFormat(const QString& filePath)
{
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        return false;
    }

    // 只读取前几行来快速检测格式
    QTextStream in(&file);
    QString firstLine = in.readLine();
    QString secondLine = in.readLine();
    file.close();

    // 检测是否是CSV格式
    if (firstLine.contains(',') || firstLine.contains('\t') || firstLine.contains(';')) {
        return true;
    }

    return false;
}

// bool DataEditorWidget::loadCsvFileWithConfig(const QString& filePath, const DataLoadConfigDialog::LoadConfig& config, QString& errorMessage)
// {
//     QFile file(filePath);
//     if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
//         errorMessage = QString("无法打开文件: %1").arg(file.errorString());
//         return false;
//     }

//     QStringList lines;
//     QString usedEncoding = config.encoding;

//     QTextStream in(&file);

// #if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
//     if (config.encoding != "UTF-8") {
//         QTextCodec* codec = QTextCodec::codecForName(config.encoding.toLocal8Bit());
//         if (codec) {
//             in.setCodec(codec);
//         }
//     }
// #else
//     if (config.encoding == "GBK" || config.encoding == "GB2312") {
//         in.setEncoding(QStringConverter::System);
//     }
// #endif

//     // 读取所有行
//     while (!in.atEnd() && lines.size() < m_maxDisplayRows) {
//         QString line = in.readLine();
//         lines.append(line);
//     }

//     file.close();

//     if (lines.isEmpty()) {
//         errorMessage = "文件为空或无法读取";
//         return false;
//     }

//     updateProgress(70, "正在解析数据格式...");

//     // 检查起始行是否有效
//     if (config.startRow > lines.size()) {
//         errorMessage = QString("起始行 %1 超出文件总行数 %2").arg(config.startRow).arg(lines.size());
//         return false;
//     }

//     if (lines.size() >= m_maxDisplayRows) {
//         m_largeFileMode = true;
//         qDebug() << "启用大文件模式，限制显示行数为" << m_maxDisplayRows;
//     }

//     // 确定表头
//     QStringList headers;
//     int dataStartIndex = config.startRow - 1; // 转换为0基索引

//     if (config.hasHeader && dataStartIndex < lines.size()) {
//         // 使用指定行作为表头
//         QStringList headerFields = splitCSVLine(lines[dataStartIndex], config.separator);
//         for (const QString& field : headerFields) {
//             QString header = field.trimmed();
//             if (header.isEmpty()) {
//                 header = QString("列%1").arg(headers.size() + 1);
//             }
//             headers.append(header);
//         }
//         dataStartIndex++; // 数据从表头的下一行开始
//     } else {
//         // 没有表头，使用第一行数据确定列数
//         if (dataStartIndex < lines.size()) {
//             QStringList firstDataFields = splitCSVLine(lines[dataStartIndex], config.separator);
//             for (int i = 0; i < firstDataFields.size(); ++i) {
//                 headers.append(QString("列%1").arg(i + 1));
//             }
//         }
//     }

//     if (headers.isEmpty()) {
//         errorMessage = "无法确定数据列结构";
//         return false;
//     }

//     m_dataModel->setColumnCount(headers.size());
//     m_dataModel->setHorizontalHeaderLabels(headers);

//     updateProgress(80, "正在加载数据...");

//     // 加载数据行
//     int rowIndex = 0;
//     for (int i = dataStartIndex; i < lines.size(); ++i) {
//         QStringList lineFields = splitCSVLine(lines[i], config.separator);

//         // 调整字段数量以匹配列数
//         while (lineFields.size() < headers.size()) {
//             lineFields.append("");
//         }
//         if (lineFields.size() > headers.size()) {
//             lineFields = lineFields.mid(0, headers.size());
//         }

//         // 插入新行
//         m_dataModel->insertRow(rowIndex);

//         // 填充数据
//         for (int col = 0; col < lineFields.size(); ++col) {
//             QString cellData = lineFields[col].trimmed();
//             QStandardItem* item = new QStandardItem(cellData);
//             item->setForeground(QBrush(QColor("#2c3e50")));
//             m_dataModel->setItem(rowIndex, col, item);
//         }
//         rowIndex++;

//         if (i % 100 == 0) {
//             updateProgress(80 + (i * 15 / lines.size()),
//                            QString("已加载 %1 行").arg(rowIndex));
//             QApplication::processEvents(); // 让界面保持响应
//         }
//     }

//     updateProgress(100, "数据加载完成");

//     qDebug() << "成功加载" << m_dataModel->rowCount() << "行数据，"
//              << m_dataModel->columnCount() << "列，使用编码:" << usedEncoding
//              << "，起始行:" << config.startRow;

//     return true;
// }

bool DataEditorWidget::loadCsvFileWithConfig(const QString& filePath, const DataLoadConfigDialog::LoadConfig& config, QString& errorMessage)
{
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        errorMessage = QString("无法打开文件: %1").arg(file.errorString());
        return false;
    }

    QStringList lines;
    QString usedEncoding = config.encoding;

    QTextStream in(&file);

#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    if (config.encoding != "UTF-8") {
        QTextCodec* codec = QTextCodec::codecForName(config.encoding.toLocal8Bit());
        if (codec) {
            in.setCodec(codec);
        }
    }
#else
    if (config.encoding == "GBK" || config.encoding == "GB2312") {
        in.setEncoding(QStringConverter::System);
    }
#endif

    // 读取所有行
    while (!in.atEnd() && lines.size() < m_maxDisplayRows) {
        QString line = in.readLine();
        lines.append(line);
    }

    file.close();

    if (lines.isEmpty()) {
        errorMessage = "文件为空或无法读取";
        return false;
    }

    updateProgress(70, "正在解析数据格式...");

    // 检查起始行是否有效
    if (config.startRow > lines.size()) {
        errorMessage = QString("起始行 %1 超出文件总行数 %2").arg(config.startRow).arg(lines.size());
        return false;
    }

    if (lines.size() >= m_maxDisplayRows) {
        m_largeFileMode = true;
        qDebug() << "启用大文件模式，限制显示行数为" << m_maxDisplayRows;
    }

    // 确定表头
    QStringList headers;
    int dataStartIndex = config.startRow - 1; // 转换为0基索引

    if (config.hasHeader && dataStartIndex < lines.size()) {
        // 使用指定行作为表头
        QStringList headerFields = splitCSVLine(lines[dataStartIndex], config.separator);
        for (const QString& field : headerFields) {
            QString header = field.trimmed();
            if (header.isEmpty()) {
                header = QString("列%1").arg(headers.size() + 1);
            }
            headers.append(header);
        }
        dataStartIndex++; // 数据从表头的下一行开始
    } else {
        // 没有表头，使用第一行数据确定列数
        if (dataStartIndex < lines.size()) {
            QStringList firstDataFields = splitCSVLine(lines[dataStartIndex], config.separator);
            for (int i = 0; i < firstDataFields.size(); ++i) {
                headers.append(QString("列%1").arg(i + 1));
            }
        }
    }

    if (headers.isEmpty()) {
        errorMessage = "无法确定数据列结构";
        return false;
    }

    m_dataModel->setColumnCount(headers.size());
    m_dataModel->setHorizontalHeaderLabels(headers);

    updateProgress(80, "正在加载数据...");

    // 加载数据行
    int rowIndex = 0;
    for (int i = dataStartIndex; i < lines.size(); ++i) {
        QStringList lineFields = splitCSVLine(lines[i], config.separator);

        // 调整字段数量以匹配列数
        while (lineFields.size() < headers.size()) {
            lineFields.append("");
        }
        if (lineFields.size() > headers.size()) {
            lineFields = lineFields.mid(0, headers.size());
        }

        // 插入新行
        m_dataModel->insertRow(rowIndex);

        // 填充数据
        for (int col = 0; col < lineFields.size(); ++col) {
            QString cellData = lineFields[col].trimmed();
            QStandardItem* item = new QStandardItem(cellData);
            item->setForeground(QBrush(QColor("#2c3e50")));
            m_dataModel->setItem(rowIndex, col, item);
        }
        rowIndex++;

        if (i % 100 == 0) {
            updateProgress(80 + (i * 15 / lines.size()),
                           QString("已加载 %1 行").arg(rowIndex));
            QApplication::processEvents(); // 让界面保持响应
        }
    }

    updateProgress(100, "数据加载完成");

    qDebug() << "成功加载" << m_dataModel->rowCount() << "行数据，"
             << m_dataModel->columnCount() << "列，使用编码:" << usedEncoding
             << "，起始行:" << config.startRow;

    return true;
}

QString DataEditorWidget::detectOptimalSeparator(const QString& filePath)
{
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        return ",";
    }

    QTextStream in(&file);
    QString sampleText;
    int lineCount = 0;
    while (!in.atEnd() && lineCount < 10) {
        sampleText += in.readLine() + "\n";
        lineCount++;
    }
    file.close();

    // 统计各种分隔符的出现频率
    QMap<QString, int> separatorCounts;
    separatorCounts[","] = sampleText.count(',');
    separatorCounts["\t"] = sampleText.count('\t');
    separatorCounts[";"] = sampleText.count(';');
    separatorCounts["|"] = sampleText.count('|');

    // 找到出现频率最高的分隔符
    QString bestSeparator = ",";
    int maxCount = 0;
    for (auto it = separatorCounts.begin(); it != separatorCounts.end(); ++it) {
        if (it.value() > maxCount) {
            maxCount = it.value();
            bestSeparator = it.key();
        }
    }

    return bestSeparator;
}

bool DataEditorWidget::loadExcelFileOptimized(const QString& filePath, QString& errorMessage)
{
    qDebug() << "尝试优化加载Excel文件:" << filePath;

    updateProgress(30, "检测文件格式...");

    // 首先快速检测是否为CSV格式的Excel文件
    if (quickDetectFileFormat(filePath)) {
        updateProgress(50, "检测到CSV格式，使用快速读取...");

        // 使用最优分隔符直接加载
        QString separator = detectOptimalSeparator(filePath);
        qDebug() << "检测到最优分隔符:" << separator;

        if (loadCSVFile(filePath, separator, errorMessage)) {
            qDebug() << "使用分隔符'" << separator << "'快速读取Excel成功";
            return true;
        }
    }

    updateProgress(60, "尝试COM组件读取...");

#ifdef Q_OS_WIN
    if (loadExcelWithCOM(filePath, errorMessage)) {
        return true;
    }
    qDebug() << "COM方式失败，尝试通用方式:" << errorMessage;
#endif

    updateProgress(80, "使用通用方式读取...");
    return loadExcelAsCSV(filePath, errorMessage);
}

bool DataEditorWidget::loadCsvFile(const QString& filePath, QString& errorMessage)
{
    updateProgress(30, "检测最佳分隔符...");

    // 使用优化的分隔符检测
    QString bestSeparator = detectOptimalSeparator(filePath);

    updateProgress(50, QString("使用分隔符 '%1' 读取数据...").arg(bestSeparator));

    if (loadCSVFile(filePath, bestSeparator, errorMessage)) {
        qDebug() << "使用最优分隔符'" << bestSeparator << "'成功读取CSV文件";
        return true;
    }

    // 如果最优分隔符失败，尝试其他分隔符
    QStringList otherSeparators = {",", "\t", ";", "|"};
    otherSeparators.removeOne(bestSeparator);

    for (const QString& separator : otherSeparators) {
        updateProgress(60 + otherSeparators.indexOf(separator) * 10,
                       QString("尝试分隔符 '%1'...").arg(separator));

        if (loadCSVFile(filePath, separator, errorMessage)) {
            qDebug() << "使用分隔符'" << separator << "'成功读取CSV文件";
            return true;
        }
    }

    errorMessage = "无法确定CSV文件的分隔符格式";
    return false;
}

bool DataEditorWidget::loadCSVFile(const QString& filePath, const QString& separator, QString& errorMessage)
{
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        errorMessage = QString("无法打开文件: %1").arg(file.errorString());
        return false;
    }

    QStringList lines;
    QString usedEncoding = "UTF-8";

#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    QStringList encodings = {"UTF-8", "GBK"};

    for (const QString& encoding : encodings) {
        file.seek(0);
        QTextStream in(&file);

        if (encoding != "UTF-8") {
            QTextCodec* codec = QTextCodec::codecForName(encoding.toLocal8Bit());
            if (codec) {
                in.setCodec(codec);
            }
        }

        lines.clear();
        while (!in.atEnd() && lines.size() < m_maxDisplayRows) {
            QString line = in.readLine().trimmed();
            if (!line.isEmpty()) {
                lines.append(line);
            }
        }

        if (!lines.isEmpty()) {
            usedEncoding = encoding;
            break;
        }
    }
#else
    QTextStream in(&file);
    while (!in.atEnd() && lines.size() < m_maxDisplayRows) {
        QString line = in.readLine().trimmed();
        if (!line.isEmpty()) {
            lines.append(line);
        }
    }

    if (lines.isEmpty()) {
        file.seek(0);
        QTextStream in2(&file);
        in2.setEncoding(QStringConverter::Latin1);

        while (!in2.atEnd() && lines.size() < m_maxDisplayRows) {
            QString line = in2.readLine().trimmed();
            if (!line.isEmpty()) {
                lines.append(line);
            }
        }
        usedEncoding = "Latin-1";
    }
#endif

    file.close();

    if (lines.isEmpty()) {
        errorMessage = "文件为空或无法读取";
        return false;
    }

    updateProgress(70, "正在解析数据格式...");

    if (lines.size() >= m_maxDisplayRows) {
        m_largeFileMode = true;
        qDebug() << "启用大文件模式，限制显示行数为" << m_maxDisplayRows;
    }

    QString firstLine = lines.first();
    QStringList fields = splitCSVLine(firstLine, separator);

    if (fields.size() < 2) {
        return false;
    }

    // 快速验证数据一致性（只检查前5行）
    int expectedFields = fields.size();
    int validLines = 0;

    for (int i = 0; i < qMin(5, lines.size()); ++i) {
        QStringList lineFields = splitCSVLine(lines[i], separator);
        if (lineFields.size() == expectedFields) {
            validLines++;
        }
    }

    if (validLines < qMin(5, lines.size()) * 0.6) {
        return false;
    }

    updateProgress(80, "正在加载数据...");

    QStringList headers = fields;
    for (int i = 0; i < headers.size(); ++i) {
        headers[i] = headers[i].trimmed();
        if (headers[i].isEmpty()) {
            headers[i] = QString("列%1").arg(i + 1);
        }
    }

    m_dataModel->setColumnCount(headers.size());
    m_dataModel->setHorizontalHeaderLabels(headers);

    // 检查第一行是否为表头
    bool firstRowIsHeader = false;
    for (const QString& field : fields) {
        bool isNumber;
        field.toDouble(&isNumber);
        if (!isNumber && !field.isEmpty()) {
            firstRowIsHeader = true;
            break;
        }
    }

    int dataStartRow = firstRowIsHeader ? 1 : 0;
    int rowIndex = 0;

    // 批量插入行以提高性能
    int totalDataRows = lines.size() - dataStartRow;
    m_dataModel->setRowCount(totalDataRows);

    for (int i = dataStartRow; i < lines.size(); ++i) {
        QStringList lineFields = splitCSVLine(lines[i], separator);

        while (lineFields.size() < headers.size()) {
            lineFields.append("");
        }

        if (lineFields.size() > headers.size()) {
            lineFields = lineFields.mid(0, headers.size());
        }

        for (int col = 0; col < lineFields.size(); ++col) {
            QString cellData = lineFields[col].trimmed();
            QStandardItem* item = new QStandardItem(cellData);
            item->setForeground(QBrush(QColor("#2c3e50")));
            m_dataModel->setItem(rowIndex, col, item);
        }
        rowIndex++;

        if (i % 500 == 0) {
            updateProgress(80 + (i * 15 / lines.size()),
                           QString("已加载 %1/%2 行").arg(i).arg(lines.size()));
            QApplication::processEvents(); // 让界面保持响应
        }
    }

    updateProgress(100, "数据加载完成");

    qDebug() << "成功加载" << m_dataModel->rowCount() << "行数据，"
             << m_dataModel->columnCount() << "列，使用编码:" << usedEncoding;

    return true;
}

QStringList DataEditorWidget::splitCSVLine(const QString& line, const QString& separator)
{
    QStringList result;
    QString current;
    bool inQuotes = false;

    for (int i = 0; i < line.length(); ++i) {
        QChar ch = line.at(i);

        if (ch == '"') {
            inQuotes = !inQuotes;
        } else if (!inQuotes && line.mid(i, separator.length()) == separator) {
            result.append(current.trimmed());
            current.clear();
            i += separator.length() - 1;
        } else {
            current.append(ch);
        }
    }

    result.append(current.trimmed());
    return result;
}

bool DataEditorWidget::loadExcelFile(const QString& filePath, QString& errorMessage)
{
    return loadExcelFileOptimized(filePath, errorMessage);
}

bool DataEditorWidget::loadExcelAsCSV(const QString& filePath, QString& errorMessage)
{
    QStringList separators = {",", "\t", ";", "|"};

    for (const QString& separator : separators) {
        if (loadCSVFile(filePath, separator, errorMessage)) {
            qDebug() << "使用分隔符'" << separator << "'成功读取文件";
            return true;
        }
    }

    errorMessage = "无法以任何CSV格式读取此Excel文件。请尝试将Excel文件另存为CSV格式后重新加载。";
    return false;
}

bool DataEditorWidget::loadJsonFile(const QString& filePath, QString& errorMessage)
{
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly)) {
        errorMessage = QString("无法打开文件: %1").arg(file.errorString());
        return false;
    }

    QByteArray data = file.readAll();
    file.close();

    QJsonParseError parseError;
    QJsonDocument doc = QJsonDocument::fromJson(data, &parseError);

    if (parseError.error != QJsonParseError::NoError) {
        errorMessage = QString("JSON解析错误: %1").arg(parseError.errorString());
        return false;
    }

    if (doc.isArray()) {
        QJsonArray array = doc.array();
        if (array.isEmpty()) {
            errorMessage = "JSON文件中没有数据";
            return false;
        }

        QJsonObject firstObj = array.first().toObject();
        QStringList headers = firstObj.keys();

        m_dataModel->setColumnCount(headers.size());
        m_dataModel->setHorizontalHeaderLabels(headers);

        for (int i = 0; i < array.size(); ++i) {
            QJsonObject obj = array[i].toObject();
            m_dataModel->insertRow(i);

            for (int col = 0; col < headers.size(); ++col) {
                QString key = headers[col];
                QString value = obj[key].toString();

                QStandardItem* item = new QStandardItem(value);
                item->setForeground(QBrush(QColor("#2c3e50")));
                m_dataModel->setItem(i, col, item);
            }
        }

        return true;
    }

    errorMessage = "不支持的JSON格式";
    return false;
}

#ifdef Q_OS_WIN
bool DataEditorWidget::loadExcelWithCOM(const QString& filePath, QString& errorMessage)
{
    try {
        QAxObject excel("Excel.Application");
        if (excel.isNull()) {
            errorMessage = "无法创建Excel.Application对象，请确保已安装Microsoft Excel";
            return false;
        }

        excel.setProperty("Visible", false);
        excel.setProperty("DisplayAlerts", false);

        QAxObject* workbooks = excel.querySubObject("Workbooks");
        QAxObject* workbook = workbooks->querySubObject("Open(const QString&)", filePath);

        if (!workbook) {
            errorMessage = "无法打开Excel文件";
            return false;
        }

        QAxObject* worksheets = workbook->querySubObject("Worksheets");
        QAxObject* worksheet = worksheets->querySubObject("Item(int)", 1);
        QAxObject* usedRange = worksheet->querySubObject("UsedRange");

        if (!usedRange) {
            errorMessage = "Excel文件中没有数据";
            workbook->dynamicCall("Close()");
            excel.dynamicCall("Quit()");
            return false;
        }

        QAxObject* rows = usedRange->querySubObject("Rows");
        QAxObject* columns = usedRange->querySubObject("Columns");

        int rowCount = rows->property("Count").toInt();
        int columnCount = columns->property("Count").toInt();

        if (rowCount == 0 || columnCount == 0) {
            errorMessage = "Excel文件中没有有效数据";
            workbook->dynamicCall("Close()");
            excel.dynamicCall("Quit()");
            return false;
        }

        m_dataModel->setRowCount(rowCount > 1 ? rowCount-1 : 0);
        m_dataModel->setColumnCount(columnCount);

        // 设置表头
        QStringList headers;
        for (int col = 1; col <= columnCount; ++col) {
            QAxObject* cell = worksheet->querySubObject("Cells(int,int)", 1, col);
            QString headerText = cell ? cell->property("Value").toString() : QString("列%1").arg(col);
            if (headerText.isEmpty()) {
                headerText = QString("列%1").arg(col);
            }
            headers.append(headerText);
        }
        m_dataModel->setHorizontalHeaderLabels(headers);

        // 读取数据
        for (int row = 2; row <= rowCount; ++row) {
            for (int col = 1; col <= columnCount; ++col) {
                QAxObject* cell = worksheet->querySubObject("Cells(int,int)", row, col);
                QString value = cell ? cell->property("Value").toString() : "";
                QStandardItem* item = new QStandardItem(value);
                item->setForeground(QBrush(QColor("#2c3e50")));
                m_dataModel->setItem(row-2, col-1, item);
            }
        }

        workbook->dynamicCall("Close()");
        excel.dynamicCall("Quit()");

        return true;

    } catch (...) {
        errorMessage = "读取Excel文件时发生未知错误";
        return false;
    }
}
#endif

void DataEditorWidget::onSave()
{
    if (m_currentFilePath.isEmpty() || m_currentFileType.isEmpty()) {
        showStyledMessageBox("保存失败", "没有加载文件，无法保存", QMessageBox::Warning);
        return;
    }

    showAnimatedProgress("保存文件", "正在保存数据...");

    bool saveSuccess = false;
    QString lowerType = m_currentFileType.toLower();

    if (lowerType == "excel") {
        saveSuccess = saveExcelFile(m_currentFilePath);
    } else if (lowerType == "txt" || lowerType == "csv") {
        saveSuccess = saveCsvFile(m_currentFilePath);
    } else if (lowerType == "json") {
        saveSuccess = saveJsonFile(m_currentFilePath);
    }

    hideAnimatedProgress();

    if (saveSuccess) {
        updateStatus("文件保存成功", "success");
        m_dataModified = false;
        showStyledMessageBox("保存成功", "文件已成功保存。", QMessageBox::Information);
        emitDataChanged();
    } else {
        updateStatus("文件保存失败", "error");
        showStyledMessageBox("保存失败", "保存文件时出错。", QMessageBox::Critical);
    }
}

void DataEditorWidget::onExport()
{
    if (!hasData()) {
        showStyledMessageBox("导出失败", "没有数据可供导出。", QMessageBox::Warning);
        return;
    }

    QString filter = "CSV Files (*.csv);;Excel Files (*.xlsx);;JSON Files (*.json);;PDF Files (*.pdf);;HTML Files (*.html)";

    QString saveFilePath = QFileDialog::getSaveFileName(this, "导出试井数据", QString(), filter);

    if (saveFilePath.isEmpty()) {
        return;
    }

    showAnimatedProgress("导出文件", "正在导出数据...");

    bool saveSuccess = false;
    QString extension = QFileInfo(saveFilePath).suffix().toLower();

    if (extension == "xlsx") {
        saveSuccess = saveExcelFile(saveFilePath);
    } else if (extension == "csv") {
        saveSuccess = saveCsvFile(saveFilePath);
    } else if (extension == "json") {
        saveSuccess = saveJsonFile(saveFilePath);
    } else if (extension == "pdf") {
        saveSuccess = exportToPdf(saveFilePath);
    } else if (extension == "html") {
        saveSuccess = exportToHtml(saveFilePath);
    } else {
        saveSuccess = saveCsvFile(saveFilePath);
    }

    hideAnimatedProgress();

    if (saveSuccess) {
        showStyledMessageBox("导出成功", QString("文件已成功导出到: %1").arg(saveFilePath), QMessageBox::Information);
    } else {
        showStyledMessageBox("导出失败", "导出文件时出错。", QMessageBox::Critical);
    }
}

// ============================================================================
// 数据处理槽函数实现
// ============================================================================

void DataEditorWidget::onDefineColumns()
{
    if (!hasData()) {
        showStyledMessageBox("列定义", "请先加载数据文件", QMessageBox::Information);
        return;
    }

    QStringList columnNames;
    for (int i = 0; i < m_dataModel->columnCount(); ++i) {
        columnNames.append(m_dataModel->headerData(i, Qt::Horizontal).toString());
    }

    ColumnDefinitionDialog dialog(columnNames, m_columnDefinitions, this);
    if (dialog.exec() == QDialog::Accepted) {
        m_columnDefinitions = dialog.getColumnDefinitions();

        // 更新列标题
        updateColumnHeaders();

        updateStatus("列定义已更新", "success");
        emit columnDefinitionsChanged();
        emitDataChanged();
    }
}

void DataEditorWidget::updateColumnHeaders()
{
    if (!m_dataModel) return;

    for (int i = 0; i < m_columnDefinitions.size() && i < m_dataModel->columnCount(); ++i) {
        const ColumnDefinition& def = m_columnDefinitions[i];

        // 直接使用组合的列名（类型\单位）
        QString newHeaderText = def.name;

        // 更新模型中的列标题
        m_dataModel->setHeaderData(i, Qt::Horizontal, newHeaderText);

        // 应用其他列定义属性
        applyColumnDefinition(i, def);
    }

    // 刷新表格显示
    ui->dataTableView->update();
    optimizeColumnWidths();
}

void DataEditorWidget::onDataClean()
{
    if (!hasData()) {
        showStyledMessageBox("数据清理", "请先加载数据文件", QMessageBox::Information);
        return;
    }

    DataCleaningDialog dialog(this);
    if (dialog.exec() == QDialog::Accepted) {
        DataCleaningDialog::CleaningOptions options = dialog.getCleaningOptions();

        showAnimatedProgress("数据清理", "正在清理数据...");

        int cleanedCount = 0;

        if (options.removeEmptyRows) {
            removeEmptyRows();
            cleanedCount++;
            updateProgress(20, "删除空行...");
        }

        if (options.removeEmptyColumns) {
            removeEmptyColumns();
            cleanedCount++;
            updateProgress(40, "删除空列...");
        }

        if (options.removeDuplicates) {
            removeDuplicates();
            cleanedCount++;
            updateProgress(60, "删除重复行...");
        }

        if (options.fillMissingValues) {
            fillMissingValues(options.fillMethod);
            cleanedCount++;
            updateProgress(80, "填充缺失值...");
        }

        if (options.removeOutliers) {
            removeOutliers(options.outlierThreshold);
            cleanedCount++;
            updateProgress(90, "删除异常值...");
        }

        if (options.standardizeFormat) {
            standardizeDataFormat();
            cleanedCount++;
            updateProgress(100, "标准化格式...");
        }

        hideAnimatedProgress();

        if (cleanedCount > 0) {
            updateStatus("数据清理完成", "success");
            m_dataModified = true;
            emitDataChanged();
            showStyledMessageBox("数据清理",
                                 QString("数据清理完成，执行了 %1 项清理操作").arg(cleanedCount),
                                 QMessageBox::Information);
        } else {
            showStyledMessageBox("数据清理", "未选择任何清理操作", QMessageBox::Information);
        }
    }
}

void DataEditorWidget::onDataStatistics()
{
    if (!hasData()) {
        showStyledMessageBox("统计分析", "没有数据可供分析", QMessageBox::Information);
        return;
    }

    showAnimatedProgress("数据统计", "正在计算统计信息...");

    QList<DataStatistics> statistics = calculateAllStatistics();

    hideAnimatedProgress();

    QString statisticsText = "试井数据统计分析结果:\n\n";

    for (const DataStatistics& stat : statistics) {
        statisticsText += QString("列名: %1\n").arg(stat.columnName);
        statisticsText += QString("数据类型: %1\n").arg(stat.dataType);
        statisticsText += QString("总计数据: %1\n").arg(stat.dataCount);
        statisticsText += QString("有效数据: %1\n").arg(stat.validCount);
        statisticsText += QString("无效数据: %1\n").arg(stat.invalidCount);

        if (stat.dataType == "数值型" && stat.validCount > 0) {
            statisticsText += QString("最小值: %1 %2\n").arg(formatNumber(stat.minimum)).arg(stat.unit);
            statisticsText += QString("最大值: %1 %2\n").arg(formatNumber(stat.maximum)).arg(stat.unit);
            statisticsText += QString("平均值: %1 %2\n").arg(formatNumber(stat.average)).arg(stat.unit);
            statisticsText += QString("中位数: %1 %2\n").arg(formatNumber(stat.median)).arg(stat.unit);
            statisticsText += QString("标准差: %1 %2\n").arg(formatNumber(stat.standardDeviation)).arg(stat.unit);
        }

        statisticsText += "\n" + QString("-").repeated(50) + "\n\n";
    }

    showStyledMessageBox("试井数据统计分析", "统计分析完成", QMessageBox::Information, statisticsText);

    emit statisticsCalculated(statistics);
}

// ============================================================================
// 搜索功能实现
// ============================================================================

void DataEditorWidget::onSearchTextChanged()
{
    m_currentSearchText = ui->searchLineEdit->text();
    m_searchTimer->stop();
    m_searchTimer->start();
}

void DataEditorWidget::onSearchData()
{
    QString searchText = m_currentSearchText.trimmed();

    if (searchText.isEmpty()) {
        clearDataFilter();
        updateStatus("就绪", "success");
    } else {
        applyDataFilter(searchText);
        int matchCount = m_proxyModel->rowCount();
        updateStatus(QString("找到 %1 条匹配记录").arg(matchCount), "info");
        emit searchCompleted(matchCount);
    }
}

void DataEditorWidget::applyDataFilter(const QString& filterText)
{
    if (m_proxyModel) {
        m_proxyModel->setFilterWildcard(filterText);
    }
}

void DataEditorWidget::clearDataFilter()
{
    if (m_proxyModel) {
        m_proxyModel->setFilterWildcard("");
    }
}

// ============================================================================
// 核心数据处理功能实现（简化版本）
// ============================================================================

DataStatistics DataEditorWidget::calculateColumnStatistics(int column) const
{
    DataStatistics stats;

    if (!m_dataModel || column < 0 || column >= m_dataModel->columnCount()) {
        return stats;
    }

    stats.columnName = m_dataModel->headerData(column, Qt::Horizontal).toString();

    // 获取列定义中的单位
    if (column < m_columnDefinitions.size()) {
        stats.unit = m_columnDefinitions[column].unit;
    }

    stats.dataCount = m_dataModel->rowCount();
    stats.validCount = 0;
    stats.invalidCount = 0;

    QList<double> numericValues;
    QStringList textValues;

    for (int row = 0; row < m_dataModel->rowCount(); ++row) {
        QStandardItem* item = m_dataModel->item(row, column);
        QString value = item ? item->text().trimmed() : "";

        if (value.isEmpty()) {
            stats.invalidCount++;
            continue;
        }

        bool isNumeric;
        double numValue = value.toDouble(&isNumeric);

        if (isNumeric) {
            numericValues.append(numValue);
            stats.validCount++;
        } else {
            textValues.append(value);
            stats.validCount++;
        }
    }

    if (numericValues.size() > textValues.size()) {
        stats.dataType = "数值型";

        if (!numericValues.isEmpty()) {
            std::sort(numericValues.begin(), numericValues.end());

            stats.minimum = numericValues.first();
            stats.maximum = numericValues.last();

            double sum = 0;
            for (double val : numericValues) {
                sum += val;
            }
            stats.average = sum / numericValues.size();

            int size = numericValues.size();
            if (size % 2 == 0) {
                stats.median = (numericValues[size/2 - 1] + numericValues[size/2]) / 2.0;
            } else {
                stats.median = numericValues[size/2];
            }

            double variance = 0;
            for (double val : numericValues) {
                variance += std::pow(val - stats.average, 2);
            }
            stats.standardDeviation = std::sqrt(variance / numericValues.size());
        }
    } else {
        stats.dataType = "文本型";
        stats.minimum = 0;
        stats.maximum = 0;
        stats.average = 0;
        stats.median = 0;
        stats.standardDeviation = 0;
    }

    return stats;
}

QList<DataStatistics> DataEditorWidget::calculateAllStatistics() const
{
    QList<DataStatistics> allStats;

    if (!m_dataModel) {
        return allStats;
    }

    for (int col = 0; col < m_dataModel->columnCount(); ++col) {
        allStats.append(calculateColumnStatistics(col));
    }

    return allStats;
}

// ============================================================================
// 简化的数据清理功能
// ============================================================================

void DataEditorWidget::removeEmptyRows()
{
    if (!m_dataModel) return;

    QList<int> emptyRows;

    for (int row = 0; row < m_dataModel->rowCount(); ++row) {
        bool isEmpty = true;
        for (int col = 0; col < m_dataModel->columnCount(); ++col) {
            QStandardItem* item = m_dataModel->item(row, col);
            if (item && !item->text().trimmed().isEmpty()) {
                isEmpty = false;
                break;
            }
        }
        if (isEmpty) {
            emptyRows.append(row);
        }
    }

    if (!emptyRows.isEmpty()) {
        std::sort(emptyRows.begin(), emptyRows.end(), std::greater<int>());

        m_undoStack->beginMacro("删除空行");
        for (int row : emptyRows) {
            QStringList rowData;
            for (int col = 0; col < m_dataModel->columnCount(); ++col) {
                QStandardItem* item = m_dataModel->item(row, col);
                rowData.append(item ? item->text() : "");
            }

            RowEditCommand* command = new RowEditCommand(m_dataModel, RowEditCommand::Delete, row, rowData);
            m_undoStack->push(command);
        }
        m_undoStack->endMacro();
    }
}

void DataEditorWidget::removeEmptyColumns()
{
    if (!m_dataModel) return;

    QList<int> emptyColumns;

    for (int col = 0; col < m_dataModel->columnCount(); ++col) {
        bool isEmpty = true;
        for (int row = 0; row < m_dataModel->rowCount(); ++row) {
            QStandardItem* item = m_dataModel->item(row, col);
            if (item && !item->text().trimmed().isEmpty()) {
                isEmpty = false;
                break;
            }
        }
        if (isEmpty) {
            emptyColumns.append(col);
        }
    }

    if (!emptyColumns.isEmpty()) {
        std::sort(emptyColumns.begin(), emptyColumns.end(), std::greater<int>());

        m_undoStack->beginMacro("删除空列");
        for (int col : emptyColumns) {
            QStandardItem* headerItem = m_dataModel->horizontalHeaderItem(col);
            QString headerName = headerItem ? headerItem->text() : QString("列%1").arg(col + 1);

            QStringList columnData;
            for (int row = 0; row < m_dataModel->rowCount(); ++row) {
                QStandardItem* item = m_dataModel->item(row, col);
                columnData.append(item ? item->text() : "");
            }

            ColumnEditCommand* command = new ColumnEditCommand(m_dataModel, ColumnEditCommand::Delete, col, headerName, columnData);
            m_undoStack->push(command);
        }
        m_undoStack->endMacro();
    }
}

void DataEditorWidget::removeDuplicates()
{
    if (!m_dataModel) return;

    QSet<QString> uniqueRows;
    QList<int> duplicateRows;

    for (int row = 0; row < m_dataModel->rowCount(); ++row) {
        QStringList rowData;
        for (int col = 0; col < m_dataModel->columnCount(); ++col) {
            QStandardItem* item = m_dataModel->item(row, col);
            rowData.append(item ? item->text().trimmed() : "");
        }

        QString rowSignature = rowData.join("|");
        if (uniqueRows.contains(rowSignature)) {
            duplicateRows.append(row);
        } else {
            uniqueRows.insert(rowSignature);
        }
    }

    if (!duplicateRows.isEmpty()) {
        std::sort(duplicateRows.begin(), duplicateRows.end(), std::greater<int>());

        m_undoStack->beginMacro("删除重复行");
        for (int row : duplicateRows) {
            QStringList rowData;
            for (int col = 0; col < m_dataModel->columnCount(); ++col) {
                QStandardItem* item = m_dataModel->item(row, col);
                rowData.append(item ? item->text() : "");
            }

            RowEditCommand* command = new RowEditCommand(m_dataModel, RowEditCommand::Delete, row, rowData);
            m_undoStack->push(command);
        }
        m_undoStack->endMacro();
    }
}

void DataEditorWidget::fillMissingValues(const QString& method)
{
    if (!m_dataModel) return;

    for (int col = 0; col < m_dataModel->columnCount(); ++col) {
        QList<double> numericValues;
        QList<int> validIndices;

        // 收集有效的数值
        for (int row = 0; row < m_dataModel->rowCount(); ++row) {
            QStandardItem* item = m_dataModel->item(row, col);
            if (item) {
                bool ok;
                double value = item->text().toDouble(&ok);
                if (ok) {
                    numericValues.append(value);
                    validIndices.append(row);
                }
            }
        }

        if (numericValues.isEmpty()) continue;

        // 填充缺失值
        for (int row = 0; row < m_dataModel->rowCount(); ++row) {
            QStandardItem* item = m_dataModel->item(row, col);
            if (!item || item->text().trimmed().isEmpty()) {
                QString fillValue;

                if (method == "zero") {
                    fillValue = "0";
                } else if (method == "average") {
                    double sum = 0;
                    for (double val : numericValues) {
                        sum += val;
                    }
                    fillValue = QString::number(sum / numericValues.size(), 'f', 3);
                } else if (method == "interpolation") {
                    // 简单的线性插值
                    if (!validIndices.isEmpty()) {
                        fillValue = QString::number(numericValues.first(), 'f', 3);
                    }
                } else if (method == "forward") {
                    // 前值填充
                    for (int prevRow = row - 1; prevRow >= 0; --prevRow) {
                        QStandardItem* prevItem = m_dataModel->item(prevRow, col);
                        if (prevItem && !prevItem->text().trimmed().isEmpty()) {
                            fillValue = prevItem->text();
                            break;
                        }
                    }
                }

                if (!fillValue.isEmpty()) {
                    if (!item) {
                        item = new QStandardItem(fillValue);
                        m_dataModel->setItem(row, col, item);
                    } else {
                        item->setText(fillValue);
                    }
                    item->setForeground(QBrush(QColor("#6c757d"))); // 标记为填充值
                }
            }
        }
    }
}

void DataEditorWidget::removeOutliers(double threshold)
{
    if (!m_dataModel) return;

    for (int col = 0; col < m_dataModel->columnCount(); ++col) {
        QList<double> values;
        QList<int> validRows;

        // 收集数值数据
        for (int row = 0; row < m_dataModel->rowCount(); ++row) {
            QStandardItem* item = m_dataModel->item(row, col);
            if (item) {
                bool ok;
                double value = item->text().toDouble(&ok);
                if (ok) {
                    values.append(value);
                    validRows.append(row);
                }
            }
        }

        if (values.size() < 3) continue; // 数据太少，跳过

        // 计算平均值和标准差
        double sum = 0;
        for (double val : values) {
            sum += val;
        }
        double mean = sum / values.size();

        double variance = 0;
        for (double val : values) {
            variance += std::pow(val - mean, 2);
        }
        double stdDev = std::sqrt(variance / values.size());

        // 标记异常值
        QList<int> outlierRows;
        for (int i = 0; i < values.size(); ++i) {
            if (std::abs(values[i] - mean) > threshold * stdDev) {
                outlierRows.append(validRows[i]);
            }
        }

        // 删除异常值（从后往前删除）
        if (!outlierRows.isEmpty()) {
            std::sort(outlierRows.begin(), outlierRows.end(), std::greater<int>());

            for (int row : outlierRows) {
                QStandardItem* item = m_dataModel->item(row, col);
                if (item) {
                    item->setText(""); // 清空异常值
                }
            }
        }
    }
}

void DataEditorWidget::standardizeDataFormat()
{
    if (!m_dataModel) return;

    for (int col = 0; col < m_dataModel->columnCount(); ++col) {
        for (int row = 0; row < m_dataModel->rowCount(); ++row) {
            QStandardItem* item = m_dataModel->item(row, col);
            if (!item) continue;

            QString text = item->text().trimmed();
            if (text.isEmpty()) continue;

            // 根据列定义标准化格式
            if (col < m_columnDefinitions.size()) {
                const ColumnDefinition& def = m_columnDefinitions[col];

                // 数值类型标准化
                if (def.type == WellTestColumnType::Pressure ||
                    def.type == WellTestColumnType::Temperature ||
                    def.type == WellTestColumnType::FlowRate ||
                    def.type == WellTestColumnType::Time) {

                    bool ok;
                    double value = text.toDouble(&ok);
                    if (ok) {
                        QString formatted = QString::number(value, 'f', def.decimalPlaces);
                        item->setText(formatted);
                    }
                }
            }
        }
    }
}

// ============================================================================
// 保存和导出功能实现
// ============================================================================

bool DataEditorWidget::saveExcelFile(const QString& filePath)
{
    // 简化版本：保存为CSV格式
    QString csvPath = filePath;
    if (filePath.endsWith(".xlsx", Qt::CaseInsensitive) || filePath.endsWith(".xls", Qt::CaseInsensitive)) {
        csvPath = filePath + ".csv";
    }
    return saveCsvFile(csvPath);
}

bool DataEditorWidget::saveCsvFile(const QString& filePath)
{
    if (!m_dataModel) {
        return false;
    }

    QFile file(filePath);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        return false;
    }

    QTextStream out(&file);

#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    out.setCodec("UTF-8");
#else
    out.setEncoding(QStringConverter::Utf8);
#endif

    // 写入表头
    QStringList headers;
    for (int col = 0; col < m_dataModel->columnCount(); ++col) {
        headers.append(m_dataModel->headerData(col, Qt::Horizontal).toString());
    }
    out << headers.join(',') << "\n";

    // 写入数据
    for (int row = 0; row < m_dataModel->rowCount(); ++row) {
        QStringList fields;
        for (int col = 0; col < m_dataModel->columnCount(); ++col) {
            QStandardItem* item = m_dataModel->item(row, col);
            QString text = item ? item->text() : "";

            if (text.contains(',') || text.contains('"') || text.contains('\n')) {
                text = '"' + text.replace('"', "\"\"") + '"';
            }
            fields.append(text);
        }
        out << fields.join(',') << "\n";
    }

    file.close();
    return true;
}

bool DataEditorWidget::saveJsonFile(const QString& filePath)
{
    if (!m_dataModel) {
        return false;
    }

    QJsonArray jsonArray;

    QStringList headers;
    for (int col = 0; col < m_dataModel->columnCount(); ++col) {
        headers.append(m_dataModel->headerData(col, Qt::Horizontal).toString());
    }

    for (int row = 0; row < m_dataModel->rowCount(); ++row) {
        QJsonObject jsonObject;

        for (int col = 0; col < m_dataModel->columnCount(); ++col) {
            QStandardItem* item = m_dataModel->item(row, col);
            QString value = item ? item->text() : "";

            bool isNumber;
            double numValue = value.toDouble(&isNumber);

            if (isNumber) {
                jsonObject[headers[col]] = numValue;
            } else {
                jsonObject[headers[col]] = value;
            }
        }

        jsonArray.append(jsonObject);
    }

    QJsonDocument doc(jsonArray);

    QFile file(filePath);
    if (!file.open(QIODevice::WriteOnly)) {
        return false;
    }

    file.write(doc.toJson());
    file.close();

    return true;
}

bool DataEditorWidget::exportToPdf(const QString& filePath)
{
    if (!m_dataModel) {
        return false;
    }

    QString htmlContent = "<html><head><meta charset='utf-8'>";
    htmlContent += "<style>";
    htmlContent += "body { font-family: 'Microsoft YaHei', Arial, sans-serif; margin: 20px; }";
    htmlContent += "h1 { color: #2c3e50; text-align: center; }";
    htmlContent += "table { border-collapse: collapse; width: 100%; margin: 20px 0; }";
    htmlContent += "th, td { border: 1px solid #e1e8ed; padding: 8px; text-align: left; }";
    htmlContent += "th { background-color: #f8f9fa; font-weight: bold; }";
    htmlContent += "tr:nth-child(even) { background-color: #f9f9f9; }";
    htmlContent += ".stats { margin-top: 20px; padding: 15px; background-color: #e3f2fd; border-radius: 6px; }";
    htmlContent += "</style></head><body>";

    htmlContent += QString("<h1>试井数据报告 - %1</h1>").arg(QFileInfo(m_currentFilePath).baseName());

    // 添加基本信息
    htmlContent += "<div class='stats'>";
    htmlContent += QString("<strong>数据概览：</strong> %1 行 × %2 列<br>")
                       .arg(m_dataModel->rowCount())
                       .arg(m_dataModel->columnCount());
    htmlContent += QString("<strong>生成时间：</strong> %1<br>")
                       .arg(QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss"));
    htmlContent += QString("<strong>文件路径：</strong> %1")
                       .arg(m_currentFilePath);
    htmlContent += "</div>";

    htmlContent += "<table>";

    // 添加表头
    htmlContent += "<thead><tr>";
    for (int col = 0; col < m_dataModel->columnCount(); ++col) {
        QString header = m_dataModel->headerData(col, Qt::Horizontal).toString();
        htmlContent += QString("<th>%1</th>").arg(header.toHtmlEscaped());
    }
    htmlContent += "</tr></thead>";

    // 添加数据行（限制显示行数以避免PDF过大）
    htmlContent += "<tbody>";
    int maxRows = qMin(500, m_dataModel->rowCount());
    for (int row = 0; row < maxRows; ++row) {
        htmlContent += "<tr>";
        for (int col = 0; col < m_dataModel->columnCount(); ++col) {
            QStandardItem* item = m_dataModel->item(row, col);
            QString text = item ? item->text() : "";
            htmlContent += QString("<td>%1</td>").arg(text.toHtmlEscaped());
        }
        htmlContent += "</tr>";
    }
    htmlContent += "</tbody>";

    htmlContent += "</table>";

    if (m_dataModel->rowCount() > maxRows) {
        htmlContent += QString("<p><em>注：为了控制文件大小，仅显示前 %1 行数据。</em></p>").arg(maxRows);
    }

    htmlContent += "</body></html>";

    QTextDocument document;
    document.setHtml(htmlContent);

    QPrinter printer(QPrinter::HighResolution);
    printer.setOutputFormat(QPrinter::PdfFormat);
    printer.setOutputFileName(filePath);
    printer.setPageMargins(QMarginsF(15, 15, 15, 15), QPageLayout::Millimeter);

    document.print(&printer);

    return QFileInfo(filePath).exists();
}

bool DataEditorWidget::exportToHtml(const QString& filePath)
{
    if (!m_dataModel) {
        return false;
    }

    QFile file(filePath);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        return false;
    }

    QTextStream out(&file);
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
    out.setCodec("UTF-8");
#else
    out.setEncoding(QStringConverter::Utf8);
#endif

    out << "<!DOCTYPE html>\n";
    out << "<html lang='zh-CN'>\n";
    out << "<head>\n";
    out << "<meta charset='UTF-8'>\n";
    out << "<meta name='viewport' content='width=device-width, initial-scale=1.0'>\n";
    out << QString("<title>试井数据 - %1</title>\n").arg(QFileInfo(m_currentFilePath).baseName());
    out << "<style>\n";
    out << "body { font-family: 'Microsoft YaHei', Arial, sans-serif; margin: 20px; background-color: #f8f9fa; }\n";
    out << ".container { max-width: 1200px; margin: 0 auto; background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }\n";
    out << "h1 { color: #2c3e50; text-align: center; margin-bottom: 30px; }\n";
    out << "table { border-collapse: collapse; width: 100%; margin-top: 20px; }\n";
    out << "th, td { border: 1px solid #e1e8ed; padding: 10px; text-align: left; }\n";
    out << "th { background: linear-gradient(to bottom, #f8f9fa, #e9ecef); color: #495057; font-weight: 600; }\n";
    out << "tr:nth-child(even) { background-color: #f8f9fa; }\n";
    out << "tr:hover { background-color: #e3f2fd; }\n";
    out << ".stats { margin-bottom: 20px; padding: 15px; background-color: #e3f2fd; border-radius: 6px; }\n";
    out << "</style>\n";
    out << "</head>\n";
    out << "<body>\n";
    out << "<div class='container'>\n";
    out << QString("<h1>试井数据 - %1</h1>\n").arg(QFileInfo(m_currentFilePath).baseName());

    // 统计信息
    out << "<div class='stats'>\n";
    out << QString("<strong>数据概览：</strong> %1 行 × %2 列 | ")
               .arg(m_dataModel->rowCount())
               .arg(m_dataModel->columnCount());
    out << QString("<strong>生成时间：</strong> %1<br>")
               .arg(QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss"));
    out << QString("<strong>文件路径：</strong> %1")
               .arg(m_currentFilePath);
    out << "</div>\n";

    out << "<table>\n";

    // 表头
    out << "<thead><tr>\n";
    for (int col = 0; col < m_dataModel->columnCount(); ++col) {
        QString header = m_dataModel->headerData(col, Qt::Horizontal).toString();
        out << QString("<th>%1</th>\n").arg(header.toHtmlEscaped());
    }
    out << "</tr></thead>\n";

    // 数据
    out << "<tbody>\n";
    for (int row = 0; row < m_dataModel->rowCount(); ++row) {
        out << "<tr>\n";
        for (int col = 0; col < m_dataModel->columnCount(); ++col) {
            QStandardItem* item = m_dataModel->item(row, col);
            QString text = item ? item->text() : "";
            out << QString("<td>%1</td>\n").arg(text.toHtmlEscaped());
        }
        out << "</tr>\n";
    }
    out << "</tbody>\n";

    out << "</table>\n";
    out << "</div>\n";
    out << "</body>\n";
    out << "</html>\n";

    file.close();
    return true;
}

// ============================================================================
// 数据验证功能实现
// ============================================================================

ValidationResult DataEditorWidget::validateData() const
{
    ValidationResult result;
    result.isValid = true;
    result.totalRows = m_dataModel ? m_dataModel->rowCount() : 0;
    result.validRows = 0;
    result.errorRows = 0;

    if (!m_dataModel) {
        result.errors.append("没有加载数据");
        result.isValid = false;
        return result;
    }

    // 验证每一行
    for (int row = 0; row < m_dataModel->rowCount(); ++row) {
        bool hasError = false;
        bool isEmpty = true;

        for (int col = 0; col < m_dataModel->columnCount(); ++col) {
            QStandardItem* item = m_dataModel->item(row, col);
            QString value = item ? item->text().trimmed() : "";

            if (!value.isEmpty()) {
                isEmpty = false;

                // 验证列定义
                if (col < m_columnDefinitions.size()) {
                    QStringList columnErrors;
                    if (!validateColumnData(col, m_columnDefinitions[col], columnErrors)) {
                        hasError = true;
                        QString columnName = m_columnDefinitions[col].name;
                        result.columnErrors[columnName].append(columnErrors);
                    }
                }
            }
        }

        if (isEmpty) {
            result.warnings.append(QString("第%1行为空行").arg(row + 1));
        } else if (hasError) {
            result.errorRows++;
        } else {
            result.validRows++;
        }
    }

    result.isValid = result.errors.isEmpty();
    return result;
}

bool DataEditorWidget::validateColumnData(int columnIndex, const ColumnDefinition& definition, QStringList& errors) const
{
    if (!m_dataModel || columnIndex < 0 || columnIndex >= m_dataModel->columnCount()) {
        return false;
    }

    bool isValid = true;
    int emptyCount = 0;

    for (int row = 0; row < m_dataModel->rowCount(); ++row) {
        QStandardItem* item = m_dataModel->item(row, columnIndex);
        QString value = item ? item->text().trimmed() : "";

        if (value.isEmpty()) {
            emptyCount++;
            if (definition.isRequired) {
                errors.append(QString("第%1行缺少必需数据").arg(row + 1));
                isValid = false;
            }
            continue;
        }

        // 数值范围验证
        if (definition.type == WellTestColumnType::Pressure ||
            definition.type == WellTestColumnType::Temperature ||
            definition.type == WellTestColumnType::FlowRate ||
            definition.type == WellTestColumnType::Time) {

            bool ok;
            double numValue = value.toDouble(&ok);

            if (!ok) {
                errors.append(QString("第%1行数据格式错误，应为数值").arg(row + 1));
                isValid = false;
            } else {
                if (numValue < definition.minValue || numValue > definition.maxValue) {
                    errors.append(QString("第%1行数值超出范围 [%2, %3]")
                                      .arg(row + 1)
                                      .arg(definition.minValue)
                                      .arg(definition.maxValue));
                    isValid = false;
                }
            }
        }
    }

    // 检查必需列的数据完整性
    if (definition.isRequired && emptyCount > m_dataModel->rowCount() * 0.5) {
        errors.append(QString("必需列'%1'有超过50%的数据缺失").arg(definition.name));
        isValid = false;
    }

    return isValid;
}

bool DataEditorWidget::isNumericData(const QString& data) const
{
    bool ok;
    data.toDouble(&ok);
    return ok;
}

bool DataEditorWidget::isDateTimeData(const QString& data) const
{
    QStringList formats = {
        "yyyy-MM-dd", "yyyy/MM/dd", "dd/MM/yyyy", "dd-MM-yyyy",
        "yyyy-MM-dd hh:mm:ss", "yyyy/MM/dd hh:mm:ss",
        "dd/MM/yyyy hh:mm:ss", "dd-MM-yyyy hh:mm:ss"
    };

    for (const QString& format : formats) {
        QDateTime dt = QDateTime::fromString(data, format);
        if (dt.isValid()) {
            return true;
        }
    }

    return false;
}

QStringList DataEditorWidget::detectDataType(int column) const
{
    if (!m_dataModel || column < 0 || column >= m_dataModel->columnCount()) {
        return QStringList() << "未知";
    }

    QSet<QString> types;

    for (int row = 0; row < m_dataModel->rowCount(); ++row) {
        QStandardItem* item = m_dataModel->item(row, column);
        QString value = item ? item->text().trimmed() : "";

        if (value.isEmpty()) {
            continue;
        }

        if (isNumericData(value)) {
            types.insert("数值型");
        } else if (isDateTimeData(value)) {
            types.insert("日期时间型");
        } else {
            types.insert("文本型");
        }
    }

    if (types.size() > 1) {
        return QStringList() << "混合型";
    } else if (types.size() == 1) {
        QStringList result;
        for (const QString& type : types) {
            result.append(type);
        }
        return result;
    } else {
        return QStringList() << "空";
    }
}

// ============================================================================
// 列定义管理功能实现
// ============================================================================

void DataEditorWidget::setColumnDefinitions(const QList<ColumnDefinition>& definitions)
{
    m_columnDefinitions = definitions;

    // 应用列定义
    for (int i = 0; i < m_columnDefinitions.size() && i < m_dataModel->columnCount(); ++i) {
        applyColumnDefinition(i, m_columnDefinitions[i]);
    }

    emit columnDefinitionsChanged();
}

QList<ColumnDefinition> DataEditorWidget::getColumnDefinitions() const
{
    return m_columnDefinitions;
}



void DataEditorWidget::applyColumnDefinition(int columnIndex, const ColumnDefinition& definition)
{
    if (!m_dataModel || columnIndex < 0 || columnIndex >= m_dataModel->columnCount()) {
        return;
    }

    // 应用数据格式
    for (int row = 0; row < m_dataModel->rowCount(); ++row) {
        QStandardItem* item = m_dataModel->item(row, columnIndex);
        if (!item) continue;

        // 根据类型格式化数值
        if (definition.type == WellTestColumnType::Pressure ||
            definition.type == WellTestColumnType::Temperature ||
            definition.type == WellTestColumnType::FlowRate ||
            definition.type == WellTestColumnType::Time) {

            bool ok;
            double value = item->text().toDouble(&ok);
            if (ok) {
                QString formatted = QString::number(value, 'f', definition.decimalPlaces);
                item->setText(formatted);
            }
        }

        // 设置颜色标记
        if (definition.isRequired) {
            item->setBackground(QBrush(QColor("#fff3cd"))); // 淡黄色背景表示必需
        }
    }
}

// ============================================================================
// 撤销重做功能实现
// ============================================================================

void DataEditorWidget::undo()
{
    if (m_undoStack && m_undoStack->canUndo()) {
        m_undoStack->undo();
        m_dataModified = true;
        updateStatus("已撤销操作", "info");
        updateDataInfo();
        emitDataChanged();
    }
}

void DataEditorWidget::redo()
{
    if (m_undoStack && m_undoStack->canRedo()) {
        m_undoStack->redo();
        m_dataModified = true;
        updateStatus("已重做操作", "info");
        updateDataInfo();
        emitDataChanged();
    }
}

bool DataEditorWidget::canUndo() const
{
    return m_undoStack && m_undoStack->canUndo();
}

bool DataEditorWidget::canRedo() const
{
    return m_undoStack && m_undoStack->canRedo();
}

// ============================================================================
// 数据模型变化处理
// ============================================================================

void DataEditorWidget::onCellDataChanged(QStandardItem* item)
{
    if (!item || !m_undoStack) {
        return;
    }

    m_dataModified = true;
    updateStatus("数据已修改", "warning");
    updateDataInfo();
    emitDataChanged();
}

void DataEditorWidget::onModelDataChanged(const QModelIndex& topLeft, const QModelIndex& bottomRight)
{
    Q_UNUSED(topLeft)
    Q_UNUSED(bottomRight)

    m_dataModified = true;
    updateStatus("数据已修改", "warning");
    updateDataInfo();
    emitDataChanged();
}

// ============================================================================
// UI更新和辅助方法实现
// ============================================================================

void DataEditorWidget::updateStatus(const QString& message, const QString& type)
{
    ui->statusLabel->setText(message);

    QString indicatorStyle;
    if (type == "success") {
        indicatorStyle = "background-color: #28a745;";
    } else if (type == "warning") {
        indicatorStyle = "background-color: #fd7e14;";
    } else if (type == "error") {
        indicatorStyle = "background-color: #fd7e14;";
    } else {
        indicatorStyle = "background-color: #4a90e2;";
    }

    ui->statusIndicator->setStyleSheet(QString("QLabel { %1 border-radius: 5px; }").arg(indicatorStyle));
}

void DataEditorWidget::updateDataInfo()
{
    if (!m_dataModel) {
        ui->dataInfoLabel->setText("无数据");
        return;
    }

    QString info = QString("%1行 × %2列")
                       .arg(m_dataModel->rowCount())
                       .arg(m_dataModel->columnCount());

    if (m_largeFileMode) {
        info += " (大文件模式)";
    }

    if (m_dataModified) {
        info += " *";
    }

    ui->dataInfoLabel->setText(info);
}

void DataEditorWidget::setButtonsEnabled(bool enabled)
{
    ui->btnSave->setEnabled(enabled);
    ui->btnExport->setEnabled(enabled);
    ui->btnDefineColumns->setEnabled(enabled);
    ui->btnTimeConvert->setEnabled(enabled);
    ui->btnPressureDropCalc->setEnabled(enabled);
    ui->btnPressureDerivativeCalc->setEnabled(enabled);
    ui->btnDataClean->setEnabled(enabled);
    ui->btnDataStatistics->setEnabled(enabled);
}

void DataEditorWidget::showAnimatedProgress(const QString& title, const QString& message)
{
    if (!m_progressDialog) {
        m_progressDialog = new AnimatedProgressDialog(title, message, this);
    } else {
        m_progressDialog->setWindowTitle(title);
        m_progressDialog->setMessage(message);
    }

    m_progressDialog->show();
    QApplication::processEvents();
}

void DataEditorWidget::hideAnimatedProgress()
{
    if (m_progressDialog) {
        m_progressDialog->hide();
    }
}

void DataEditorWidget::updateProgress(int value, const QString& message)
{
    if (m_progressDialog) {
        m_progressDialog->setProgress(value);
        if (!message.isEmpty()) {
            m_progressDialog->setMessage(message);
        }
        QApplication::processEvents();
    }
}

void DataEditorWidget::clearData()
{
    if (m_dataModel) {
        m_dataModel->clear();
    }

    if (m_undoStack) {
        m_undoStack->clear();
    }

    m_currentFilePath = "";
    m_currentFileType = "";
    ui->filePathLineEdit->clear();
    ui->searchLineEdit->clear();

    clearDataFilter();

    m_columnDefinitions.clear();
    updateStatus("就绪", "success");
    setButtonsEnabled(false);

    m_dataModified = false;
    m_largeFileMode = false;

    updateDataInfo();
    emitDataChanged();
}

void DataEditorWidget::applyColumnStyles()
{
    if (!m_dataModel) return;

    QColor textColor("#2c3e50");

    for (int row = 0; row < m_dataModel->rowCount(); ++row) {
        for (int col = 0; col < m_dataModel->columnCount(); ++col) {
            QStandardItem* item = m_dataModel->item(row, col);
            if (item) {
                item->setForeground(QBrush(textColor));
            }
        }
    }
}

void DataEditorWidget::optimizeColumnWidths()
{
    if (!ui->dataTableView || !m_dataModel) {
        return;
    }

    try {
        ui->dataTableView->resizeColumnsToContents();

        QHeaderView* header = ui->dataTableView->horizontalHeader();
        if (!header) {
            return;
        }

        header->setDefaultSectionSize(100); // 设置默认列宽
        header->setMinimumSectionSize(60);   // 设置最小列宽

        for (int i = 0; i < header->count(); ++i) {
            try {
                int width = header->sectionSize(i);
                if (width > 200) {
                    header->resizeSection(i, 200);  // 设置最大列宽
                } else if (width < 80) {
                    header->resizeSection(i, 80);   // 设置最小列宽
                }
            } catch (...) {
                qDebug() << "调整第" << i << "列宽度时出错";
                continue;
            }
        }
    } catch (...) {
        qDebug() << "优化列宽时发生异常";
    }
}

void DataEditorWidget::optimizeTableDisplay()
{
    if (!ui->dataTableView) return;

    // 设置表格的行高和列宽
    QHeaderView* verticalHeader = ui->dataTableView->verticalHeader();
    QHeaderView* horizontalHeader = ui->dataTableView->horizontalHeader();

    // 设置行高 - 调整为合适的大小
    verticalHeader->setDefaultSectionSize(24);  // 设置默认行高
    verticalHeader->setMinimumSectionSize(20);   // 设置最小行高

    // 设置列宽
    horizontalHeader->setDefaultSectionSize(100); // 设置默认列宽
    horizontalHeader->setMinimumSectionSize(60); // 设置最小列宽

    // 设置表格的其他显示属性
    ui->dataTableView->setAlternatingRowColors(true);
    ui->dataTableView->setShowGrid(true);
    ui->dataTableView->setGridStyle(Qt::SolidLine);

    // 确保表格内容能够完全显示
    ui->dataTableView->setHorizontalScrollMode(QAbstractItemView::ScrollPerPixel);
    ui->dataTableView->setVerticalScrollMode(QAbstractItemView::ScrollPerPixel);

    // 设置选择模式
    ui->dataTableView->setSelectionBehavior(QAbstractItemView::SelectItems);
    ui->dataTableView->setSelectionMode(QAbstractItemView::ExtendedSelection);
}

// ============================================================================
// 选择和交互方法实现
// ============================================================================

int DataEditorWidget::getSelectedRow() const
{
    QModelIndexList selection = ui->dataTableView->selectionModel()->selectedIndexes();
    if (selection.isEmpty()) {
        return -1;
    }
    QModelIndex sourceIndex = m_proxyModel->mapToSource(selection.first());
    return sourceIndex.row();
}

int DataEditorWidget::getSelectedColumn() const
{
    QModelIndexList selection = ui->dataTableView->selectionModel()->selectedIndexes();
    if (selection.isEmpty()) {
        return -1;
    }
    QModelIndex sourceIndex = m_proxyModel->mapToSource(selection.first());
    return sourceIndex.column();
}

QList<int> DataEditorWidget::getSelectedRows() const
{
    QList<int> rows;
    QModelIndexList selection = ui->dataTableView->selectionModel()->selectedIndexes();

    for (const QModelIndex& index : selection) {
        QModelIndex sourceIndex = m_proxyModel->mapToSource(index);
        int row = sourceIndex.row();
        if (!rows.contains(row)) {
            rows.append(row);
        }
    }

    return rows;
}

QList<int> DataEditorWidget::getSelectedColumns() const
{
    QList<int> columns;
    QModelIndexList selection = ui->dataTableView->selectionModel()->selectedIndexes();

    for (const QModelIndex& index : selection) {
        QModelIndex sourceIndex = m_proxyModel->mapToSource(index);
        int column = sourceIndex.column();
        if (!columns.contains(column)) {
            columns.append(column);
        }
    }

    return columns;
}

bool DataEditorWidget::checkDataModifiedAndPrompt()
{
    if (m_dataModified) {
        QMessageBox msgBox;
        msgBox.setWindowTitle("保存更改");
        msgBox.setText("当前数据已被修改，是否保存更改？");
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No | QMessageBox::Cancel);
        msgBox.setDefaultButton(QMessageBox::Yes);

        int result = msgBox.exec();

        if (result == QMessageBox::Yes) {
            onSave();
            return true;
        } else if (result == QMessageBox::No) {
            return true;
        } else {
            return false;
        }
    }

    return true;
}

// ============================================================================
// 辅助方法实现
// ============================================================================

void DataEditorWidget::emitDataChanged()
{
    // 使用更安全的方式发射信号
    if (this && !isVisible()) {
        // 如果窗口不可见，延迟发射信号
        QTimer::singleShot(500, this, [this]() {
            if (this) {
                try {
                    emit dataChanged();
                } catch (...) {
                    qDebug() << "延迟发射dataChanged信号时出错";
                }
            }
        });
    } else {
        // 直接发射信号
        try {
            emit dataChanged();
        } catch (...) {
            qDebug() << "发射dataChanged信号时出错";
        }
    }
}

QString DataEditorWidget::formatNumber(double number, int precision) const
{
    return QString::number(number, 'f', precision);
}

void DataEditorWidget::showStyledMessageBox(const QString& title, const QString& text,
                                            QMessageBox::Icon icon, const QString& detailedText)
{
    QMessageBox msgBox;
    msgBox.setWindowTitle(title);
    msgBox.setText(text);
    msgBox.setIcon(icon);

    if (!detailedText.isEmpty()) {
        msgBox.setDetailedText(detailedText);
    }

    msgBox.setStyleSheet(R"(
        QMessageBox {
            background-color: #ffffff;
            color: #2c3e50;
            font-family: "Microsoft YaHei", "微软雅黑", Arial, sans-serif;
        }
        QMessageBox QLabel {
            color: #2c3e50;
            font-size: 13px;
            padding: 10px;
        }
        QMessageBox QPushButton {
            background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                        stop:0 #4a90e2, stop:1 #357abd);
            color: white;
            border: none;
            border-radius: 6px;
            padding: 8px 20px;
            font-weight: bold;
            min-width: 80px;
            font-size: 12px;
        }
        QMessageBox QPushButton:hover {
            background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                        stop:0 #357abd, stop:1 #2a628a);
        }
        QMessageBox QPushButton:pressed {
            background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                        stop:0 #2a628a, stop:1 #1e4a6b);
        }
    )");

    msgBox.exec();
}

// ============================================================================
// 压力导数计算功能实现（简化版）
// ============================================================================

// 设置压力导数计算器
void DataEditorWidget::setupPressureDerivativeCalculator()
{
    m_pressureDerivativeCalculator = new PressureDerivativeCalculator(this);

    // 连接进度信号
    connect(m_pressureDerivativeCalculator, &PressureDerivativeCalculator::progressUpdated,
            this, [this](int progress, const QString& message) {
                if (m_progressDialog) {
                    m_progressDialog->setProgress(progress);
                    m_progressDialog->setMessage(message);
                    QApplication::processEvents();
                }
            });

    // 连接计算完成信号
    connect(m_pressureDerivativeCalculator, &PressureDerivativeCalculator::calculationCompleted,
            this, [this](const PressureDerivativeResult& result) {
                emit pressureDerivativeCalculated(result);
            });
}

// 压力导数计算槽函数
void DataEditorWidget::onPressureDerivativeCalc()
{
    if (!hasData()) {
        showStyledMessageBox("压力导数计算", "请先加载数据文件", QMessageBox::Information);
        return;
    }

    if (m_dataModel->rowCount() < 3) {
        showStyledMessageBox("压力导数计算", "数据行数不足（至少需要3行数据）", QMessageBox::Warning);
        return;
    }

    // 自动检测列
    PressureDerivativeConfig config = m_pressureDerivativeCalculator->autoDetectColumns(m_dataModel);

    // 检查是否找到压力列
    if (config.pressureColumnIndex == -1) {
        showStyledMessageBox("压力导数计算", "未找到压力列，请确保数据中包含压力数据列", QMessageBox::Warning);
        return;
    }

    // 检查是否找到时间列
    if (config.timeColumnIndex == -1) {
        showStyledMessageBox("压力导数计算", "未找到时间列，请确保数据中包含时间数据列", QMessageBox::Warning);
        return;
    }

    // 获取压力单位
    QStandardItem* pressureHeader = m_dataModel->horizontalHeaderItem(config.pressureColumnIndex);
    if (pressureHeader) {
        QString headerText = pressureHeader->text();
        if (headerText.contains("MPa")) {
            config.pressureUnit = "MPa";
        } else if (headerText.contains("kPa")) {
            config.pressureUnit = "kPa";
        } else if (headerText.contains("psi")) {
            config.pressureUnit = "psi";
        } else {
            config.pressureUnit = "MPa";
        }
    }

    // 显示计算进度
    showAnimatedProgress("压力导数计算", "正在计算压力导数...");

    // 执行计算
    PressureDerivativeResult result = m_pressureDerivativeCalculator->calculatePressureDerivative(m_dataModel, config);

    hideAnimatedProgress();

    if (result.success) {
        updateStatus(QString("压力导数计算完成 - 已添加列: %1").arg(result.columnName), "success");
        m_dataModified = true;
        emitDataChanged();

        // 添加列定义
        ColumnDefinition newColumnDef;
        newColumnDef.name = result.columnName;
        newColumnDef.type = WellTestColumnType::PressureDerivative;
        newColumnDef.unit = config.pressureUnit;
        newColumnDef.description = "压力导数";
        newColumnDef.isRequired = false;
        newColumnDef.minValue = -999999;
        newColumnDef.maxValue = 999999;
        newColumnDef.decimalPlaces = 6;

        if (result.addedColumnIndex < m_columnDefinitions.size()) {
            m_columnDefinitions.insert(result.addedColumnIndex, newColumnDef);
        } else {
            m_columnDefinitions.append(newColumnDef);
        }

        showStyledMessageBox("压力导数计算完成",
                             QString("压力导数计算成功完成！\n"
                                     "新增列：%1\n"
                                     "处理行数：%2")
                                 .arg(result.columnName)
                                 .arg(result.processedRows),
                             QMessageBox::Information);
    } else {
        updateStatus("压力导数计算失败", "error");
        showStyledMessageBox("压力导数计算失败", result.errorMessage, QMessageBox::Warning);
    }
}

// 使用配置计算压力导数
PressureDerivativeResult DataEditorWidget::calculatePressureDerivativeWithConfig(const PressureDerivativeConfig& config)
{
    if (!m_pressureDerivativeCalculator) {
        PressureDerivativeResult result;
        result.success = false;
        result.errorMessage = "压力导数计算器未初始化";
        return result;
    }

    return m_pressureDerivativeCalculator->calculatePressureDerivative(m_dataModel, config);
}

// 获取默认压力导数配置
PressureDerivativeConfig DataEditorWidget::getDefaultPressureDerivativeConfig()
{
    PressureDerivativeConfig config;
    config.timeUnit = "h";
    config.pressureUnit = "MPa";

    if (m_pressureDerivativeCalculator && m_dataModel) {
        PressureDerivativeConfig autoConfig = m_pressureDerivativeCalculator->autoDetectColumns(m_dataModel);
        config.pressureColumnIndex = autoConfig.pressureColumnIndex;
        config.timeColumnIndex = autoConfig.timeColumnIndex;
    }

    return config;
}
