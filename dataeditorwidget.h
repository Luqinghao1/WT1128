#ifndef DATAEDITORWIDGET_H
#define DATAEDITORWIDGET_H

#include <QWidget>
#include <QString>
#include <QTableView>
#include <QStandardItemModel>
#include <QFile>
#include <QInputDialog>
#include <QMessageBox>
#include <QFileDialog>
#include <QSortFilterProxyModel>
#include <QUndoStack>
#include <QUndoCommand>
#include <QTimer>
#include <QProgressDialog>
#include <QTextDocument>
#include <QTextDocumentWriter>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QSet>
#include <QList>
#include <QDialog>
#include <QComboBox>
#include <QSpinBox>
#include <QDateTimeEdit>
#include <QCheckBox>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QFormLayout>
#include <QGroupBox>
#include <QLabel>
#include <QPushButton>
#include <QProgressBar>
#include <QMovie>
#include <QPropertyAnimation>
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
#include <QMenu>
#include <QAction>
#include <QPoint>
#include <QTime>
#include <QDate>
#include <QDateTime>

// Qt6兼容性处理
#if QT_VERSION >= QT_VERSION_CHECK(6, 0, 0)
#include <QStringConverter>
#else
#include <QTextCodec>
#endif

// 新增：压力导数计算器头文件
#include "PressureDerivativeCalculator.h"

namespace Ui {
class DataEditorWidget;
}

// 试井数据列类型枚举 - 新增日期、时刻类型和压力导数类型
enum class WellTestColumnType {
    SerialNumber,       // 序号
    Date,               // 日期 - 新增
    Time,               // 时间
    TimeOfDay,          // 时刻 - 新增
    Pressure,           // 压力
    Temperature,        // 温度
    FlowRate,           // 流量
    Depth,              // 深度
    Viscosity,          // 粘度
    Density,            // 密度
    Permeability,       // 渗透率
    Porosity,           // 孔隙度
    WellRadius,         // 井半径
    SkinFactor,         // 表皮系数
    Distance,           // 距离
    Volume,             // 体积
    PressureDrop,       // 压降
    PressureDerivative, // 压力导数 - 新增
    Custom              // 自定义
};

// 单位系统枚举
enum class UnitSystem {
    Metric,         // 公制
    Imperial,       // 英制
    Oilfield        // 油田单位
};

// 数据列定义结构
struct ColumnDefinition {
    QString name;                           // 列名
    WellTestColumnType type;               // 数据类型
    QString unit;                          // 单位
    QString description;                   // 描述
    bool isRequired;                       // 是否必需
    double minValue;                       // 最小值
    double maxValue;                       // 最大值
    int decimalPlaces;                     // 小数位数

    ColumnDefinition() :
        type(WellTestColumnType::Custom),
        isRequired(false),
        minValue(-999999),
        maxValue(999999),
        decimalPlaces(3) {}
};

// 数据统计结构
struct DataStatistics {
    QString columnName;
    int dataCount;
    int validCount;
    int invalidCount;
    double minimum;
    double maximum;
    double average;
    double median;
    double standardDeviation;
    QString dataType;
    QString unit;
};

// 数据验证结果结构
struct ValidationResult {
    bool isValid;
    QStringList errors;
    QStringList warnings;
    int totalRows;
    int validRows;
    int errorRows;
    QMap<QString, QStringList> columnErrors;
};

// 修改后的时间转换配置结构
struct TimeConversionConfig {
    int dateColumnIndex;         // 日期列索引（-1表示不使用日期列）
    int timeColumnIndex;         // 时刻列索引（-1表示不使用时刻列）
    int sourceTimeColumnIndex;   // 源时间列索引（用于兼容旧版本）
    QString outputUnit;          // 输出单位（s, m, h）
    QString newColumnName;       // 新列名称
    bool useDateAndTime;         // 是否使用日期+时刻模式
};


// 时间转换结果结构
struct TimeConversionResult {
    bool success;
    QString errorMessage;
    int addedColumnIndex;
    QString columnName;
    int processedRows;
};

// 压降计算结果结构
struct PressureDropResult {
    bool success;
    QString errorMessage;
    int addedColumnIndex;
    QString columnName;
    int processedRows;
};

// 撤销命令基类
class DataEditCommand : public QUndoCommand
{
public:
    DataEditCommand(QStandardItemModel* model, QUndoCommand* parent = nullptr);
    virtual ~DataEditCommand() = default;

protected:
    QStandardItemModel* m_model;
};

// 单元格编辑命令
class CellEditCommand : public DataEditCommand
{
public:
    CellEditCommand(QStandardItemModel* model, int row, int column,
                    const QString& oldValue, const QString& newValue,
                    QUndoCommand* parent = nullptr);
    void undo() override;
    void redo() override;

private:
    int m_row;
    int m_column;
    QString m_oldValue;
    QString m_newValue;
};

// 行操作命令
class RowEditCommand : public DataEditCommand
{
public:
    enum Operation { Insert, Delete };

    RowEditCommand(QStandardItemModel* model, Operation op, int row,
                   const QStringList& rowData = QStringList(),
                   QUndoCommand* parent = nullptr);
    void undo() override;
    void redo() override;

private:
    Operation m_operation;
    int m_row;
    QStringList m_rowData;
};

// 列操作命令
class ColumnEditCommand : public DataEditCommand
{
public:
    enum Operation { Insert, Delete };

    ColumnEditCommand(QStandardItemModel* model, Operation op, int column,
                      const QString& headerName = QString(),
                      const QStringList& columnData = QStringList(),
                      QUndoCommand* parent = nullptr);
    void undo() override;
    void redo() override;

private:
    Operation m_operation;
    int m_column;
    QString m_headerName;
    QStringList m_columnData;
};

// 数据读取配置对话框
class DataLoadConfigDialog : public QDialog
{
    Q_OBJECT

public:
    struct LoadConfig {
        int startRow;           // 开始读取的行号（从1开始）
        bool hasHeader;         // 是否有表头
        QString encoding;       // 编码格式
        QString separator;      // 分隔符
    };

    explicit DataLoadConfigDialog(const QString& filePath, QWidget* parent = nullptr);
    LoadConfig getLoadConfig() const;

private slots:
    void onPreviewClicked();
    void onHasHeaderChanged(bool hasHeader);

private:
    void setupUI();
    void loadFilePreview();
    QString detectEncoding(const QString& filePath);
    QString detectSeparator(const QString& filePath);

    QString m_filePath;
    LoadConfig m_config;

    QSpinBox* m_startRowSpin;
    QCheckBox* m_hasHeaderCheck;
    QComboBox* m_encodingCombo;
    QComboBox* m_separatorCombo;
    QTextEdit* m_previewText;
    QPushButton* m_previewButton;
};

// 列定义对话框
class ColumnDefinitionDialog : public QDialog
{
    Q_OBJECT

public:
    explicit ColumnDefinitionDialog(const QStringList& columnNames,
                                    const QList<ColumnDefinition>& definitions = QList<ColumnDefinition>(),
                                    QWidget* parent = nullptr);

    QList<ColumnDefinition> getColumnDefinitions() const;

private slots:
    void onTypeChanged(int index);
    void onUnitChanged(int index);
    void onCustomTypeTextChanged(const QString& text);
    void onCustomTypeChanged(const QString& text);
    void onCustomUnitChanged(const QString& text);
    void onResetClicked();

private:
    void setupUI();
    void updateUnitsForType(WellTestColumnType type, QComboBox* unitCombo);
    void updatePreviewLabel(int index);
    QString getDefaultUnit(WellTestColumnType type, UnitSystem system = UnitSystem::Metric);
    void loadPresetDefinitions();

    QStringList m_columnNames;
    QList<ColumnDefinition> m_definitions;
    QList<QComboBox*> m_typeComboBoxes;
    QList<QComboBox*> m_unitComboBoxes;
    QList<QLineEdit*> m_customTypeEdits;
    QList<QLineEdit*> m_customUnitEdits;
    QList<QCheckBox*> m_requiredChecks;
    QList<QLabel*> m_previewLabels;
    QVBoxLayout* m_mainLayout;
    QFormLayout* m_formLayout;
};

// 修改后的时间转换对话框
class TimeConversionDialog : public QDialog
{
    Q_OBJECT

public:
    explicit TimeConversionDialog(const QStringList& columnNames, QWidget* parent = nullptr);

    TimeConversionConfig getConversionConfig() const;

private slots:
    void onPreviewClicked();
    void onConversionModeChanged();

private:
    void setupUI();
    void updateUIForMode();
    QString previewConversion(const QString& sampleDateInput, const QString& sampleTimeInput, const QString& unit);

    QStringList m_columnNames;
    TimeConversionConfig m_config;

    QRadioButton* m_dateTimeRadio;
    QRadioButton* m_timeOnlyRadio;
    QComboBox* m_dateColumnCombo;
    QComboBox* m_timeColumnCombo;
    QComboBox* m_sourceColumnCombo;
    QComboBox* m_outputUnitCombo;
    QLineEdit* m_newColumnNameEdit;
    QLabel* m_previewLabel;
    QPushButton* m_previewButton;
};

// 数据清理对话框
class DataCleaningDialog : public QDialog
{
    Q_OBJECT

public:
    explicit DataCleaningDialog(QWidget* parent = nullptr);

    struct CleaningOptions {
        bool removeEmptyRows;
        bool removeEmptyColumns;
        bool removeDuplicates;
        bool fillMissingValues;
        bool removeOutliers;
        bool standardizeFormat;
        QString fillMethod;  // "zero", "interpolation", "average"
        double outlierThreshold;
    };

    CleaningOptions getCleaningOptions() const;

private:
    void setupUI();

    QCheckBox* m_removeEmptyRowsCheck;
    QCheckBox* m_removeEmptyColumnsCheck;
    QCheckBox* m_removeDuplicatesCheck;
    QCheckBox* m_fillMissingValuesCheck;
    QCheckBox* m_removeOutliersCheck;
    QCheckBox* m_standardizeFormatCheck;
    QComboBox* m_fillMethodCombo;
    QSpinBox* m_outlierThresholdSpin;
};

// 动画进度对话框
class AnimatedProgressDialog : public QDialog
{
    Q_OBJECT

public:
    explicit AnimatedProgressDialog(const QString& title, const QString& message, QWidget* parent = nullptr);

    void setProgress(int value);
    void setMessage(const QString& message);
    void setMaximum(int maximum);

protected:
    void closeEvent(QCloseEvent* event) override;

private:
    void setupUI();
    void setupAnimation();

    QLabel* m_iconLabel;
    QLabel* m_messageLabel;
    QProgressBar* m_progressBar;
    QMovie* m_loadingMovie;
    QPropertyAnimation* m_fadeAnimation;
    QGraphicsOpacityEffect* m_opacityEffect;
};

// ============================================================================
// 主数据编辑器类
// ============================================================================
class DataEditorWidget : public QWidget
{
    Q_OBJECT

public:
    explicit DataEditorWidget(QWidget *parent = nullptr);
    ~DataEditorWidget();

    // 加载并显示数据
    void loadData(const QString& filePath, const QString& fileType);
    void loadDataWithConfig(const QString& filePath, const QString& fileType, const DataLoadConfigDialog::LoadConfig& config);

    // 获取数据模型和文件信息的方法
    QStandardItemModel* getDataModel() const { return m_dataModel; }
    QString getCurrentFileName() const { return m_currentFilePath; }
    QString getCurrentFileType() const { return m_currentFileType; }
    bool hasData() const { return m_dataModel && m_dataModel->rowCount() > 0 && m_dataModel->columnCount() > 0; }

    // 数据处理功能
    DataStatistics calculateColumnStatistics(int column) const;
    QList<DataStatistics> calculateAllStatistics() const;
    ValidationResult validateData() const;
    void applyDataFilter(const QString& filterText);
    void clearDataFilter();

    // 撤销重做功能
    void undo();
    void redo();
    bool canUndo() const;
    bool canRedo() const;

    // 列定义管理
    void setColumnDefinitions(const QList<ColumnDefinition>& definitions);
    QList<ColumnDefinition> getColumnDefinitions() const;

    // 修改后的时间转换功能
    TimeConversionResult convertTimeColumn(const TimeConversionConfig& config);

    // 压降计算功能
    PressureDropResult calculatePressureDrop();

    // 新增：压力导数计算功能
    PressureDerivativeResult calculatePressureDerivativeWithConfig(const PressureDerivativeConfig& config);
    PressureDerivativeConfig getDefaultPressureDerivativeConfig();

    // 列标题更新功能
    void updateColumnHeaders();

signals:
    // 文件更换信号
    void fileChanged(const QString& filePath, const QString& fileType);

    // 数据变化信号
    void dataChanged();

    // 数据处理信号
    void statisticsCalculated(const QList<DataStatistics>& statistics);
    void dataValidated(const ValidationResult& result);
    void searchCompleted(int matchCount);
    void columnDefinitionsChanged();
    void timeConversionCompleted(const TimeConversionResult& result);
    void pressureDropCalculated(const PressureDropResult& result);

    // 新增：压力导数计算完成信号
    void pressureDerivativeCalculated(const PressureDerivativeResult& result);

private slots:
    // 文件操作槽函数
    void onOpenFile();
    void onSave();
    void onExport();

    // 数据处理槽函数
    void onDefineColumns();
    void onTimeConvert();
    void onDataClean();
    void onDataStatistics();
    void onPressureDropCalc();

    // 新增：压力导数计算槽函数
    void onPressureDerivativeCalc();

    // 搜索槽函数
    void onSearchTextChanged();
    void onSearchData();

    // 模型数据变化槽函数
    void onCellDataChanged(QStandardItem* item);
    void onModelDataChanged(const QModelIndex& topLeft, const QModelIndex& bottomRight);

    // 右键菜单槽函数
    void onTableContextMenuRequested(const QPoint& pos);
    void onAddRowAbove();
    void onAddRowBelow();
    void onDeleteSelectedRows();
    void onAddColumnLeft();
    void onAddColumnRight();
    void onDeleteSelectedColumns();

private:
    Ui::DataEditorWidget *ui;

    // 数据模型和代理
    QStandardItemModel* m_dataModel;
    QSortFilterProxyModel* m_proxyModel;

    // 撤销重做栈
    QUndoStack* m_undoStack;

    // 当前加载的文件信息
    QString m_currentFilePath;
    QString m_currentFileType;

    // 数据状态
    bool m_dataModified;
    QString m_currentSearchText;
    QTimer* m_searchTimer;

    // 列定义
    QList<ColumnDefinition> m_columnDefinitions;

    // 进度对话框
    AnimatedProgressDialog* m_progressDialog;

    // 数据缓存（用于大文件处理）
    bool m_largeFileMode;
    int m_maxDisplayRows;

    // 右键菜单相关
    QMenu* m_contextMenu;
    QAction* m_addRowAboveAction;
    QAction* m_addRowBelowAction;
    QAction* m_deleteRowsAction;
    QAction* m_addColumnLeftAction;
    QAction* m_addColumnRightAction;
    QAction* m_deleteColumnsAction;
    QPoint m_lastContextMenuPos;

    // 新增：压力导数计算器
    PressureDerivativeCalculator* m_pressureDerivativeCalculator;

    // 初始化方法
    void init();
    void setupModels();
    void setupConnections();
    void setupUI();
    void setupContextMenu();

    // 新增：设置压力导数计算器
    void setupPressureDerivativeCalculator();

    // 文件读取方法 - 优化后的方法
    bool loadExcelFile(const QString& filePath, QString& errorMessage);
    bool loadCsvFile(const QString& filePath, QString& errorMessage);
    bool loadJsonFile(const QString& filePath, QString& errorMessage);

    // 优化的Excel读取方法
    bool loadExcelFileOptimized(const QString& filePath, QString& errorMessage);
    bool quickDetectFileFormat(const QString& filePath);
    QString detectOptimalSeparator(const QString& filePath);

    // 新增：带配置的文件读取方法
    bool loadCsvFileWithConfig(const QString& filePath, const DataLoadConfigDialog::LoadConfig& config, QString& errorMessage);
    bool loadExcelFileWithConfig(const QString& filePath, const DataLoadConfigDialog::LoadConfig& config, QString& errorMessage);

#ifdef Q_OS_WIN
    bool loadExcelWithCOM(const QString& filePath, QString& errorMessage);
#endif

    bool loadExcelAsCSV(const QString& filePath, QString& errorMessage);
    bool loadCSVFile(const QString& filePath, const QString& separator, QString& errorMessage);
    QStringList splitCSVLine(const QString& line, const QString& separator);

    // 文件保存方法
    bool saveExcelFile(const QString& filePath);
    bool saveCsvFile(const QString& filePath);
    bool saveJsonFile(const QString& filePath);
    bool exportToPdf(const QString& filePath);
    bool exportToHtml(const QString& filePath);

    // 数据处理方法
    void removeEmptyRows();
    void removeEmptyColumns();
    void removeDuplicates();
    void fillMissingValues(const QString& method = "interpolation");
    void removeOutliers(double threshold = 2.0);
    void standardizeDataFormat();

    // 压降计算相关方法 - 优化的压降计算
    int findPressureColumn() const;
    int findTimeColumn() const;
    QString getPressureUnit() const;
    bool isValidPressureData(const QString& data) const;

    // 新增：压力导数计算相关方法
    // （自动检测列，无需配置对话框）

    // 修改后的时间转换相关方法
    QTime parseTimeString(const QString& timeStr) const;
    QDate parseDateString(const QString& dateStr) const;
    QDateTime combineDateAndTime(const QDate& date, const QTime& time) const;
    double calculateDateTimeDifference(const QDateTime& baseDateTime, const QDateTime& currentDateTime, const QString& unit) const;
    double calculateTimeDifference(const QTime& time1, const QTime& time2, const QString& unit) const;
    double convertTimeToUnit(double seconds, const QString& unit) const;
    bool isValidTimeFormat(const QString& timeStr) const;
    bool isValidDateFormat(const QString& dateStr) const;

    // 数据验证方法
    bool isNumericData(const QString& data) const;
    bool isDateTimeData(const QString& data) const;
    bool isValidRange(double value, double min, double max) const;
    QStringList detectDataType(int column) const;

    // UI更新方法
    void updateStatus(const QString& message, const QString& type = "info");
    void updateDataInfo();
    void setButtonsEnabled(bool enabled);
    void showAnimatedProgress(const QString& title, const QString& message);
    void hideAnimatedProgress();
    void updateProgress(int value, const QString& message = QString());

    // 数据处理方法
    void clearData();
    void applyColumnStyles();
    void optimizeColumnWidths();
    void optimizeTableDisplay();

    // 选择和交互方法
    int getSelectedRow() const;
    int getSelectedColumn() const;
    QList<int> getSelectedRows() const;
    QList<int> getSelectedColumns() const;
    bool checkDataModifiedAndPrompt();

    // 右键菜单辅助方法
    int getRowFromPosition(const QPoint& pos) const;
    int getColumnFromPosition(const QPoint& pos) const;

    // 辅助方法
    void emitDataChanged();
    QString formatNumber(double number, int precision = 3) const;
    void showStyledMessageBox(const QString& title, const QString& text,
                              QMessageBox::Icon icon, const QString& detailedText = "");

    // 列定义辅助方法
    ColumnDefinition getDefaultColumnDefinition(const QString& columnName);
    void applyColumnDefinition(int columnIndex, const ColumnDefinition& definition);
    bool validateColumnData(int columnIndex, const ColumnDefinition& definition, QStringList& errors) const;
};

#endif // DATAEDITORWIDGET_H
