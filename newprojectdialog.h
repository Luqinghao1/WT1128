#ifndef NEWPROJECTDIALOG_H
#define NEWPROJECTDIALOG_H

#include <QDialog>
#include <QDateTime>
#include <QListWidgetItem>
#include <QStackedWidget>
#include <QMap>
#include <QDir>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QMessageBox>
#include <QCloseEvent>
#include <QPropertyAnimation>
#include <QGraphicsEffect>
#include <QTimer>
#include <QFileDialog>
#include <QStandardPaths>
#include <QJsonDocument>
#include <QJsonObject>

namespace Ui {
class NewProjectDialog;
}

class NewProjectDialog : public QDialog
{
    Q_OBJECT

public:
    explicit NewProjectDialog(QWidget *parent = nullptr);
    ~NewProjectDialog();

    // 项目信息结构体，用于存储表单数据
    struct ProjectInfo {
        // 基本信息
        QString projectName;
        QString wellName;
        QString fieldName;
        QString engineerName;
        QDateTime creationDate;
        QString projectPath;        // 新增：项目文件路径

        // 井型和储层类型
        QString wellType;           // 水平井、垂直井、定向井
        QString reservoirType;      // 页岩油、致密油、常规油藏等
        QString fractureType;       // 多段压裂、单段压裂、未压裂

        // 测试信息
        QString testType;
        QString analysisType;
        QDate testDate;
        double testDuration;
        QString testDurationUnit;

        // 井筒参数 - 统一使用国际标准单位制
        double wellRadius;          // m
        double horizontalLength;    // m
        double perforationTopDepth; // m
        double perforationBottomDepth; // m
        double payZone;             // m
        double skinFactor;

        // 压裂参数
        int fractureStages;         // 压裂段数
        double fractureSpacing;     // m
        double fractureHalfLength;  // m
        double fractureConductivity; // mD·m

        // 井下设备
        double tuningSize;          // mm
        double casingSize;          // mm
        QString completionType;
        double bhpGaugeDepth;       // m

        // 储层参数
        QString reservoirModel;
        double initialPressure;     // MPa
        double porosity;            // 小数
        double matrixPorosity;
        double fracturePorosity;
        double reservoirTemp;       // °C
        double permeability;        // mD
        double matrixPermeability;
        double fracturePermeability;
        double rockCompressibility; // 1/MPa

        // 边界条件
        QString boundaryType;
        double distanceToBoundary;  // m
        double anisotropyRatio;
        double formationDip;        // 度

        // 流体属性
        QString fluidSystem;
        double referencePressure;   // MPa
        bool hasOil;
        bool hasGas;
        bool hasWater;
        double viscosity;           // mPa·s
        double fvf;                 // m³/m³
        double fluidCompressibility; // 1/MPa
        double totalCompressibility; // 1/MPa

        // 分析设置
        QString analysisModel;
        QString flowRegime;
        QString solutionMethod;
        QString optimizationMethod;
        QString pressureDifferenceDefinition;
        QString timeFormat;
        QString primaryPlotType;
        bool showGrid;
        bool autoAdjustScale;
        bool generateReport;
        bool saveParameters;

        // 高级设置
        bool enableDualPorosity;
        bool enableTriplePorosity;
        bool enableTransientFlow;
        bool enableNonDarcyFlow;
        QString interferenceModel;
    };

    // 获取项目信息
    ProjectInfo getProjectInfo() const;

signals:
    // 项目创建信号
    void projectCreated(const ProjectInfo &projectInfo);

protected:
    // 重写关闭事件
    void closeEvent(QCloseEvent *event) override;
    void showEvent(QShowEvent *event) override;

private slots:
    // 导航列表项被选中
    void onNavigationItemClicked(QListWidgetItem *item);

    // 预设按钮点击
    void applyShaleOilPreset();       // 页岩油预设
    void applyTightOilPreset();       // 致密油预设
    void applyConventionalOilPreset(); // 常规油藏预设


    // 流体相选择变化
    void onFluidPhaseChanged();

    // 分析类型变化
    void onAnalysisTypeChanged();

    // 井型变化
    void onWellTypeChanged();

    // 储层类型变化
    void onReservoirTypeChanged();

    // 压裂类型变化
    void onFractureTypeChanged();

    // 总压缩系数计算
    void calculateTotalCompressibility();

    // 压裂参数计算
    void calculateFractureParameters();

    // 按钮响应函数
    void onNextButtonClicked();
    void onBackButtonClicked();
    void onCancelButtonClicked();
    void onHelpButtonClicked();

    // 验证当前页面数据
    bool validateCurrentPage();
    bool validateCurrentPageSilently();  // 新增：静默验证函数

    // 数据变化时自动保存
    void onDataChanged();

    // 动画完成
    void onAnimationFinished();

    // 项目路径选择
    void onSelectProjectPath();

private:
    void init();
    void setupConnections();
    void populateDefaultValues();
    void setupAnimations();

    // 更新状态信息
    void updateStatusLabel();
    void updateProgressBar();

    // 收集所有数据
    void collectAllData();

    // 表单验证 - 带错误提示的版本
    bool validateBasicInfo();
    bool validateWellParams();
    bool validateReservoirParams();
    bool validateFluidProperties();
    bool validateAnalysisSettings();
    bool validateFractureParams();

    // 表单验证 - 静默版本（不显示错误对话框）
    bool validateBasicInfoSilently();
    bool validateWellParamsSilently();
    bool validateReservoirParamsSilently();
    bool validateFluidPropertiesSilently();
    bool validateAnalysisSettingsSilently();
    bool validateFractureParamsSilently();

    // 设置控件可见性
    void updateControlsVisibility();

    // 界面动效
    void animatePageTransition(int fromPage, int toPage);

    // 项目文件操作
    bool createProjectFile();
    QString generateProjectFileName();
    bool saveProjectToFile(const QString &filePath);

private:
    Ui::NewProjectDialog *ui;

    // 导航页面映射
    QMap<QString, int> m_pageMap;

    // 当前页面索引
    int m_currentPageIndex;

    // 项目信息
    ProjectInfo m_projectInfo;

    // 表单是否有修改
    bool m_isModified;

    // 是否正在应用预设
    bool m_isApplyingPreset;

    // 是否正在初始化
    bool m_isInitializing;

    // 动画相关
    QPropertyAnimation *m_fadeAnimation;
    QTimer *m_validationTimer;

    // 页面验证状态
    QMap<int, bool> m_pageValidationStatus;
};

#endif // NEWPROJECTDIALOG_H
