#ifndef PLOTTINGWIDGET_H
#define PLOTTINGWIDGET_H

#include <QWidget>
#include <QTimer>
#include <QVector>
#include <QPointF>
#include <QPainter>
#include <QMouseEvent>
#include <QWheelEvent>
#include <QMenu>
#include <QAction>
#include <QFileDialog>
#include <QColorDialog>
#include <QInputDialog>
#include <QGroupBox>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QTabWidget>
#include <QPushButton>
#include <QComboBox>
#include <QCheckBox>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QTextEdit>
#include <QSplitter>
#include <QScrollArea>
#include <QSlider>
#include <QProgressBar>
#include <QStandardItemModel>
#include <QMessageBox>
#include <QLineEdit>
#include <QListWidget>
#include <QDialog>
#include <QTableWidget>
#include <QTableWidgetItem>
#include <QSvgGenerator>
#include <QMainWindow>
#include <QMdiArea>
#include <QMdiSubWindow>
#include <cmath>

namespace Ui {
class PlottingWidget;
}

// 线型枚举
enum class LineStyle {
    Solid,      // 实线
    Dash,       // 虚线
    Dot,        // 点线
    DashDot,    // 点划线
    DashDotDot  // 双点划线
};

// 坐标轴类型枚举
enum class AxisType {
    Linear,     // 常规坐标系
    Logarithmic // 对数坐标系
};

// 曲线类型枚举
enum class CurveType {
    Normal,     // 普通曲线
    Step        // 阶梯状曲线
};

// 添加线型转换函数声明
Qt::PenStyle lineStyleToQt(LineStyle style);
QString lineStyleToString(LineStyle style);
LineStyle stringToLineStyle(const QString &str);

struct WellTestData {
    QVector<double> time;           // 时间数据
    QVector<double> pressure;      // 压力数据
    QVector<double> pressureDerivative; // 压力导数数据
    QVector<double> deltaTime;     // 时间间隔
    QVector<double> deltaPressure; // 压力变化
    QString wellName;              // 井名
    QString testType;              // 测试类型
    QDateTime testDate;            // 测试日期
};

// 表格数据结构，用于存储导入的Excel/TXT数据
struct TableData {
    QStringList headers;                    // 列标题
    QVector<QVector<double>> columns;       // 列数据
    QString fileName;                       // 文件名
    int rowCount;                          // 行数
};

// 曲线数据结构，支持多曲线
struct CurveData {
    QString name;                          // 曲线名称
    QColor color;                          // 曲线颜色
    QVector<double> xData;                 // X轴数据
    QVector<double> yData;                 // Y轴数据
    bool visible;                          // 是否可见
    int lineWidth;                         // 线宽
    int pointSize;                         // 点大小
    QString xLabel;                        // X轴标签
    QString yLabel;                        // Y轴标签
    QString xUnit;                         // X轴单位
    QString yUnit;                         // Y轴单位
    QString curveType;                     // 曲线类型
    LineStyle lineStyle;                   // 线型
    AxisType xAxisType;                    // X轴类型
    AxisType yAxisType;                    // Y轴类型
    CurveType drawType;                    // 绘制类型（普通/阶梯）

    CurveData() : visible(true), lineWidth(2), pointSize(4), curveType("自定义"),
        lineStyle(LineStyle::Solid), xAxisType(AxisType::Linear), yAxisType(AxisType::Linear),
        drawType(CurveType::Normal) {}
};

struct PlotSettings {
    bool showGrid;
    bool logScaleX;
    bool logScaleY;
    QColor backgroundColor;
    QColor gridColor;
    QColor textColor;
    int lineWidth;
    int pointSize;
    QString xAxisTitle;
    QString yAxisTitle;
    QString plotTitle;
    double xMin, xMax, yMin, yMax;
    bool autoScale;
    bool showLegend;              // 是否显示图例
    QPointF legendPosition;       // 图例位置
    AxisType xAxisType;           // X轴类型
    AxisType yAxisType;           // Y轴类型
};

// 双图窗口类（用于压力产量联合显示）- 修改后版本
class DualPlotWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit DualPlotWindow(const QString &windowTitle, QWidget *parent = nullptr);
    ~DualPlotWindow();

    void addPressureCurve(const CurveData &curve);
    void addProductionCurve(const CurveData &curve);
    void updatePlots();
    void setAxisSettings(const QString &xTitle, const QString &pressureTitle, const QString &productionTitle);

protected:
    void resizeEvent(QResizeEvent *event) override;
    bool eventFilter(QObject *obj, QEvent *event) override;

private slots:
    // 鼠标事件处理槽函数
    void onPressureMousePress(QMouseEvent *event);
    void onPressureMouseMove(QMouseEvent *event);
    void onPressureMouseRelease(QMouseEvent *event);
    void onPressureWheel(QWheelEvent *event);

    void onProductionMousePress(QMouseEvent *event);
    void onProductionMouseMove(QMouseEvent *event);
    void onProductionMouseRelease(QMouseEvent *event);
    void onProductionWheel(QWheelEvent *event);

private:
    QWidget *m_centralWidget;
    QWidget *m_pressurePlotWidget;
    QWidget *m_productionPlotWidget;
    QSplitter *m_splitter;

    // 绘图数据
    QVector<CurveData> m_pressureCurves;
    QVector<CurveData> m_productionCurves;
    PlotSettings m_pressureSettings;
    PlotSettings m_productionSettings;
    QRect m_pressurePlotArea;
    QRect m_productionPlotArea;

    // 交互状态 - 分别为压力图和产量图
    bool m_pressureDragging;
    bool m_productionDragging;
    bool m_pressureSelecting;
    bool m_productionSelecting;
    bool m_pressurePanning;
    bool m_productionPanning;
    bool m_pressureLegendDragging;
    bool m_productionLegendDragging;

    QPoint m_lastPressureMousePos;
    QPoint m_lastProductionMousePos;
    QPoint m_pressureSelectionStart;
    QPoint m_productionSelectionStart;
    QRect m_pressureSelectionRect;
    QRect m_productionSelectionRect;
    QPoint m_pressureLegendOffset;
    QPoint m_productionLegendOffset;

    // 同步设置
    bool m_syncZoom;       // 是否同步缩放
    bool m_syncPan;        // 是否同步平移

    void setupUI();
    void initializePlotSettings();
    void paintPressurePlot();
    void paintProductionPlot();
    void synchronizeXAxis();
    void calculatePressureBounds();
    void calculateProductionBounds();

    // 绘图辅助函数
    void drawGridOnWidget(QPainter &painter, const QRect &plotArea, const PlotSettings &settings);
    void drawAxesOnWidget(QPainter &painter, const QRect &plotArea, const PlotSettings &settings);
    void drawLegendOnWidget(QPainter &painter, const QVector<CurveData> &curves,
                            const QRect &plotArea, const QPoint &offset);
    void drawCurveOnWidget(QPainter &painter, const CurveData &curve,
                           const QRect &plotArea, const PlotSettings &settings);

    // 坐标转换函数
    QPointF pressureDataToPixel(const QPointF &dataPoint);
    QPointF pixelToPressureData(const QPointF &pixelPoint);
    QPointF productionDataToPixel(const QPointF &dataPoint);
    QPointF pixelToProductionData(const QPointF &pixelPoint);

    // 缩放和平移函数
    void zoomPressureAtPoint(const QPointF &point, double factor);
    void zoomProductionAtPoint(const QPointF &point, double factor);
    void panPressureView(const QPointF &delta);
    void panProductionView(const QPointF &delta);
};

// 独立绘图窗口类
class PlotWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit PlotWindow(const QString &windowTitle, QWidget *parent = nullptr);
    ~PlotWindow();

    void addCurve(const CurveData &curve);
    void updatePlot();
    void setAxisSettings(bool logX, bool logY, const QString &xTitle, const QString &yTitle);
    void setPlotTitle(const QString &title);

protected:
    void paintEvent(QPaintEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void wheelEvent(QWheelEvent *event) override;
    void contextMenuEvent(QContextMenuEvent *event) override;
    bool eventFilter(QObject *obj, QEvent *event) override;

private slots:
    void onExportPlot();
    void onToggleGrid();
    void onToggleLegend();
    void onResetZoom();
    void onZoomIn();
    void onZoomOut();

private:
    QWidget *m_centralWidget;
    QWidget *m_plotWidget;

    // 绘图数据
    QVector<CurveData> m_curves;
    PlotSettings m_plotSettings;
    QRect m_plotArea;
    QRect m_legendArea;

    // 交互状态
    bool m_isDragging;
    bool m_isSelecting;
    bool m_isPanning;
    bool m_isDraggingLegend;
    QPoint m_lastMousePos;
    QPoint m_selectionStart;
    QRect m_selectionRect;
    QPoint m_legendOffset;

    // 标记和注释
    QVector<QPointF> m_markers;
    QVector<QPair<QPointF, QString>> m_annotations;

    // 右键菜单
    QMenu* m_contextMenu;
    QAction* m_exportAction;
    QAction* m_gridAction;
    QAction* m_legendAction;
    QAction* m_zoomInAction;
    QAction* m_zoomOutAction;
    QAction* m_resetZoomAction;

    void setupUI();
    void setupContextMenu();
    void initializePlotSettings();

    // 绘图函数
    void paintPlotArea(QPainter &painter);
    void drawBackground(QPainter &painter);
    void drawGrid(QPainter &painter);
    void drawAxes(QPainter &painter);
    void drawCurves(QPainter &painter);
    void drawCurve(QPainter &painter, const CurveData &curve);
    void drawStepCurve(QPainter &painter, const CurveData &curve);
    void drawLegend(QPainter &painter);
    void drawMarkers(QPainter &painter);
    void drawAnnotations(QPainter &painter);
    void drawSelection(QPainter &painter);

    // 坐标转换
    QPointF dataToPixel(const QPointF &dataPoint);
    QPointF pixelToData(const QPointF &pixelPoint);

    // 数据边界计算
    void calculateDataBounds();
    QPair<double, double> calculateOptimalRange(double min, double max, bool isLog);

    // 缩放和平移
    void zoomAtPoint(const QPointF &point, double factor);
    void panView(const QPointF &delta);

    // 辅助函数
    double transformToLogScale(double value, bool isLog);
    double transformFromLogScale(double value, bool isLog);
    QVector<double> generateOptimizedAxisLabels(double min, double max, AxisType axisType);
    QString formatAxisLabel(double value, bool isLog);
    QString formatScientific(double value, int decimals = 2);
    bool isValidDataPoint(double x, double y);
};

class PlottingWidget : public QWidget
{
    Q_OBJECT

public:
    explicit PlottingWidget(QWidget *parent = nullptr);
    ~PlottingWidget();

    // 数据管理
    void setWellTestData(const WellTestData &data);
    void addWellTestData(const WellTestData &data);
    void clearAllData();
    void removeDataSet(int index);

    // 设置表格数据
    void setTableData(const TableData &data);
    void setTableDataFromModel(QStandardItemModel* model, const QString &fileName = "");

    // 多曲线管理
    void addCurve(const CurveData &curve);
    void removeCurve(int index);
    void removeCurve(const QString &name);
    void updateCurve(int index, const CurveData &curve);
    void setCurveVisible(int index, bool visible);
    void setCurveVisible(const QString &name, bool visible);
    QVector<CurveData> getAllCurves() const;
    int getCurveCount() const;

    // 绘图控制
    void updatePlot();
    void resetZoom();
    void zoomIn();
    void zoomOut();
    void zoomToFit();

    // 新增单独缩放功能
    void zoomXIn();
    void zoomXOut();
    void zoomYIn();
    void zoomYOut();

    void exportPlot(const QString &fileName, const QString &format = "PNG");

    // 设置管理
    void setPlotSettings(const PlotSettings &settings);
    PlotSettings getPlotSettings() const;
    void resetToDefaultSettings();

    // 数据分析
    void performLogLogAnalysis();
    void performSemiLogAnalysis();
    void performCartesianAnalysis();
    void performDerivativeAnalysis();
    void performModelMatching();

signals:
    void dataPointClicked(double x, double y);
    void zoomChanged(double xMin, double xMax, double yMin, double yMax);
    void analysisCompleted(const QString &analysisType, const QMap<QString, double> &results);
    void plotExported(const QString &fileName);
    void curveAdded(const QString &curveName);
    void curveRemoved(const QString &curveName);

private slots:
    void onColorSettingsChanged();
    void onExportPlot();
    void onCurveSelectionChanged();
    void onMarkerAdded();
    void onAnnotationAdded();
    void onRemoveAllMarkers();
    void onRemoveLastMarker();
    void onRemoveAllAnnotations();
    void onAddCurve();
    void onEditCurve();
    void onRemoveCurve();

    // 修改后的快速绘图槽函数
    void onPressureProdDataPlot();  // 新的压力产量联合绘图
    void onPressureDerivativePlot();

protected:
    void paintEvent(QPaintEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void mouseDoubleClickEvent(QMouseEvent *event) override;
    void wheelEvent(QWheelEvent *event) override;
    void contextMenuEvent(QContextMenuEvent *event) override;
    bool eventFilter(QObject *obj, QEvent *event) override;

private:
    Ui::PlottingWidget *ui;
    void applyDialogStyle(QDialog *dialog);

    // 数据存储
    QVector<WellTestData> m_wellTestDataSets;
    WellTestData m_currentData;

    // 表格数据存储
    TableData m_tableData;
    bool m_hasTableData;

    // 多曲线数据存储
    QVector<CurveData> m_curves;

    // 绘图设置
    PlotSettings m_plotSettings;
    QRect m_plotArea;
    QRect m_legendArea;           // 图例区域

    // 交互状态
    bool m_isDragging;
    bool m_isSelecting;
    bool m_isPanning;
    bool m_isDraggingLegend;      // 图例拖拽状态
    QPoint m_lastMousePos;
    QPoint m_selectionStart;
    QRect m_selectionRect;
    QPoint m_legendDragStart;     // 图例拖拽起始点
    QPoint m_legendOffset;        // 图例偏移量

    // 缩放和平移 - 调整成员变量声明顺序
    double m_zoomFactor;
    double m_zoomFactorX; // X轴独立缩放
    double m_zoomFactorY; // Y轴独立缩放
    QPointF m_viewCenter;
    QPointF m_panOffset;

    // UI组件引用
    QCheckBox* m_showGridCheck;
    QCheckBox* m_showLegendCheck;
    QPushButton* m_gridColorBtn;

    // 曲线管理控件
    QListWidget* m_curvesListWidget;

    // 右键菜单
    QMenu* m_contextMenu;
    QMenu* m_zoomMenu;
    QMenu* m_dataMenu;
    QAction* m_addMarkerAction;
    QAction* m_addAnnotationAction;
    QAction* m_removeAllMarkersAction;
    QAction* m_removeLastMarkerAction;
    QAction* m_removeAllAnnotationsAction;
    QAction* m_zoomInAction;
    QAction* m_zoomOutAction;
    QAction* m_zoomFitAction;
    QAction* m_resetZoomAction;

    // 新增单独缩放菜单项
    QAction* m_zoomXInAction;
    QAction* m_zoomXOutAction;
    QAction* m_zoomYInAction;
    QAction* m_zoomYOutAction;

    // 标记和注释
    QVector<QPointF> m_markers;
    QVector<QPair<QPointF, QString>> m_annotations;

    // 坐标显示
    QLabel* m_coordinateLabel;

    // 绘图窗口管理
    QVector<PlotWindow*> m_plotWindows;
    QVector<DualPlotWindow*> m_dualPlotWindows;

    // 初始化函数
    void initializeUI();
    void setupContextMenu();
    void setupConnections();
    void setupDefaultSettings();

    // 绘图函数
    void paintPlotArea(QPaintEvent *event);
    void drawBackground(QPainter &painter);
    void drawGrid(QPainter &painter);
    void drawAxes(QPainter &painter);
    void drawAllCurves(QPainter &painter);
    void drawCurve(QPainter &painter, const CurveData &curve);
    void drawStepCurve(QPainter &painter, const CurveData &curve);
    void drawLegend(QPainter &painter);
    void drawMarkers(QPainter &painter);
    void drawAnnotations(QPainter &painter);
    void drawSelection(QPainter &painter);
    void drawNoDataMessage(QPainter &painter);
    void drawCoordinates(QPainter &painter);

    // 坐标转换函数
    QPointF dataToPixel(const QPointF &dataPoint);
    QPointF pixelToData(const QPointF &pixelPoint);

    // 对数坐标转换函数
    double transformToLogScale(double value, bool isLog);
    double transformFromLogScale(double value, bool isLog);

    // 坐标轴标签生成函数 - 优化后的版本，参考Saphir软件
    QVector<double> generateOptimizedAxisLabels(double min, double max, AxisType axisType);

    // 格式化坐标轴标签函数 - 新增，参考Saphir软件风格
    QString formatAxisLabel(double value, bool isLog);

    // 计算函数
    void calculateDataBounds();

    // 曲线管理函数
    void updateCurvesList();

    // 鼠标缩放函数
    void zoomAtPoint(const QPointF &point, double factor);
    void panView(const QPointF &delta);

    // 辅助函数
    void updateControlsFromSettings();
    QString formatScientific(double value, int decimals = 2);

    // 数据处理函数
    bool isValidDataPoint(double x, double y);

    // 数据范围计算优化函数
    QPair<double, double> calculateOptimalRange(double min, double max, bool isLog);

    // 绘图窗口创建函数
    PlotWindow* createPlotWindow(const QString &title, const QString &dataType);
    DualPlotWindow* createDualPlotWindow(const QString &title);

    // 对话框函数
    void showDataSelectionDialog(const QString &plotType);
    void showPressureProdDataDialog();  // 新的压力产量联合对话框

    // 数据创建函数
    CurveData createCurveFromTableData(int xColumn, int yColumn, const QString &curveName,
                                       const QColor &color, AxisType xAxisType, AxisType yAxisType,
                                       const QString &xLabel = "", const QString &yLabel = "",
                                       const QString &xUnit = "", const QString &yUnit = "",
                                       int lineWidth = 2, int pointSize = 4);

    CurveData createProductionCurve(int timeIndex, int productionIndex, const QString &curveName,
                                    const QColor &color, AxisType timeAxisType, AxisType productionAxisType,
                                    const QString &timeLabel = "", const QString &productionLabel = "",
                                    const QString &timeUnit = "", const QString &productionUnit = "",
                                    const QString &curveTypeStr = "时间vs产量", int lineWidth = 2, int pointSize = 4);

    // 阶梯状曲线数据创建函数
    CurveData createStepProductionCurve(int timeIndex, int productionIndex, const QString &curveName,
                                        const QColor &color, AxisType timeAxisType, AxisType productionAxisType,
                                        const QString &timeLabel = "", const QString &productionLabel = "",
                                        const QString &timeUnit = "", const QString &productionUnit = "",
                                        int lineWidth = 2, int pointSize = 4);
};

// 美化的图例绘制函数（可以作为全局辅助函数）
void drawBeautifulLegend(QPainter &painter,
                         const QVector<CurveData> &curves,
                         const QRect &plotArea,
                         const QPoint &offset);

#endif // PLOTTINGWIDGET_H
