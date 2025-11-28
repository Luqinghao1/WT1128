#ifndef MODELWIDGET2_H
#define MODELWIDGET2_H

#include <QWidget>
#include <QTimer>
#include <QFuture>
#include <QFutureWatcher>
#include <QtConcurrent/QtConcurrent>
#include <QPainter>
#include <QPaintEvent>
#include <QMouseEvent>
#include <QWheelEvent>

namespace Ui {
class ModelWidget2;
}

// 计算参数结构体 - 有限导流模型
struct CalculationParametersFD {
    double omega = 0.0155211;     // 储容比
    double S = 0.810876;          // 表皮因子
    double cD = 8.08669e-8;       // 井筒存储系数
    double lambda = 0.083277;     // 窜流系数
    int mf = 3;                   // 裂缝条数
    int nf = 5;                   // 离散段数
    double Xf = 193.0;            // 裂缝长度
    double yy = 295.0;            // 裂缝间距控制
    double y = 2758.0;            // 水平井长度
    int N = 4;                    // Stefest算法参数
    double CFD = 0.9;             // 有限导流系数
    double kpd = 0.04;            // 压力敏感因子

    // 曲线拟合参数
    bool enableSmoothing = true;     // 是否启用平滑
    int fittingPoints = 100;          // 拟合后的点数
    int originalPoints = 100;         // 原始计算点数
};

// 计算结果结构体 - 有限导流模型
struct CalculationResultFD {
    QVector<double> tD;       // 时间
    QVector<double> pd;       // 压力（原始）
    QVector<double> pd1;      // 压力（考虑压力敏感）
    QVector<double> dpd;      // 压力导数
    QVector<double> td_dpd;   // 压力导数对应的时间点

    // 平滑后的数据
    QVector<double> tD_smooth;
    QVector<double> pd_smooth;
    QVector<double> dpd_smooth;
    QVector<double> td_dpd_smooth;

    bool isValid = false;
    bool hasSmoothedData = false;
    QString errorMessage;
};

// 简化的样条插值类（与modelwidget1相同）
class SimpleSplineFD {
public:
    void setData(const QVector<double>& x, const QVector<double>& y);
    double interpolate(double x) const;
    QVector<double> interpolateRange(const QVector<double>& xNew) const;
    bool isValid() const { return m_isValid; }

private:
    void computeCoefficients();
    int findSegment(double x) const;

    QVector<double> m_x, m_y;
    QVector<double> m_a, m_b, m_c, m_d;
    bool m_isValid = false;
};

// 自定义图表绘制组件（与modelwidget1相同但独立）
class LogLogChartWidgetFD : public QWidget
{
    Q_OBJECT

public:
    explicit LogLogChartWidgetFD(QWidget *parent = nullptr);

    void setData(const QVector<double>& xData, const QVector<double>& yData1,
                 const QVector<double>& yData2, const QVector<double>& xData2 = QVector<double>());
    void clearData();
    void resetView();
    void autoFitData();
    void setShowOriginalData(bool show) { m_showOriginalData = show; update(); }
    void setShowSmoothedCurve(bool show) { m_showSmoothedCurve = show; update(); }

protected:
    void paintEvent(QPaintEvent *event) override;
    void mousePressEvent(QMouseEvent *event) override;
    void mouseMoveEvent(QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void wheelEvent(QWheelEvent *event) override;

private:
    void drawAxis(QPainter& painter, const QRect& plotRect);
    void drawData(QPainter& painter, const QRect& plotRect);
    void drawLegend(QPainter& painter, const QRect& plotRect);
    QPointF dataToPixel(double x, double y, const QRect& plotRect);

    // 原始数据
    QVector<double> m_xData;
    QVector<double> m_yData1;
    QVector<double> m_yData2;
    QVector<double> m_xData2;

    QString m_title;

    // 显示范围
    double m_xMin, m_xMax, m_yMin, m_yMax;
    bool m_hasData;

    // 显示选项
    bool m_showOriginalData;
    bool m_showSmoothedCurve;

    // 鼠标交互
    bool m_isDragging;
    QPoint m_lastMousePos;
};

class ModelWidget2 : public QWidget
{
    Q_OBJECT

public:
    explicit ModelWidget2(QWidget *parent = nullptr);
    ~ModelWidget2();

private slots:
    void onCalculateClicked();
    void onParameterChanged();
    void onExportResults();
    void onResetParameters();
    void onCalculationFinished();
    void onResetView();
    void onFitToData();
    void onShowOriginalDataChanged(bool show);
    void onShowSmoothedCurveChanged(bool show);
    void onSmoothingEnabledChanged(bool enabled);

signals:
    void calculationCompleted(const QString &analysisType, const QMap<QString, double> &results);

private:
    void resetParametersToDefault();
    void updateChart(const CalculationResultFD& result);
    void updateResultText(const CalculationResultFD& result);
    void setCalculationInProgress(bool inProgress);

    // 计算函数
    CalculationResultFD performCalculation(const CalculationParametersFD& params);
    void smoothData(CalculationResultFD& result, const CalculationParametersFD& params);

    // 数学计算函数
    QVector<double> solveLinearSystem(const QVector<QVector<double>>& A, const QVector<double>& b);
    double stefestCoefficient(int i, int N);
    double factorial(int n);
    double flaplace(double z, const CalculationParametersFD& params);
    double e_function(double z, int i, int j, int k, int v, int mf, int nf,
                      double omega, double lambda, double Xf, double yy, double y,
                      double CFD, double deltaL);
    double f_function(int j, int nf, double Xf, double y);

    // 数值积分和贝塞尔函数
    double integralBesselK0(double XDkv, double YDkv, double yDij,
                            double fz, double xDij, double xDij1);
    double besselK0(double x);
    double gaussQuadrature(double XDkv, double YDkv, double yDij, double fz,
                           double a, double b);

    // 压力导数计算（改进版本）
    void computePressureDerivative(const QVector<double>& tD, const QVector<double>& pd,
                                   double cD, QVector<double>& dpd, QVector<double>& td_dpd);

    // 数据平滑
    QVector<double> movingAverage(const QVector<double>& data, int windowSize);

private:
    Ui::ModelWidget2 *ui;
    LogLogChartWidgetFD* m_chartWidget;
    CalculationParametersFD m_currentParams;
    QFutureWatcher<CalculationResultFD>* m_calculationWatcher;
    QTimer* m_progressTimer;
    bool m_isCalculating;
    CalculationResultFD m_lastResult;
};

#endif // MODELWIDGET2_H
