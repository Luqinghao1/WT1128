#ifndef MODELWIDGET3_H
#define MODELWIDGET3_H

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
class ModelWidget3;
}

// 计算参数结构体 - 分段多簇压裂水平井双重孔隙介质页岩油藏渗流模型
struct CalculationParametersMix {
    double omega1 = 0.05;         // 储容比1
    double omega2 = 0.05;         // 储容比2
    double S = 0.1;               // 表皮因子
    double cD = 1e-8;             // 井筒存储系数
    double lambda1 = 1e-1;        // 窜流系数1
    double lambda2 = 1e-1;        // 窜流系数2
    int mf1 = 2;                  // 裂缝条数1
    int mf2 = 2;                  // 裂缝条数2
    int nf = 5;                   // 离散段数
    double Xf1 = 40.0;            // 裂缝长度1
    double Xf2 = 40.0;            // 裂缝长度2
    double yy1 = 70.0;            // 裂缝间距控制1
    double yy2 = 70.0;            // 裂缝间距控制2
    double y = 800.0;             // 水平井长度
    int N = 4;                    // Stefest算法参数
    double CFD1 = 0.4;            // 有限导流系数1
    double CFD2 = 0.4;            // 有限导流系数2
    double kpd = 0.045;           // 压力敏感因子

    // 曲线拟合参数
    bool enableSmoothing = true;     // 是否启用平滑
    int fittingPoints = 100;
    int originalPoints = 100;
};

// 计算结果结构体 - 分段多簇压裂水平井双重孔隙介质页岩油藏渗流模型
struct CalculationResultMix {
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

// 简化的样条插值类
class SimpleSplineMix {
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

// 自定义图表绘制组件
class LogLogChartWidgetMix : public QWidget
{
    Q_OBJECT

public:
    explicit LogLogChartWidgetMix(QWidget *parent = nullptr);

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

class ModelWidget3 : public QWidget
{
    Q_OBJECT

public:
    explicit ModelWidget3(QWidget *parent = nullptr);
    ~ModelWidget3();

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
    void updateChart(const CalculationResultMix& result);
    void updateResultText(const CalculationResultMix& result);
    void setCalculationInProgress(bool inProgress);

    // 计算函数
    CalculationResultMix performCalculation(const CalculationParametersMix& params);
    void smoothData(CalculationResultMix& result, const CalculationParametersMix& params);

    // 数学计算函数
    QVector<double> solveLinearSystem(const QVector<QVector<double>>& A, const QVector<double>& b);
    double stefestCoefficient(int i, int N);
    double factorial(int n);
    double flaplace(double z, const CalculationParametersMix& params);

    // 混合模型的专用函数 - 修正函数签名
    double e1_function(double z, int i, int j, int k, int v, int mf1, int nf,
                       double omega1, double lambda1, double Xf1, double yy1, double y, double CFD1);
    double e2_function(double z, int i, int j, int k, int v, int mf1, int mf2, int nf,
                       double omega2, double lambda2, double Xf1, double Xf2, double yy1, double yy2, double y, double CFD2);
    double e3_function(double z, int i, int j, int k, int v, int mf1, int mf2, int nf,
                       double omega1, double lambda1, double Xf1, double Xf2, double yy1, double yy2, double y, double CFD1);
    double e4_function(double z, int i, int j, int k, int v, int mf2, int nf,
                       double omega2, double lambda2, double Xf1, double Xf2, double yy2, double y, double CFD2);
    // 修正f1和f2函数签名 - 添加了i参数以与MATLAB一致
    double f1_function(int i, int v, int nf, double Xf1, double y);
    double f2_function(int i, int v, int nf, double Xf1, double Xf2, double y);

    // 数值积分和贝塞尔函数
    double integralBesselK0(double XDkv, double YDkv, double yDij,
                            double fz, double xDij, double xDij1);
    double besselK0(double x);
    double gaussQuadrature(double XDkv, double YDkv, double yDij, double fz,
                           double a, double b);

    // 压力导数计算
    void computePressureDerivative(const QVector<double>& tD, const QVector<double>& pd,
                                   double cD, QVector<double>& dpd, QVector<double>& td_dpd);

    // 数据平滑
    QVector<double> movingAverage(const QVector<double>& data, int windowSize);

private:
    Ui::ModelWidget3 *ui;
    LogLogChartWidgetMix* m_chartWidget;
    CalculationParametersMix m_currentParams;
    QFutureWatcher<CalculationResultMix>* m_calculationWatcher;
    QTimer* m_progressTimer;
    bool m_isCalculating;
    CalculationResultMix m_lastResult;
};

#endif // MODELWIDGET3_H
