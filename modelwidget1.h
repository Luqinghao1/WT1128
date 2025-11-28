#ifndef MODELWIDGET1_H
#define MODELWIDGET1_H

#include <QWidget>
#include <QPainter>
#include <QVector>
#include <QMap>
#include <QMouseEvent>
#include <QWheelEvent>
#include <functional>

namespace Ui {
class ModelWidget1;
}

// =========================================================
// LogLogChartWidget: 绘图组件
// =========================================================
class LogLogChartWidget : public QWidget
{
    Q_OBJECT

public:
    explicit LogLogChartWidget(QWidget *parent = nullptr);

    void setData(const QVector<double>& xData, const QVector<double>& yData1,
                 const QVector<double>& yData2, const QVector<double>& xData2 = QVector<double>());
    void clearData();
    void resetView();
    void autoFitData();
    void setShowOriginalData(bool show) { m_showOriginalData = show; update(); }

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

    QVector<double> m_xData;
    QVector<double> m_yData1;
    QVector<double> m_yData2;
    QVector<double> m_xData2;

    QString m_title;
    double m_xMin, m_xMax, m_yMin, m_yMax;
    bool m_hasData;
    bool m_showOriginalData;
    bool m_isDragging;
    QPoint m_lastMousePos;
};

// =========================================================
// ModelWidget1: 数学核心
// =========================================================
class ModelWidget1 : public QWidget
{
    Q_OBJECT

public:
    explicit ModelWidget1(QWidget *parent = nullptr);
    ~ModelWidget1();

private slots:
    void onCalculateClicked();
    void onResetParameters();
    void onExportResults();
    void onResetView();
    void onFitToData();

signals:
    void calculationCompleted(const QString &analysisType, const QMap<QString, double> &results);

private:
    void runCalculation();

    double flaplace(double s, double cD, double S, int mf, int nf,
                    double omega, double lambda, double Xf, double yy, double y, int N);

    double getStefestVi(int i, int N);
    double factorial(int n);
    double besselK0(double x);
    QVector<double> gaussElimination(QVector<QVector<double>> A, QVector<double> b);

    double integral_bessel(double XDkv, double YDkv, double yDij, double fz, double a, double b);

    double adaptiveGauss(std::function<double(double)> f, double a, double b, double eps, int depth, int maxDepth);
    double gauss15(std::function<double(double)> f, double a, double b);

private:
    Ui::ModelWidget1 *ui;
    LogLogChartWidget *m_chartWidget;

    QVector<double> res_tD;
    QVector<double> res_pD;
    QVector<double> res_dpD;
    QVector<double> res_td_dpD;
};

#endif // MODELWIDGET1_H
