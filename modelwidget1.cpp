#include "modelwidget1.h"
#include "ui_modelwidget1.h"
#include <cmath>
#include <algorithm>
#include <QDebug>
#include <QMessageBox>
#include <QFileDialog>
#include <QTextStream>
#include <QDateTime>
#include <QCoreApplication> // 用于 processEvents

// =========================================================
// LogLogChartWidget 实现 (保持不变)
// =========================================================

LogLogChartWidget::LogLogChartWidget(QWidget *parent)
    : QWidget(parent)
    , m_xMin(1e1), m_xMax(1e15)
    , m_yMin(1e-2), m_yMax(1e1)
    , m_hasData(false)
    , m_showOriginalData(true)
    , m_isDragging(false)
{
    setMinimumSize(600, 400);
    setStyleSheet("QWidget { background-color: white; border: 1px solid gray; }");
    setMouseTracking(true);
}

void LogLogChartWidget::setData(const QVector<double>& xData, const QVector<double>& yData1,
                                const QVector<double>& yData2, const QVector<double>& xData2)
{
    m_xData = xData;
    m_yData1 = yData1;
    m_yData2 = yData2;
    m_xData2 = xData2.isEmpty() ? xData : xData2;
    m_hasData = !xData.isEmpty() && !yData1.isEmpty();
    if (m_hasData) autoFitData();
    update();
}

void LogLogChartWidget::clearData() {
    m_hasData = false; m_xData.clear(); m_yData1.clear(); m_yData2.clear(); m_xData2.clear();
    resetView(); update();
}

void LogLogChartWidget::resetView() {
    if (m_hasData) autoFitData();
    else { m_xMin = 1e1; m_xMax = 1e15; m_yMin = 1e-2; m_yMax = 1e1; }
    update();
}

void LogLogChartWidget::autoFitData() {
    if (!m_hasData) return;
    QVector<double> allX, allY;
    for (double x : m_xData) if (x > 0 && std::isfinite(x)) allX.append(x);
    for (double x : m_xData2) if (x > 0 && std::isfinite(x)) allX.append(x);
    for (double y : m_yData1) if (y > 0 && std::isfinite(y)) allY.append(y);
    for (double y : m_yData2) if (y > 0 && std::isfinite(y)) allY.append(y);
    if (allX.isEmpty() || allY.isEmpty()) return;

    auto [xMinIt, xMaxIt] = std::minmax_element(allX.begin(), allX.end());
    auto [yMinIt, yMaxIt] = std::minmax_element(allY.begin(), allY.end());
    double logXMin = log10(*xMinIt), logXMax = log10(*xMaxIt);
    double logYMin = log10(*yMinIt), logYMax = log10(*yMaxIt);
    double xMargin = (logXMax - logXMin) * 0.05;
    double yMargin = (logYMax - logYMin) * 0.05;

    m_xMin = pow(10.0, logXMin - xMargin); m_xMax = pow(10.0, logXMax + xMargin);
    m_yMin = pow(10.0, logYMin - yMargin); m_yMax = pow(10.0, logYMax + yMargin);
}

void LogLogChartWidget::mousePressEvent(QMouseEvent *event) {
    if (event->button() == Qt::LeftButton) {
        m_isDragging = true; m_lastMousePos = event->pos();
        setCursor(Qt::ClosedHandCursor);
    }
}

void LogLogChartWidget::mouseMoveEvent(QMouseEvent *event) {
    if (m_isDragging) {
        QPoint delta = event->pos() - m_lastMousePos;
        m_lastMousePos = event->pos();
        QRect plotRect = rect().adjusted(80, 50, -50, -80);
        double logXRange = log10(m_xMax) - log10(m_xMin);
        double logYRange = log10(m_yMax) - log10(m_yMin);
        double deltaLogX = -delta.x() * logXRange / plotRect.width();
        double deltaLogY = delta.y() * logYRange / plotRect.height();
        m_xMin *= pow(10.0, deltaLogX); m_xMax *= pow(10.0, deltaLogX);
        m_yMin *= pow(10.0, deltaLogY); m_yMax *= pow(10.0, deltaLogY);
        update();
    } else setCursor(Qt::OpenHandCursor);
}

void LogLogChartWidget::mouseReleaseEvent(QMouseEvent *event) {
    if (event->button() == Qt::LeftButton) {
        m_isDragging = false; setCursor(Qt::OpenHandCursor);
    }
}

void LogLogChartWidget::wheelEvent(QWheelEvent *event) {
    double factor = (event->angleDelta().y() > 0) ? 1.15 : 1.0/1.15;
    m_xMin *= (1.0/factor); m_xMax *= (1.0/factor);
    m_yMin *= (1.0/factor); m_yMax *= (1.0/factor);
    update();
}

void LogLogChartWidget::paintEvent(QPaintEvent *) {
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);
    QRect plotRect = rect().adjusted(80, 50, -50, -80);
    painter.fillRect(rect(), Qt::white);
    painter.fillRect(plotRect, QColor(250, 250, 250));

    drawAxis(painter, plotRect);
    if (m_hasData) { drawData(painter, plotRect); drawLegend(painter, plotRect); }
    else {
        painter.setPen(Qt::gray);
        painter.drawText(plotRect, Qt::AlignCenter, "点击'开始计算'生成压力分析图表");
    }
}

void LogLogChartWidget::drawAxis(QPainter& painter, const QRect& plotRect) {
    painter.setPen(QPen(Qt::black, 2)); painter.drawRect(plotRect);
    painter.setPen(QPen(Qt::lightGray, 1)); painter.setFont(QFont("Arial", 9));

    double logXMin = log10(m_xMin), logXMax = log10(m_xMax);
    for (int exp = floor(logXMin); exp <= ceil(logXMax); exp++) {
        double x = pow(10.0, exp);
        if (x >= m_xMin && x <= m_xMax) {
            QPointF p = dataToPixel(x, m_yMin, plotRect);
            if (p.x() >= plotRect.left() && p.x() <= plotRect.right()) {
                painter.drawLine(p.x(), plotRect.bottom(), p.x(), plotRect.top());
                painter.setPen(Qt::black);
                painter.drawText(QRect(p.x()-25, plotRect.bottom()+5, 50, 20), Qt::AlignCenter, QString("1e%1").arg(exp));
                painter.setPen(Qt::lightGray);
            }
        }
    }
    double logYMin = log10(m_yMin), logYMax = log10(m_yMax);
    for (int exp = floor(logYMin); exp <= ceil(logYMax); exp++) {
        double y = pow(10.0, exp);
        if (y >= m_yMin && y <= m_yMax) {
            QPointF p = dataToPixel(m_xMin, y, plotRect);
            if (p.y() >= plotRect.top() && p.y() <= plotRect.bottom()) {
                painter.drawLine(plotRect.left(), p.y(), plotRect.right(), p.y());
                painter.setPen(Qt::black);
                painter.drawText(QRect(plotRect.left()-50, p.y()-10, 45, 20), Qt::AlignRight|Qt::AlignVCenter, QString("1e%1").arg(exp));
                painter.setPen(Qt::lightGray);
            }
        }
    }
    painter.setPen(Qt::black); painter.setFont(QFont("Arial", 11, QFont::Bold));
    painter.drawText(plotRect.center().x()-30, plotRect.bottom()+40, "tD/CD");
    painter.save();
    painter.translate(plotRect.left()-60, plotRect.center().y());
    painter.rotate(-90);
    painter.drawText(-50, 0, "PD & dPD");
    painter.restore();
}

void LogLogChartWidget::drawData(QPainter& painter, const QRect& plotRect) {
    auto drawCurve = [&](const QVector<double>& x, const QVector<double>& y, QColor color) {
        if (x.isEmpty() || y.isEmpty()) return;
        painter.setPen(QPen(color, 2));
        QVector<QPointF> points;
        for (int i=0; i<qMin(x.size(), y.size()); ++i) {
            if (x[i]>0 && y[i]>0) {
                QPointF p = dataToPixel(x[i], y[i], plotRect);
                points.append(p);
            }
        }
        for (int i=1; i<points.size(); ++i) {
            if (plotRect.contains(points[i-1].toPoint()) || plotRect.contains(points[i].toPoint()))
                painter.drawLine(points[i-1], points[i]);
        }
        if (m_showOriginalData) {
            painter.setPen(QPen(color, 1)); painter.setBrush(color);
            for (auto& p : points) if (plotRect.contains(p.toPoint())) painter.drawEllipse(p, 2, 2);
        }
    };
    drawCurve(m_xData, m_yData1, Qt::red);
    drawCurve(m_xData2, m_yData2, Qt::blue);
}

void LogLogChartWidget::drawLegend(QPainter& painter, const QRect& plotRect) {
    painter.setFont(QFont("Arial", 10));
    int x = plotRect.right()-100, y = plotRect.top()+20;
    painter.setPen(QPen(Qt::red, 2)); painter.drawLine(x, y, x+20, y);
    painter.setPen(Qt::black); painter.drawText(x+25, y+5, "PD");
    if (!m_yData2.isEmpty()) {
        y += 20;
        painter.setPen(QPen(Qt::blue, 2)); painter.drawLine(x, y, x+20, y);
        painter.setPen(Qt::black); painter.drawText(x+25, y+5, "dPD");
    }
}

QPointF LogLogChartWidget::dataToPixel(double x, double y, const QRect& plotRect) {
    x = qMax(x, 1e-20); y = qMax(y, 1e-20);
    double lx = log10(x), ly = log10(y);
    double lxmin = log10(qMax(m_xMin, 1e-20)), lxmax = log10(qMax(m_xMax, 1e-20));
    double lymin = log10(qMax(m_yMin, 1e-20)), lymax = log10(qMax(m_yMax, 1e-20));
    double px = plotRect.left() + (lx - lxmin)/(lxmax - lxmin)*plotRect.width();
    double py = plotRect.bottom() - (ly - lymin)/(lymax - lymin)*plotRect.height();
    return QPointF(px, py);
}

// =========================================================
// ModelWidget1 主逻辑
// =========================================================

ModelWidget1::ModelWidget1(QWidget *parent) : QWidget(parent), ui(new Ui::ModelWidget1) {
    ui->setupUi(this);
    m_chartWidget = new LogLogChartWidget(this);
    QVBoxLayout *layout = new QVBoxLayout(ui->chartTab);
    layout->addWidget(m_chartWidget);
    layout->setContentsMargins(0,0,0,0);

    connect(ui->calculateButton, &QPushButton::clicked, this, &ModelWidget1::onCalculateClicked);
    connect(ui->resetButton, &QPushButton::clicked, this, &ModelWidget1::onResetParameters);
    connect(ui->exportButton, &QPushButton::clicked, this, &ModelWidget1::onExportResults);
    connect(ui->resetViewButton, &QPushButton::clicked, this, &ModelWidget1::onResetView);
    connect(ui->fitToDataButton, &QPushButton::clicked, this, &ModelWidget1::onFitToData);

    onResetParameters();
}

ModelWidget1::~ModelWidget1() { delete ui; }

void ModelWidget1::onResetParameters() {
    ui->omegaSpinBox->setValue(0.05);
    ui->sSpinBox->setValue(1.0);
    ui->cDSpinBox->setValue(1e-8);
    ui->lambdaSpinBox->setValue(0.1);
    ui->mfSpinBox->setValue(3);
    ui->nfSpinBox->setValue(5);
    ui->xfSpinBox->setValue(40);
    ui->yySpinBox->setValue(70);
    ui->ySpinBox->setValue(1000);
    ui->nSpinBox->setValue(4);
}

void ModelWidget1::onResetView() { m_chartWidget->resetView(); }
void ModelWidget1::onFitToData() { m_chartWidget->autoFitData(); }

void ModelWidget1::onCalculateClicked() {
    ui->calculateButton->setEnabled(false);
    ui->calculateButton->setText("计算中...");

    // 强制刷新一次界面
    QCoreApplication::processEvents();

    runCalculation();

    ui->calculateButton->setEnabled(true);
    ui->calculateButton->setText("开始计算");
    ui->exportButton->setEnabled(true);
    ui->resetViewButton->setEnabled(true);
    ui->fitToDataButton->setEnabled(true);
    ui->tabWidget->setCurrentIndex(0);
}

// ------------------------------------------------------------------------
// 数学计算核心 (性能优化版)
// ------------------------------------------------------------------------

void ModelWidget1::runCalculation() {
    double omega = ui->omegaSpinBox->value();
    double S = ui->sSpinBox->value();
    double cD = ui->cDSpinBox->value();
    double lambda = ui->lambdaSpinBox->value();
    int mf = ui->mfSpinBox->value();
    int nf = ui->nfSpinBox->value();
    double Xf = ui->xfSpinBox->value();
    double yy = ui->yySpinBox->value();
    double y_len = ui->ySpinBox->value();
    int N = ui->nSpinBox->value();

    int steps = 100;
    res_tD.clear(); res_pD.clear(); res_dpD.clear(); res_td_dpD.clear();
    res_tD.resize(steps); res_pD.resize(steps);

    for (int k = 0; k < steps; ++k) {
        double exp = -7.0 + (6.0 - (-7.0)) * k / (steps - 1);
        res_tD[k] = pow(10.0, exp);
    }

    // 主循环：计算压力
    for (int k = 0; k < steps; ++k) {
        // 防止卡死：每处理一个时间点刷新一次界面
        if (k % 2 == 0) QCoreApplication::processEvents();

        double pd_val = 0.0;
        double current_tD = res_tD[k];
        double ln2 = log(2.0);

        for (int m = 1; m <= N; ++m) {
            double s = m * ln2 / current_tD;
            double L = flaplace(s, cD, S, mf, nf, omega, lambda, Xf, yy, y_len, N);
            double Vi = getStefestVi(m, N);
            pd_val += Vi * ln2 * L / current_tD;
        }
        res_pD[k] = pd_val;
    }

    // 计算压力导数 (MATLAB: td=tD/cD; dpd=td(2:end).*diff(pd)./diff(td))
    QVector<double> plot_tD_CD;
    for(double val : res_tD) plot_tD_CD.append(val / cD);

    for (int k = 0; k < steps - 1; ++k) {
        double t_next = plot_tD_CD[k+1];
        double t_curr = plot_tD_CD[k];
        double p_next = res_pD[k+1];
        double p_curr = res_pD[k];

        if (t_next != t_curr) {
            double val = t_next * (p_next - p_curr) / (t_next - t_curr);
            res_dpD.append(val);
            res_td_dpD.append(t_next);
        }
    }

    m_chartWidget->setData(plot_tD_CD, res_pD, res_dpD, res_td_dpD);

    QString resultText = QString("计算完成\n点数: %1\n").arg(steps);
    resultText += "tD/CD\t\tPD\t\tdPD\n";
    for(int i=0; i<10 && i<res_pD.size(); ++i) {
        double dp = (i < res_dpD.size()) ? res_dpD[i] : 0.0;
        resultText += QString("%1\t%2\t%3\n").arg(plot_tD_CD[i],0,'e',4).arg(res_pD[i],0,'e',4).arg(dp,0,'e',4);
    }
    ui->resultTextEdit->setText(resultText);

    QMap<QString, double> rMap; rMap["PointCount"] = steps;
    emit calculationCompleted("MFHW_DualPorosity", rMap);
}

double ModelWidget1::flaplace(double z, double cD, double S, int mf, int nf,
                              double omega, double lambda, double Xf, double yy, double y, int)
{
    int size = mf * 2 * nf;
    QVector<QVector<double>> E(size + 1, QVector<double>(size + 1));
    QVector<double> F(size + 1, 0.0);

    double fz = (z * (omega * z * (1.0 - omega) + lambda)) / (lambda + (1.0 - omega) * z);

    for (int i = 1; i <= mf; ++i) {
        double yij = yy + (y - 2*yy)/(mf-1.0) * (i-1);
        double yDij = yij / y;

        for (int j = 1; j <= 2 * nf; ++j) {
            int row = (i - 1) * 2 * nf + j - 1;

            double xij = (j - nf - 1.0) / nf * Xf;
            double xDij = xij / y;
            double xij1 = (j - nf) / (double)nf * Xf;
            double xDij1 = xij1 / y;

            for (int k = 1; k <= mf; ++k) {
                double Ykv = yy + (y - 2*yy)/(mf-1.0) * (k-1);
                double YDkv = Ykv / y;

                for (int v = 1; v <= 2 * nf; ++v) {
                    int col = (k - 1) * 2 * nf + v - 1;
                    double Xkv = (2.0*v - 2*nf - 1.0) / (2.0*nf) * Xf;
                    double XDkv = Xkv / y;

                    double val = integral_bessel(XDkv, YDkv, yDij, fz, xDij, xDij1);
                    E[row][col] = val;
                }
            }
            E[row][size] = -1.0;
            E[size][row] = xDij1 - xDij;
        }
    }
    E[size][size] = 0.0;
    F[size] = 1.0 / z;

    QVector<double> res = gaussElimination(E, F);
    if(res.size() > size) return res[size];
    return 0.0;
}

// ------------------------------------------------------------------------
// 高斯-勒让德积分 (优化版：远场近似 + 递归限深)
// ------------------------------------------------------------------------
static const double GL15_X[] = { 0.0,
                                0.2011940939974345, 0.3941513470775634, 0.5709721726085388, 0.7244177313601701,
                                0.8482065834104272, 0.9372985251687639, 0.9879925180204854 };
static const double GL15_W[] = { 0.2025782419255613,
                                0.1984314853271116, 0.1861610000155622, 0.1662692058169939, 0.1395706779049514,
                                0.1071592204671719, 0.0703660474881081, 0.0307532419961173 };

double ModelWidget1::gauss15(std::function<double(double)> f, double a, double b) {
    double halfLen = 0.5 * (b - a);
    double center = 0.5 * (a + b);
    double sum = GL15_W[0] * f(center);
    for (int i = 1; i < 8; ++i) {
        double dx = halfLen * GL15_X[i];
        sum += GL15_W[i] * (f(center - dx) + f(center + dx));
    }
    return sum * halfLen;
}

double ModelWidget1::adaptiveGauss(std::function<double(double)> f, double a, double b, double eps, int depth, int maxDepth) {
    double c = (a + b) / 2.0;
    double v1 = gauss15(f, a, b);
    double v2 = gauss15(f, a, c) + gauss15(f, c, b);

    // 如果达到最大深度，直接返回细分结果，防止卡死
    if (depth >= maxDepth) return v2;

    // 稍微放宽容差以提高速度 (1e-6 工程上足够)
    if (std::abs(v1 - v2) < 1e-10 * std::abs(v2) + eps) return v2;

    return adaptiveGauss(f, a, c, eps/2, depth+1, maxDepth) + adaptiveGauss(f, c, b, eps/2, depth+1, maxDepth);
}

double ModelWidget1::integral_bessel(double XDkv, double YDkv, double yDij, double fz, double a, double b) {
    double sqrt_fz = sqrt(fz);
    double dist_y_sq = (YDkv - yDij) * (YDkv - yDij);

    auto func = [=](double xwD) -> double {
        double dist_x = XDkv - xwD;
        double arg = sqrt(dist_x * dist_x + dist_y_sq) * sqrt_fz;
        return besselK0(arg);
    };

    // 性能优化：远场判断
    // 如果场点到源线段中心的距离远大于线段长度的一半，则认为在远场，直接积分不递归
    double segmentLen = b - a;
    double distToCenter = sqrt(pow(XDkv - (a+b)/2.0, 2) + dist_y_sq);

    if (distToCenter > 2.0 * segmentLen) {
        // 远场：使用高阶高斯积分 (30点)
        return gauss15(func, a, (a+b)/2.0) + gauss15(func, (a+b)/2.0, b);
    }

    // 奇异性处理：源和场重合/极近
    if (dist_y_sq < 1e-16) {
        if (XDkv > a + 1e-9 && XDkv < b - 1e-9) {
            // 奇异点在内部：拆分 + 限制深度(12层足够)
            return adaptiveGauss(func, a, XDkv, 1e-7, 0, 12) + adaptiveGauss(func, XDkv, b, 1e-7, 0, 12);
        }
    }

    // 近场非奇异：限制深度
    return adaptiveGauss(func, a, b, 1e-7, 0, 12);
}

// ------------------------------------------------------------------------
// 通用数学工具
// ------------------------------------------------------------------------

QVector<double> ModelWidget1::gaussElimination(QVector<QVector<double>> A, QVector<double> b) {
    int n = b.size();
    for(int i=0; i<n; ++i) A[i].append(b[i]);

    for (int k = 0; k < n - 1; ++k) {
        if (std::abs(A[k][k]) < 1e-12) continue;
        for (int i = k + 1; i < n; ++i) {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j < n + 1; ++j) A[i][j] -= factor * A[k][j];
        }
    }
    QVector<double> x(n);
    if (std::abs(A[n-1][n-1]) > 1e-12) x[n-1] = A[n-1][n] / A[n-1][n-1];
    for (int i = n - 2; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) sum += A[i][j] * x[j];
        if (std::abs(A[i][i]) > 1e-12) x[i] = (A[i][n] - sum) / A[i][i];
    }
    return x;
}

double ModelWidget1::besselK0(double x) {
    if (x <= 1e-15) return 50.0;
#if __cplusplus >= 201703L
    return std::cyl_bessel_k(0, x);
#else
    if (x <= 2.0) {
        double y = x * x / 4.0;
        return (-log(x / 2.0) * std::cyl_bessel_i(0, x)) +
               (0.42278420 + y * (0.23069756 + y * (0.03488590 + y * (0.00262698 + y * (0.00010750 + y * 0.00000740)))));
    } else {
        double y = 2.0 / x;
        return (exp(-x) / sqrt(x)) * (1.25331414 + y * (-0.07832358 + y * (0.02189568 + y * (-0.01062446 + y * (0.00587872 + y * (-0.00251540 + y * 0.00053208))))));
    }
#endif
}

double ModelWidget1::factorial(int n) {
    if (n <= 1) return 1.0;
    double res = 1.0;
    for (int i = 2; i <= n; ++i) res *= i;
    return res;
}

double ModelWidget1::getStefestVi(int i, int N) {
    double sum = 0.0;
    int k_start = (i + 1) / 2;
    int k_end = std::min(i, N / 2);
    for (int k = k_start; k <= k_end; ++k) {
        double num = pow(k, N / 2.0) * factorial(2 * k);
        double den = factorial(N / 2 - k) * factorial(k) * factorial(k - 1) * factorial(i - k) * factorial(2 * k - i);
        sum += num / den;
    }
    return ((i + N / 2) % 2 == 0 ? 1.0 : -1.0) * sum;
}

void ModelWidget1::onExportResults() {
    if (res_tD.isEmpty()) return;
    QString path = QFileDialog::getSaveFileName(this, "导出CSV", "", "CSV Files (*.csv)");
    if (path.isEmpty()) return;
    QFile f(path);
    if (f.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QTextStream out(&f);
        double cD = ui->cDSpinBox->value();
        out << "tD,tD/CD,PD,dPD\n";
        for (int i = 0; i < res_tD.size(); ++i) {
            double dp = (i < res_dpD.size()) ? res_dpD[i] : 0.0;
            out << res_tD[i] << "," << res_tD[i]/cD << "," << res_pD[i] << "," << dp << "\n";
        }
        f.close();
        QMessageBox::information(this, "导出成功", "文件已保存");
    }
}
