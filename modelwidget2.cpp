#include "modelwidget2.h"
#include "ui_modelwidget2.h"
#include <QMessageBox>
#include <QFileDialog>
#include <QTextStream>
#include <QDebug>
#include <QDateTime>
#include <cmath>
#include <algorithm>
#include <numeric>

// ===============================
// SimpleSplineFD 实现
// ===============================

void SimpleSplineFD::setData(const QVector<double>& x, const QVector<double>& y)
{
    if (x.size() != y.size() || x.size() < 3) {
        m_isValid = false;
        return;
    }

    // 确保x是递增的
    QVector<int> indices(x.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](int i, int j) {
        return x[i] < x[j];
    });

    m_x.resize(x.size());
    m_y.resize(y.size());
    for (int i = 0; i < indices.size(); ++i) {
        m_x[i] = x[indices[i]];
        m_y[i] = y[indices[i]];
    }

    computeCoefficients();
    m_isValid = true;
}

void SimpleSplineFD::computeCoefficients()
{
    int n = m_x.size();
    m_a = m_y;
    m_b.resize(n);
    m_c.resize(n);
    m_d.resize(n);

    // 计算步长
    QVector<double> h(n-1);
    for (int i = 0; i < n-1; ++i) {
        h[i] = m_x[i+1] - m_x[i];
        if (h[i] <= 0) {
            m_isValid = false;
            return;
        }
    }

    // 自然边界条件的三次样条
    QVector<double> alpha(n-1);
    for (int i = 1; i < n-1; ++i) {
        alpha[i] = 3.0 * ((m_y[i+1] - m_y[i]) / h[i] - (m_y[i] - m_y[i-1]) / h[i-1]);
    }

    QVector<double> l(n), mu(n), z(n);
    l[0] = 1.0;
    mu[0] = 0.0;
    z[0] = 0.0;

    for (int i = 1; i < n-1; ++i) {
        l[i] = 2.0 * (m_x[i+1] - m_x[i-1]) - h[i-1] * mu[i-1];
        if (qAbs(l[i]) < 1e-12) {
            m_isValid = false;
            return;
        }
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i-1] * z[i-1]) / l[i];
    }

    l[n-1] = 1.0;
    z[n-1] = 0.0;
    m_c[n-1] = 0.0;

    for (int j = n-2; j >= 0; --j) {
        m_c[j] = z[j] - mu[j] * m_c[j+1];
        m_b[j] = (m_y[j+1] - m_y[j]) / h[j] - h[j] * (m_c[j+1] + 2.0 * m_c[j]) / 3.0;
        m_d[j] = (m_c[j+1] - m_c[j]) / (3.0 * h[j]);
    }
}

int SimpleSplineFD::findSegment(double x) const
{
    if (x <= m_x.first()) return 0;
    if (x >= m_x.last()) return m_x.size() - 2;

    int left = 0, right = m_x.size() - 1;
    while (left < right - 1) {
        int mid = (left + right) / 2;
        if (x < m_x[mid]) {
            right = mid;
        } else {
            left = mid;
        }
    }
    return left;
}

double SimpleSplineFD::interpolate(double x) const
{
    if (!m_isValid) return 0.0;

    int i = findSegment(x);
    double dx = x - m_x[i];

    return m_a[i] + m_b[i] * dx + m_c[i] * dx * dx + m_d[i] * dx * dx * dx;
}

QVector<double> SimpleSplineFD::interpolateRange(const QVector<double>& xNew) const
{
    QVector<double> result(xNew.size());
    for (int i = 0; i < xNew.size(); ++i) {
        result[i] = interpolate(xNew[i]);
    }
    return result;
}

// ===============================
// LogLogChartWidgetFD 实现
// ===============================

LogLogChartWidgetFD::LogLogChartWidgetFD(QWidget *parent)
    : QWidget(parent)
    , m_xMin(1e-3), m_xMax(1e13)
    , m_yMin(1e-3), m_yMax(1e1)
    , m_hasData(false)
    , m_showOriginalData(true)
    , m_showSmoothedCurve(true)
    , m_isDragging(false)
{
    setMinimumSize(600, 400);
    setStyleSheet("QWidget { background-color: white; border: 1px solid gray; }");
    setMouseTracking(true);
}

void LogLogChartWidgetFD::setData(const QVector<double>& xData, const QVector<double>& yData1,
                                  const QVector<double>& yData2, const QVector<double>& xData2)
{
    m_xData = xData;
    m_yData1 = yData1;
    m_yData2 = yData2;
    m_xData2 = xData2.isEmpty() ? xData : xData2;
    m_hasData = !xData.isEmpty() && !yData1.isEmpty();

    if (m_hasData) {
        autoFitData();
    }
    update();
}

void LogLogChartWidgetFD::clearData()
{
    m_hasData = false;
    m_xData.clear();
    m_yData1.clear();
    m_yData2.clear();
    m_xData2.clear();
    resetView();
    update();
}

void LogLogChartWidgetFD::resetView()
{
    if (m_hasData) {
        autoFitData();
    } else {
        m_xMin = 1e-3;
        m_xMax = 1e13;
        m_yMin = 1e-3;
        m_yMax = 1e1;
    }
    update();
}

void LogLogChartWidgetFD::autoFitData()
{
    if (!m_hasData) return;

    // 收集所有有效数据
    QVector<double> allX, allY;

    for (double x : m_xData) {
        if (x > 0 && std::isfinite(x)) allX.append(x);
    }
    for (double x : m_xData2) {
        if (x > 0 && std::isfinite(x)) allX.append(x);
    }
    for (double y : m_yData1) {
        if (y > 0 && std::isfinite(y)) allY.append(y);
    }
    for (double y : m_yData2) {
        if (y > 0 && std::isfinite(y)) allY.append(y);
    }

    if (allX.isEmpty() || allY.isEmpty()) return;

    // 找到数据范围
    auto [xMinIt, xMaxIt] = std::minmax_element(allX.begin(), allX.end());
    auto [yMinIt, yMaxIt] = std::minmax_element(allY.begin(), allY.end());

    double xMin = *xMinIt;
    double xMax = *xMaxIt;
    double yMin = *yMinIt;
    double yMax = *yMaxIt;

    // 在对数空间添加边距
    double logXMin = log10(xMin);
    double logXMax = log10(xMax);
    double logYMin = log10(yMin);
    double logYMax = log10(yMax);

    double xMargin = (logXMax - logXMin) * 0.05;
    double yMargin = (logYMax - logYMin) * 0.05;

    m_xMin = pow(10.0, logXMin - xMargin);
    m_xMax = pow(10.0, logXMax + xMargin);
    m_yMin = pow(10.0, logYMin - yMargin);
    m_yMax = pow(10.0, logYMax + yMargin);
}

void LogLogChartWidgetFD::mousePressEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton) {
        m_isDragging = true;
        m_lastMousePos = event->pos();
        setCursor(Qt::ClosedHandCursor);
    }
}

void LogLogChartWidgetFD::mouseMoveEvent(QMouseEvent *event)
{
    if (m_isDragging) {
        QPoint delta = event->pos() - m_lastMousePos;
        m_lastMousePos = event->pos();

        QRect plotRect = rect().adjusted(80, 50, -50, -80);

        double logXRange = log10(m_xMax) - log10(m_xMin);
        double logYRange = log10(m_yMax) - log10(m_yMin);

        double deltaLogX = -delta.x() * logXRange / plotRect.width();
        double deltaLogY = delta.y() * logYRange / plotRect.height();

        m_xMin = pow(10.0, log10(m_xMin) + deltaLogX);
        m_xMax = pow(10.0, log10(m_xMax) + deltaLogX);
        m_yMin = pow(10.0, log10(m_yMin) + deltaLogY);
        m_yMax = pow(10.0, log10(m_yMax) + deltaLogY);

        update();
    } else {
        setCursor(Qt::OpenHandCursor);
    }
}

void LogLogChartWidgetFD::mouseReleaseEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton) {
        m_isDragging = false;
        setCursor(Qt::OpenHandCursor);
    }
}

void LogLogChartWidgetFD::wheelEvent(QWheelEvent *event)
{
    const double zoomInFactor = 1.15;
    const double zoomOutFactor = 1.0 / zoomInFactor;

    double zoomFactor = (event->angleDelta().y() > 0) ? zoomInFactor : zoomOutFactor;

    QRect plotRect = rect().adjusted(80, 50, -50, -80);
    QPointF mousePos = event->position();

    double relativeX = (mousePos.x() - plotRect.left()) / plotRect.width();
    double relativeY = (plotRect.bottom() - mousePos.y()) / plotRect.height();

    relativeX = qBound(0.0, relativeX, 1.0);
    relativeY = qBound(0.0, relativeY, 1.0);

    double logXRange = log10(m_xMax) - log10(m_xMin);
    double logYRange = log10(m_yMax) - log10(m_yMin);

    double mouseLogX = log10(m_xMin) + relativeX * logXRange;
    double mouseLogY = log10(m_yMin) + relativeY * logYRange;

    double newLogXRange = logXRange / zoomFactor;
    double newLogYRange = logYRange / zoomFactor;

    m_xMin = pow(10.0, mouseLogX - relativeX * newLogXRange);
    m_xMax = pow(10.0, mouseLogX + (1.0 - relativeX) * newLogXRange);
    m_yMin = pow(10.0, mouseLogY - relativeY * newLogYRange);
    m_yMax = pow(10.0, mouseLogY + (1.0 - relativeY) * newLogYRange);

    update();
    event->accept();
}

void LogLogChartWidgetFD::paintEvent(QPaintEvent *event)
{
    Q_UNUSED(event)

    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);

    QRect plotRect = rect().adjusted(80, 50, -50, -80);

    painter.fillRect(rect(), Qt::white);
    painter.fillRect(plotRect, QColor(250, 250, 250));

    if (!m_title.isEmpty()) {
        painter.setFont(QFont("Arial", 14, QFont::Bold));
        painter.setPen(Qt::black);
        painter.drawText(rect().adjusted(0, 10, 0, 0), Qt::AlignHCenter | Qt::AlignTop, m_title);
    }

    drawAxis(painter, plotRect);

    if (m_hasData) {
        drawData(painter, plotRect);
        drawLegend(painter, plotRect);
    } else {
        painter.setFont(QFont("Arial", 12));
        painter.setPen(Qt::gray);
        painter.drawText(plotRect, Qt::AlignCenter, "点击'开始计算'生成压力分析图表");
    }
}

void LogLogChartWidgetFD::drawAxis(QPainter& painter, const QRect& plotRect)
{
    painter.setPen(QPen(Qt::black, 2));
    painter.drawRect(plotRect);

    painter.setPen(QPen(Qt::lightGray, 1));
    painter.setFont(QFont("Arial", 9));

    // X轴对数刻度
    double logXMin = log10(m_xMin);
    double logXMax = log10(m_xMax);
    double xLogSpan = logXMax - logXMin;

    int xStep = (xLogSpan > 10) ? 2 : 1;
    int startExp = (int)floor(logXMin);
    int endExp = (int)ceil(logXMax);

    for (int exp = startExp; exp <= endExp; exp += xStep) {
        double x = pow(10.0, exp);
        if (x >= m_xMin && x <= m_xMax) {
            QPointF point = dataToPixel(x, m_yMin, plotRect);
            if (point.x() >= plotRect.left() && point.x() <= plotRect.right()) {
                painter.drawLine(point.x(), plotRect.bottom(), point.x(), plotRect.top());

                painter.setPen(Qt::black);
                QString label = QString("1e%1").arg(exp);
                QRect labelRect(point.x() - 25, plotRect.bottom() + 5, 50, 20);
                painter.drawText(labelRect, Qt::AlignCenter, label);
                painter.setPen(Qt::lightGray);
            }
        }
    }

    // Y轴对数刻度
    double logYMin = log10(m_yMin);
    double logYMax = log10(m_yMax);
    double yLogSpan = logYMax - logYMin;

    int yStep = (yLogSpan > 8) ? 2 : 1;
    startExp = (int)floor(logYMin);
    endExp = (int)ceil(logYMax);

    for (int exp = startExp; exp <= endExp; exp += yStep) {
        double y = pow(10.0, exp);
        if (y >= m_yMin && y <= m_yMax) {
            QPointF point = dataToPixel(m_xMin, y, plotRect);
            if (point.y() >= plotRect.top() && point.y() <= plotRect.bottom()) {
                painter.drawLine(plotRect.left(), point.y(), plotRect.right(), point.y());

                painter.setPen(Qt::black);
                QString label = QString("1e%1").arg(exp);
                QRect labelRect(plotRect.left() - 50, point.y() - 10, 45, 20);
                painter.drawText(labelRect, Qt::AlignRight | Qt::AlignVCenter, label);
                painter.setPen(Qt::lightGray);
            }
        }
    }

    // 坐标轴标签
    painter.setPen(Qt::black);
    painter.setFont(QFont("Arial", 11, QFont::Bold));
    painter.drawText(plotRect.center().x() - 30, plotRect.bottom() + 40, "tD/CD");

    painter.save();
    painter.translate(plotRect.left() - 60, plotRect.center().y());
    painter.rotate(-90);
    painter.drawText(-50, 0, "PD & dPD/ln(tD)");
    painter.restore();
}

void LogLogChartWidgetFD::drawData(QPainter& painter, const QRect& plotRect)
{
    // 绘制压力曲线
    if (!m_xData.isEmpty() && !m_yData1.isEmpty()) {
        painter.setPen(QPen(Qt::red, 2));
        QVector<QPointF> points;

        for (int i = 0; i < qMin(m_xData.size(), m_yData1.size()); ++i) {
            if (m_xData[i] > 0 && m_yData1[i] > 0) {
                QPointF point = dataToPixel(m_xData[i], m_yData1[i], plotRect);
                if (plotRect.contains(point.toPoint())) {
                    points.append(point);
                }
            }
        }

        for (int i = 1; i < points.size(); ++i) {
            painter.drawLine(points[i-1], points[i]);
        }

        // 绘制数据点
        if (m_showOriginalData) {
            painter.setPen(QPen(Qt::red, 1));
            painter.setBrush(QBrush(Qt::red));
            for (const QPointF& point : points) {
                painter.drawEllipse(point, 3, 3);
            }
        }
    }

    // 绘制压力导数曲线
    if (!m_xData2.isEmpty() && !m_yData2.isEmpty()) {
        painter.setPen(QPen(Qt::blue, 2));
        QVector<QPointF> points;

        for (int i = 0; i < qMin(m_xData2.size(), m_yData2.size()); ++i) {
            if (m_xData2[i] > 0 && m_yData2[i] > 0) {
                QPointF point = dataToPixel(m_xData2[i], m_yData2[i], plotRect);
                if (plotRect.contains(point.toPoint())) {
                    points.append(point);
                }
            }
        }

        for (int i = 1; i < points.size(); ++i) {
            painter.drawLine(points[i-1], points[i]);
        }

        // 绘制数据点
        if (m_showOriginalData) {
            painter.setPen(QPen(Qt::blue, 1));
            painter.setBrush(QBrush(Qt::blue));
            for (const QPointF& point : points) {
                painter.drawEllipse(point, 3, 3);
            }
        }
    }
}

void LogLogChartWidgetFD::drawLegend(QPainter& painter, const QRect& plotRect)
{
    painter.setFont(QFont("Arial", 10));

    int legendX = plotRect.right() - 160;
    int legendY = plotRect.top() + 20;
    int lineSpacing = 20;

    // 压力曲线图例
    painter.setPen(QPen(Qt::red, 2));
    painter.drawLine(legendX, legendY, legendX + 30, legendY);
    painter.setPen(Qt::black);
    painter.drawText(legendX + 35, legendY + 5, "压力");

    // 压力导数图例
    if (!m_yData2.isEmpty()) {
        legendY += lineSpacing;
        painter.setPen(QPen(Qt::blue, 2));
        painter.drawLine(legendX, legendY, legendX + 30, legendY);
        painter.setPen(Qt::black);
        painter.drawText(legendX + 35, legendY + 5, "压力导数");
    }
}

QPointF LogLogChartWidgetFD::dataToPixel(double x, double y, const QRect& plotRect)
{
    x = qMax(x, 1e-15);
    y = qMax(y, 1e-15);

    double logX = log10(x);
    double logY = log10(y);

    double logXMin = log10(qMax(m_xMin, 1e-15));
    double logXMax = log10(qMax(m_xMax, 1e-15));
    double logYMin = log10(qMax(m_yMin, 1e-15));
    double logYMax = log10(qMax(m_yMax, 1e-15));

    double pixelX = plotRect.left() + (logX - logXMin) / (logXMax - logXMin) * plotRect.width();
    double pixelY = plotRect.bottom() - (logY - logYMin) / (logYMax - logYMin) * plotRect.height();

    return QPointF(pixelX, pixelY);
}

// ===============================
// ModelWidget2 实现
// ===============================

ModelWidget2::ModelWidget2(QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::ModelWidget2)
    , m_chartWidget(nullptr)
    , m_calculationWatcher(nullptr)
    , m_progressTimer(new QTimer(this))
    , m_isCalculating(false)
{
    ui->setupUi(this);

    // 创建图表组件
    m_chartWidget = new LogLogChartWidgetFD();
    QVBoxLayout* chartLayout = new QVBoxLayout(ui->chartTab);
    chartLayout->addWidget(m_chartWidget);

    resetParametersToDefault();

    m_calculationWatcher = new QFutureWatcher<CalculationResultFD>(this);
    connect(m_calculationWatcher, &QFutureWatcher<CalculationResultFD>::finished,
            this, &ModelWidget2::onCalculationFinished);

    connect(m_progressTimer, &QTimer::timeout, this, [this]() {
        static int value = 0;
        value = (value + 1) % 100;
        if (ui->progressBar) {
            ui->progressBar->setValue(value);
        }
    });

    // 连接参数变化信号
    connect(ui->omegaSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            this, &ModelWidget2::onParameterChanged);
    connect(ui->sSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            this, &ModelWidget2::onParameterChanged);
    connect(ui->cDSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            this, &ModelWidget2::onParameterChanged);
    connect(ui->lambdaSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            this, &ModelWidget2::onParameterChanged);
    connect(ui->mfSpinBox, QOverload<int>::of(&QSpinBox::valueChanged),
            this, &ModelWidget2::onParameterChanged);
    connect(ui->nfSpinBox, QOverload<int>::of(&QSpinBox::valueChanged),
            this, &ModelWidget2::onParameterChanged);
    connect(ui->xfSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            this, &ModelWidget2::onParameterChanged);
    connect(ui->yySpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            this, &ModelWidget2::onParameterChanged);
    connect(ui->ySpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            this, &ModelWidget2::onParameterChanged);
    connect(ui->nSpinBox, QOverload<int>::of(&QSpinBox::valueChanged),
            this, &ModelWidget2::onParameterChanged);
    connect(ui->cfdSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            this, &ModelWidget2::onParameterChanged);
    connect(ui->kpdSpinBox, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
            this, &ModelWidget2::onParameterChanged);
    connect(ui->enableFittingCheckBox, &QCheckBox::toggled,
            this, &ModelWidget2::onSmoothingEnabledChanged);
    connect(ui->fittingPointsSpinBox, QOverload<int>::of(&QSpinBox::valueChanged),
            this, &ModelWidget2::onParameterChanged);
    connect(ui->originalPointsSpinBox, QOverload<int>::of(&QSpinBox::valueChanged),
            this, &ModelWidget2::onParameterChanged);

    // 连接按钮信号
    connect(ui->calculateButton, &QPushButton::clicked, this, &ModelWidget2::onCalculateClicked);
    connect(ui->resetButton, &QPushButton::clicked, this, &ModelWidget2::onResetParameters);
    connect(ui->exportButton, &QPushButton::clicked, this, &ModelWidget2::onExportResults);
    connect(ui->resetViewButton, &QPushButton::clicked, this, &ModelWidget2::onResetView);
    connect(ui->fitToDataButton, &QPushButton::clicked, this, &ModelWidget2::onFitToData);
    connect(ui->showFittedCurveCheckBox, &QCheckBox::toggled,
            this, &ModelWidget2::onShowSmoothedCurveChanged);
    connect(ui->showOriginalDataCheckBox, &QCheckBox::toggled,
            this, &ModelWidget2::onShowOriginalDataChanged);
}

ModelWidget2::~ModelWidget2()
{
    if (m_calculationWatcher && m_calculationWatcher->isRunning()) {
        m_calculationWatcher->cancel();
        m_calculationWatcher->waitForFinished();
    }
    delete ui;
}

void ModelWidget2::onCalculateClicked()
{
    if (m_isCalculating) {
        QMessageBox::information(this, "提示", "计算正在进行中，请等待...");
        return;
    }

    onParameterChanged();
    setCalculationInProgress(true);

    QFuture<CalculationResultFD> future = QtConcurrent::run([this]() {
        return performCalculation(m_currentParams);
    });

    m_calculationWatcher->setFuture(future);
}

void ModelWidget2::onParameterChanged()
{
    if (!ui->omegaSpinBox) return;

    m_currentParams.omega = ui->omegaSpinBox->value();
    m_currentParams.S = ui->sSpinBox->value();
    m_currentParams.cD = ui->cDSpinBox->value();
    m_currentParams.lambda = ui->lambdaSpinBox->value();
    m_currentParams.mf = ui->mfSpinBox->value();
    m_currentParams.nf = ui->nfSpinBox->value();
    m_currentParams.Xf = ui->xfSpinBox->value();
    m_currentParams.yy = ui->yySpinBox->value();
    m_currentParams.y = ui->ySpinBox->value();
    m_currentParams.N = ui->nSpinBox->value();
    m_currentParams.CFD = ui->cfdSpinBox->value();
    m_currentParams.kpd = ui->kpdSpinBox->value();

    m_currentParams.enableSmoothing = ui->enableFittingCheckBox->isChecked();
    m_currentParams.fittingPoints = ui->fittingPointsSpinBox->value();
    m_currentParams.originalPoints = ui->originalPointsSpinBox->value();
}

void ModelWidget2::onExportResults()
{
    if (!m_lastResult.isValid || m_lastResult.tD.isEmpty()) {
        QMessageBox::warning(this, "警告", "没有可导出的计算结果，请先进行计算。");
        return;
    }

    QString fileName = QFileDialog::getSaveFileName(this, "导出计算结果",
                                                    "有限导流模型分析结果.csv",
                                                    "CSV文件 (*.csv)");
    if (fileName.isEmpty()) {
        return;
    }

    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QMessageBox::critical(this, "错误", "无法创建文件：" + fileName);
        return;
    }

    QTextStream out(&file);

    // 写入标题
    out << "tD,tD/CD,原始压力(PD),压力敏感压力(PD1),压力导数时间(td_dpd/CD),压力导数(dPD/ln(tD))";
    if (m_lastResult.hasSmoothedData) {
        out << ",tD_平滑,tD_平滑/CD,压力_平滑,压力导数时间_平滑,压力导数_平滑";
    }
    out << "\n";

    // 写入数据
    int maxRows = qMax(m_lastResult.tD.size(),
                       m_lastResult.hasSmoothedData ? m_lastResult.tD_smooth.size() : 0);

    for (int i = 0; i < maxRows; ++i) {
        // 原始数据
        if (i < m_lastResult.tD.size()) {
            double tDOverCD = m_lastResult.tD[i] / m_currentParams.cD;
            out << QString("%1,%2,%3,%4,")
                       .arg(m_lastResult.tD[i], 0, 'e', 6)
                       .arg(tDOverCD, 0, 'e', 6)
                       .arg(m_lastResult.pd[i], 0, 'e', 6)
                       .arg(m_lastResult.pd1[i], 0, 'e', 6);

            if (i < m_lastResult.td_dpd.size() && i < m_lastResult.dpd.size()) {
                out << QString("%1,%2")
                .arg(m_lastResult.td_dpd[i], 0, 'e', 6)
                    .arg(m_lastResult.dpd[i], 0, 'e', 6);
            } else {
                out << ",";
            }
        } else {
            out << ",,,,,";
        }

        // 平滑数据
        if (m_lastResult.hasSmoothedData) {
            if (i < m_lastResult.tD_smooth.size()) {
                double tDOverCD_smooth = m_lastResult.tD_smooth[i] / m_currentParams.cD;
                out << QString(",%1,%2,%3,")
                           .arg(m_lastResult.tD_smooth[i], 0, 'e', 6)
                           .arg(tDOverCD_smooth, 0, 'e', 6)
                           .arg(m_lastResult.pd_smooth[i], 0, 'e', 6);

                if (i < m_lastResult.td_dpd_smooth.size() && i < m_lastResult.dpd_smooth.size()) {
                    out << QString("%1,%2")
                    .arg(m_lastResult.td_dpd_smooth[i], 0, 'e', 6)
                        .arg(m_lastResult.dpd_smooth[i], 0, 'e', 6);
                } else {
                    out << ",";
                }
            } else {
                out << ",,,,";
            }
        }

        out << "\n";
    }

    QMessageBox::information(this, "成功", "结果已导出到：" + fileName);
}

void ModelWidget2::onResetParameters()
{
    resetParametersToDefault();
}

void ModelWidget2::onCalculationFinished()
{
    setCalculationInProgress(false);

    CalculationResultFD result = m_calculationWatcher->result();

    if (result.isValid) {
        m_lastResult = result;
        updateChart(result);
        updateResultText(result);

        ui->exportButton->setEnabled(true);
        ui->resetViewButton->setEnabled(true);
        ui->fitToDataButton->setEnabled(true);
        ui->showOriginalDataCheckBox->setEnabled(result.hasSmoothedData);
        ui->showFittedCurveCheckBox->setEnabled(result.hasSmoothedData);

        QMessageBox::information(this, "计算完成", "有限导流压力分析计算完成！");

        ui->tabWidget->setCurrentIndex(0);

        // 发出计算完成信号
        QMap<QString, double> results;
        results["点数"] = result.tD.size();
        results["最大压力"] = *std::max_element(result.pd1.begin(), result.pd1.end());
        results["最小压力"] = *std::min_element(result.pd1.begin(), result.pd1.end());
        QString analysisType = "有限导流双重孔隙介质页岩油藏渗流模型";
        emit calculationCompleted(analysisType, results);

    } else {
        QMessageBox::critical(this, "计算错误", "计算过程中发生错误：\n" + result.errorMessage);
        ui->resultTextEdit->setPlainText("计算失败：" + result.errorMessage);
    }
}

void ModelWidget2::onResetView()
{
    if (m_chartWidget) {
        m_chartWidget->resetView();
    }
}

void ModelWidget2::onFitToData()
{
    if (m_chartWidget) {
        m_chartWidget->autoFitData();
        m_chartWidget->update();
    }
}

void ModelWidget2::onShowOriginalDataChanged(bool show)
{
    if (m_chartWidget) {
        m_chartWidget->setShowOriginalData(show);
    }
}

void ModelWidget2::onShowSmoothedCurveChanged(bool show)
{
    if (m_chartWidget) {
        m_chartWidget->setShowSmoothedCurve(show);
    }
}

void ModelWidget2::onSmoothingEnabledChanged(bool enabled)
{
    ui->fittingPointsSpinBox->setEnabled(enabled);
    onParameterChanged();
}


void ModelWidget2::resetParametersToDefault()
{
    if (!ui->omegaSpinBox) return;

    ui->omegaSpinBox->setValue(0.0155211);
    ui->sSpinBox->setValue(0.810876);
    ui->cDSpinBox->setValue(8.08669e-8);
    ui->lambdaSpinBox->setValue(0.083277);
    ui->mfSpinBox->setValue(3);
    ui->nfSpinBox->setValue(5);
    ui->xfSpinBox->setValue(193.0);
    ui->yySpinBox->setValue(295.0);
    ui->ySpinBox->setValue(2758.0);
    ui->nSpinBox->setValue(4);
    ui->cfdSpinBox->setValue(0.9);
    ui->kpdSpinBox->setValue(0.04);
    ui->originalPointsSpinBox->setValue(100);
    ui->enableFittingCheckBox->setChecked(true);
    ui->fittingPointsSpinBox->setValue(100);

    onParameterChanged();
}

void ModelWidget2::updateChart(const CalculationResultFD& result)
{
    if (!m_chartWidget || result.tD.isEmpty() || result.pd1.isEmpty()) {
        return;
    }

    QString title = QString("有限导流双重孔隙介质页岩油藏渗流模型压力分析");

    // 使用平滑后的数据（如果有）
    if (result.hasSmoothedData && !result.tD_smooth.isEmpty()) {
        // 转换为tD/CD坐标
        QVector<double> tDOverCD_smooth;
        for (double tD : result.tD_smooth) {
            tDOverCD_smooth.append(tD / m_currentParams.cD);
        }

        m_chartWidget->setData(tDOverCD_smooth, result.pd_smooth,
                               result.dpd_smooth, result.td_dpd_smooth);
    } else {
        // 使用原始数据（压力敏感处理后的pd1）
        QVector<double> tDOverCD;
        for (double tD : result.tD) {
            tDOverCD.append(tD / m_currentParams.cD);
        }

        m_chartWidget->setData(tDOverCD, result.pd1, result.dpd, result.td_dpd);
    }
}

void ModelWidget2::updateResultText(const CalculationResultFD& result)
{
    if (!ui->resultTextEdit) return;

    QString text = QString("有限导流双重孔隙介质页岩油藏渗流模型压力分析计算结果\n");

    text += QString("\n计算完成时间: %1\n").arg(QDateTime::currentDateTime().toString());
    text += QString("有限导流系数 (CFD): %1\n").arg(m_currentParams.CFD);
    text += QString("压力敏感因子 (kpd): %1\n").arg(m_currentParams.kpd);

    // 添加部分数据点展示
    text += QString("\n原始计算数据前10个点：\n");
    text += QString("序号\t时间tD\t\ttD/CD\t\t原始压力PD\t压力敏感PD1\t\t压力导数dPD\n");
    text += QString("--------------------------------------------------------------------------------\n");

    int maxPoints = qMin(10, result.tD.size());
    for (int i = 0; i < maxPoints; ++i) {
        double tDOverCD = result.tD[i] / m_currentParams.cD;
        double originalPressure = (i < result.pd.size()) ? result.pd[i] : 0.0;
        double sensitivityPressure = (i < result.pd1.size()) ? result.pd1[i] : 0.0;
        double derivative = (i < result.dpd.size()) ? result.dpd[i] : 0.0;

        text += QString("%1\t%2\t%3\t%4\t%5\t%6\n")
                    .arg(i + 1)
                    .arg(result.tD[i], 0, 'e', 3)
                    .arg(tDOverCD, 0, 'e', 3)
                    .arg(originalPressure, 0, 'e', 3)
                    .arg(sensitivityPressure, 0, 'e', 3)
                    .arg(derivative, 0, 'e', 3);
    }

    if (result.tD.size() > 10) {
        text += QString("\n... (共%1个原始数据点，仅显示前10个)\n").arg(result.tD.size());
    }

    text += QString("\n完整数据请查看图表或导出CSV文件。\n");

    ui->resultTextEdit->setPlainText(text);
}

void ModelWidget2::setCalculationInProgress(bool inProgress)
{
    m_isCalculating = inProgress;
    ui->calculateButton->setEnabled(!inProgress);
    ui->calculateButton->setText(inProgress ? "计算中..." : "开始计算");
    ui->progressBar->setVisible(inProgress);

    if (inProgress) {
        m_progressTimer->start(50);
        ui->progressBar->setRange(0, 0);
    } else {
        m_progressTimer->stop();
        ui->progressBar->setRange(0, 100);
        ui->progressBar->setValue(0);
    }
}



CalculationResultFD ModelWidget2::performCalculation(const CalculationParametersFD& params)
{
    CalculationResultFD result;

    try {
        int N = params.N;
        const int numPoints = params.originalPoints;

        // 生成tD序列：logspace(-11,5,numPoints)
        QVector<double> tD;
        tD.reserve(numPoints);
        for (int i = 0; i < numPoints; ++i) {
            double exponent = -11.0 + 16.0 * i / (numPoints - 1);
            double td_val = pow(10.0, exponent);
            tD.append(td_val);
        }

        result.tD = tD;
        result.pd.resize(numPoints);
        result.pd.fill(0.0);

        // 主循环计算
        for (int k = 0; k < numPoints; ++k) {
            double current_tD = tD[k];

            for (int m = 1; m <= N; ++m) {
                double s = m * log(2.0) / current_tD;
                double L = flaplace(s, params);
                double vi = stefestCoefficient(m, N);
                result.pd[k] += vi * log(2.0) * L / current_tD;
            }
        }

        // 压力敏感处理
        result.pd1.resize(numPoints);
        for (int k = 0; k < numPoints; ++k) {
            result.pd1[k] = -1.0 / params.kpd * log(1.0 - params.kpd * result.pd[k]);
        }

        // 使用改进的压力导数计算方法
        computePressureDerivative(result.tD, result.pd1, params.cD,
                                  result.dpd, result.td_dpd);

        // 执行数据平滑（如果可用）
        if (params.enableSmoothing) {
            smoothData(result, params);
        }

        result.isValid = true;

    } catch (const std::exception& e) {
        result.isValid = false;
        result.errorMessage = QString("计算异常: %1").arg(e.what());
    } catch (...) {
        result.isValid = false;
        result.errorMessage = "未知计算错误";
    }

    return result;
}

// 改进的压力导数计算方法（直接修改原函数）
void ModelWidget2::computePressureDerivative(const QVector<double>& tD, const QVector<double>& pd,
                                             double cD, QVector<double>& dpd, QVector<double>& td_dpd)
{
    dpd.clear();
    td_dpd.clear();

    if (tD.size() < 3) return;

    // 转换为td = tD/cD
    QVector<double> td;
    td.reserve(tD.size());
    for (double tD_val : tD) {
        td.append(tD_val / cD);
    }

    // Step 1: 先对压力数据进行预平滑
    QVector<double> smoothedPd = pd;

    // 使用多通道平滑策略
    for (int pass = 0; pass < 3; ++pass) {
        QVector<double> temp = smoothedPd;

        for (int i = 1; i < smoothedPd.size() - 1; ++i) {
            double t = td[i];
            double smoothingFactor;

            if (t < 1e2) {
                smoothingFactor = 0.3;  // 早期轻度平滑
            } else if (t < 1e6) {
                smoothingFactor = 0.5;  // 中期中度平滑
            } else {
                smoothingFactor = 0.7;  // 后期强度平滑
            }

            double avgValue = (temp[i-1] + temp[i] + temp[i+1]) / 3.0;
            smoothedPd[i] = temp[i] * (1 - smoothingFactor) + avgValue * smoothingFactor;
        }
    }

    // Step 2: 在对数空间计算导数
    QVector<double> log_td, log_pd;
    for (int i = 0; i < td.size(); ++i) {
        if (td[i] > 0 && smoothedPd[i] > 0) {
            log_td.append(log(td[i]));
            log_pd.append(log(smoothedPd[i]));
        }
    }

    if (log_td.size() < 3) {
        // 回退到简单差分
        for (int i = 1; i < td.size(); ++i) {
            double td_curr = td[i];
            double td_prev = td[i-1];
            double pd_curr = smoothedPd[i];
            double pd_prev = smoothedPd[i-1];

            if (td_curr > td_prev && td_prev > 0) {
                double dpd_val = td_curr * (pd_curr - pd_prev) / (td_curr - td_prev);
                dpd.append(dpd_val);
                td_dpd.append(td_curr);
            }
        }
        return;
    }

    // Step 3: 使用改进的数值导数方法
    QVector<double> dlogP_dlogt;
    dlogP_dlogt.reserve(log_td.size());

    // 使用五点公式计算导数
    for (int i = 0; i < log_td.size(); ++i) {
        double derivative = 0.0;

        if (i == 0) {
            // 前向差分
            if (log_td.size() >= 2) {
                derivative = (log_pd[1] - log_pd[0]) / (log_td[1] - log_td[0]);
            }
        } else if (i == 1 && log_td.size() >= 3) {
            // 三点公式
            double h0 = log_td[1] - log_td[0];
            double h1 = log_td[2] - log_td[1];
            double a = -(2*h1 + h0) / (h0 * (h0 + h1));
            double b = (h0 + h1) / (h0 * h1);
            double c = -h0 / (h1 * (h0 + h1));
            derivative = a * log_pd[0] + b * log_pd[1] + c * log_pd[2];
        } else if (i >= 2 && i < log_td.size() - 2) {
            // 五点中心差分公式
            double h = log_td[i+1] - log_td[i];
            if (qAbs(log_td[i] - log_td[i-1] - h) < 0.1*h &&
                qAbs(log_td[i+2] - log_td[i+1] - h) < 0.1*h &&
                qAbs(log_td[i-1] - log_td[i-2] - h) < 0.1*h) {
                // 等间距情况
                derivative = (log_pd[i-2] - 8*log_pd[i-1] + 8*log_pd[i+1] - log_pd[i+2]) / (12*h);
            } else {
                // 非等间距，使用三点公式
                double h1 = log_td[i] - log_td[i-1];
                double h2 = log_td[i+1] - log_td[i];
                double weight1 = h2 / (h1 + h2);
                double weight2 = h1 / (h1 + h2);
                derivative = weight1 * (log_pd[i] - log_pd[i-1])/h1 +
                             weight2 * (log_pd[i+1] - log_pd[i])/h2;
            }
        } else if (i == log_td.size() - 2 && log_td.size() >= 3) {
            // 三点公式
            int n = log_td.size();
            double h0 = log_td[n-2] - log_td[n-3];
            double h1 = log_td[n-1] - log_td[n-2];
            double a = h1 / (h0 * (h0 + h1));
            double b = -(h0 + h1) / (h0 * h1);
            double c = (2*h0 + h1) / (h1 * (h0 + h1));
            derivative = a * log_pd[n-3] + b * log_pd[n-2] + c * log_pd[n-1];
        } else if (i == log_td.size() - 1) {
            // 后向差分
            if (log_td.size() >= 2) {
                int n = log_td.size();
                derivative = (log_pd[n-1] - log_pd[n-2]) / (log_td[n-1] - log_td[n-2]);
            }
        }

        dlogP_dlogt.append(derivative);
    }

    // Step 4: 对导数进行额外的平滑
    QVector<double> smoothedDerivative = dlogP_dlogt;

    for (int pass = 0; pass < 2; ++pass) {
        QVector<double> temp = smoothedDerivative;

        for (int i = 2; i < smoothedDerivative.size() - 2; ++i) {
            double sum = temp[i-2] * 0.05 +
                         temp[i-1] * 0.25 +
                         temp[i] * 0.4 +
                         temp[i+1] * 0.25 +
                         temp[i+2] * 0.05;
            smoothedDerivative[i] = sum;
        }
    }

    // Step 5: 转换回压力导数
    dpd.reserve(td.size() - 1);
    td_dpd.reserve(td.size() - 1);

    for (int i = 1; i < td.size(); ++i) {
        if (i-1 < smoothedDerivative.size()) {
            double dpd_val = smoothedPd[i] * smoothedDerivative[i-1];

            // 异常值过滤和修正
            if (dpd_val > 0 && std::isfinite(dpd_val)) {
                if (dpd.size() > 0) {
                    double ratio = dpd_val / dpd.last();
                    if (ratio > 5.0 || ratio < 0.2) {
                        if (dpd.size() > 1) {
                            double log_dpd_prev = log(dpd.last());
                            double log_dpd_prev2 = log(dpd[dpd.size()-2]);
                            double log_td_prev = log(td_dpd.last());
                            double log_td_prev2 = log(td_dpd[td_dpd.size()-2]);
                            double log_td_curr = log(td[i]);

                            double slope = (log_dpd_prev - log_dpd_prev2) / (log_td_prev - log_td_prev2);
                            double log_dpd_pred = log_dpd_prev + slope * (log_td_curr - log_td_prev);
                            dpd_val = exp(log_dpd_pred);
                        } else {
                            dpd_val = dpd.last() * 1.1;
                        }
                    }
                }

                dpd.append(dpd_val);
                td_dpd.append(td[i]);
            }
        }
    }
}

void ModelWidget2::smoothData(CalculationResultFD& result, const CalculationParametersFD& params)
{
    if (result.tD.size() < 3) {
        return;
    }

    // 在对数空间生成更密集的时间序列
    result.tD_smooth.clear();
    result.tD_smooth.reserve(params.fittingPoints);

    double logMin = log10(result.tD.first());
    double logMax = log10(result.tD.last());
    double logStep = (logMax - logMin) / (params.fittingPoints - 1);

    for (int i = 0; i < params.fittingPoints; ++i) {
        double logT = logMin + i * logStep;
        result.tD_smooth.append(pow(10.0, logT));
    }

    // 对压力敏感处理后的压力进行样条插值
    QVector<double> logTD, logPD;
    for (int i = 0; i < result.tD.size(); ++i) {
        if (result.tD[i] > 0 && result.pd1[i] > 0) {
            logTD.append(log10(result.tD[i]));
            logPD.append(log10(result.pd1[i]));
        }
    }

    SimpleSplineFD pressureSpline;
    pressureSpline.setData(logTD, logPD);

    if (!pressureSpline.isValid()) {
        result.hasSmoothedData = false;
        return;
    }

    result.pd_smooth.clear();
    result.pd_smooth.reserve(params.fittingPoints);

    for (double tD : result.tD_smooth) {
        double logT = log10(tD);
        double logP = pressureSpline.interpolate(logT);
        result.pd_smooth.append(pow(10.0, logP));
    }

    // 使用改进的方法计算平滑后数据的压力导数
    computePressureDerivative(result.tD_smooth, result.pd_smooth, params.cD,
                              result.dpd_smooth, result.td_dpd_smooth);

    result.hasSmoothedData = true;
}




QVector<double> ModelWidget2::movingAverage(const QVector<double>& data, int windowSize)
{
    if (data.size() < windowSize || windowSize < 1) {
        return data;
    }

    QVector<double> smoothed = data;
    int halfWindow = windowSize / 2;

    for (int i = halfWindow; i < data.size() - halfWindow; ++i) {
        double sum = 0.0;
        int count = 0;

        for (int j = -halfWindow; j <= halfWindow; ++j) {
            sum += data[i + j];
            count++;
        }

        if (count > 0) {
            smoothed[i] = sum / count;
        }
    }

    return smoothed;
}

// 数学计算函数实现（保持不变）...
// [以下所有数学函数保持原样，包括：
//  solveLinearSystem, stefestCoefficient, factorial, flaplace,
//  e_function, f_function, integralBesselK0, besselK0, gaussQuadrature]

QVector<double> ModelWidget2::solveLinearSystem(const QVector<QVector<double>>& A, const QVector<double>& b)
{
    int n = b.size();
    QVector<QVector<double>> Ab(n);

    // 构建增广矩阵
    for (int i = 0; i < n; ++i) {
        Ab[i] = A[i];
        Ab[i].append(b[i]);
    }

    // 高斯消元法
    for (int k = 0; k < n - 1; ++k) {
        // 部分主元选择
        int maxRow = k;
        double maxVal = qAbs(Ab[k][k]);
        for (int i = k + 1; i < n; ++i) {
            if (qAbs(Ab[i][k]) > maxVal) {
                maxVal = qAbs(Ab[i][k]);
                maxRow = i;
            }
        }

        if (maxRow != k) {
            Ab[k].swap(Ab[maxRow]);
        }

        if (qAbs(Ab[k][k]) < 1e-12) {
            continue;
        }

        // 消元
        for (int i = k + 1; i < n; ++i) {
            double factor = Ab[i][k] / Ab[k][k];
            for (int j = k; j <= n; ++j) {
                Ab[i][j] -= factor * Ab[k][j];
            }
        }
    }

    // 回代
    QVector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        if (qAbs(Ab[i][i]) < 1e-12) {
            x[i] = 0.0;
            continue;
        }

        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += Ab[i][j] * x[j];
        }
        x[i] = (Ab[i][n] - sum) / Ab[i][i];
    }

    return x;
}

double ModelWidget2::stefestCoefficient(int i, int N)
{
    // 与MATLAB完全一致的Stefest系数计算
    double sum = 0.0;
    int start = (i + 1) / 2;
    int end = qMin(i, N / 2);

    for (int k = start; k <= end; ++k) {
        double numerator = pow(k, N/2.0) * factorial(2*k);
        double denominator = factorial(N/2-k) * factorial(k) * factorial(k-1) *
                             factorial(i-k) * factorial(2*k-i);

        if (denominator > 0) {
            sum += numerator / denominator;
        }
    }

    int sign = ((N/2 + i) % 2 == 0) ? 1 : -1;
    return sign * sum;
}

double ModelWidget2::factorial(int n)
{
    if (n < 0) return 0.0;
    if (n <= 1) return 1.0;

    static QVector<double> cache = {1.0, 1.0};

    if (n < cache.size()) {
        return cache[n];
    }

    for (int i = cache.size(); i <= n; ++i) {
        double val = cache[i-1] * i;
        if (std::isinf(val)) {
            return 1e308;
        }
        cache.append(val);
    }

    return cache[n];
}

double ModelWidget2::flaplace(double z, const CalculationParametersFD& params)
{
    int mf = params.mf;
    int nf = params.nf;
    double omega = params.omega;
    double lambda = params.lambda;
    double Xf = params.Xf;
    double yy = params.yy;
    double y = params.y;
    double CFD = params.CFD;
    double S = params.S;
    double cD = params.cD;

    // 计算deltaL
    double deltaL = Xf / (nf * y);

    int matrixSize = mf * nf;

    // 构建矩阵
    QVector<QVector<double>> E1(matrixSize, QVector<double>(matrixSize, 0.0));
    QVector<double> E2(matrixSize, 0.0);
    QVector<double> E3(matrixSize, -1.0);

    // 填充E1矩阵 - 与MATLAB FD_1完全一致
    for (int i = 1; i <= mf; ++i) {
        for (int j = 1; j <= nf; ++j) {
            for (int k = 1; k <= mf; ++k) {
                for (int v = 1; v <= nf; ++v) {
                    int row_index = (i-1) * nf + j - 1;
                    int col_index = (k-1) * nf + v - 1;

                    if (row_index >= 0 && row_index < matrixSize &&
                        col_index >= 0 && col_index < matrixSize) {
                        E1[row_index][col_index] = e_function(z, i, j, k, v, mf, nf,
                                                              omega, lambda, Xf, yy, y, CFD, deltaL);
                    }
                }
            }
        }
    }

    // 填充E2向量 - 修正：每个裂缝段对应一个值
    for (int i = 1; i <= mf; ++i) {
        for (int j = 1; j <= nf; ++j) {
            int index = (i-1) * nf + j - 1;
            if (index >= 0 && index < matrixSize) {
                E2[index] = (2 * M_PI / CFD) * f_function(j, nf, Xf, y);
            }
        }
    }

    // 构建G1矩阵
    QVector<QVector<double>> G1(mf, QVector<double>(matrixSize + mf, 0.0));
    for (int i = 1; i <= mf; ++i) {
        int start_col = (i-1) * nf;
        int end_col = i * nf - 1;
        for (int j = start_col; j <= end_col; ++j) {
            if (j < matrixSize) {
                G1[i-1][j] = deltaL;
            }
        }
        G1[i-1][matrixSize + i - 1] = -1.0;
    }

    // 构建G2向量
    QVector<double> G2(mf, 0.0);

    // 构建G3向量 - 修正：最后mf个元素都是1
    QVector<double> G3(matrixSize + mf, 0.0);
    for (int i = matrixSize; i < matrixSize + mf; ++i) {
        G3[i] = 1.0;
    }

    double G4 = 0.0;

    // 组装完整矩阵
    int totalSize = matrixSize + mf + 1;
    QVector<QVector<double>> I(totalSize, QVector<double>(totalSize, 0.0));

    // E部分
    for (int i = 0; i < matrixSize; ++i) {
        for (int j = 0; j < matrixSize; ++j) {
            I[i][j] = E1[i][j];
        }
        I[i][matrixSize] = E2[i];
        I[i][totalSize-1] = E3[i];
    }

    // G部分
    for (int i = 0; i < mf; ++i) {
        for (int j = 0; j < matrixSize + mf; ++j) {
            I[matrixSize + i][j] = G1[i][j];
        }
        I[matrixSize + i][totalSize-1] = G2[i];
    }

    for (int j = 0; j < matrixSize + mf; ++j) {
        I[totalSize-1][j] = G3[j];
    }
    I[totalSize-1][totalSize-1] = G4;

    // 构建F向量
    QVector<double> F(totalSize, 0.0);
    F[totalSize-1] = 1.0 / z;

    // 求解线性系统
    QVector<double> result_matrix = solveLinearSystem(I, F);

    if (result_matrix.size() > totalSize-1) {
        double pd1 = result_matrix[totalSize-1];
        // 计算最终的压力响应
        double pwd = (z * pd1 + S) / (z + z * z * cD * (z * pd1 + S));
        return pwd;
    }

    return 0.0;
}

double ModelWidget2::e_function(double z, int i, int j, int k, int v, int mf, int nf,
                                double omega, double lambda, double Xf, double yy, double y,
                                double CFD, double deltaL)
{
    // 与MATLAB FD_1完全一致的计算
    double fz = (z * (omega * z * (1 - omega) + lambda)) / (lambda + (1 - omega) * z);

    double xij = (j - nf - 1) / (double)nf * Xf;
    double xDij = xij / y;
    double xij1 = (j - nf) / (double)nf * Xf;
    double xDij1 = xij1 / y;

    double Xkv = (2*v - 2*nf - 1) / (double)(2*nf) * Xf;
    double XDkv = Xkv / y;

    double yij = yy + (y - 2*yy) / (mf - 1) * (i - 1);
    double yDij = yij / y;

    double Ykv = yy + (y - 2*yy) / (mf - 1) * (k - 1);
    double YDkv = Ykv / y;

    double integral = integralBesselK0(XDkv, YDkv, yDij, fz, xDij, xDij1);

    // 添加有限导流项 - 与MATLAB FD_1一致
    if (i == k) {
        if (j == v) {
            integral -= M_PI / (CFD * 8) * deltaL * deltaL;
        } else {
            integral -= M_PI / CFD * qAbs(j - v) * deltaL * deltaL;
        }
    }

    return integral;
}

double ModelWidget2::f_function(int j, int nf, double Xf, double y)
{
    // 与MATLAB完全一致的计算 - 返回裂缝段长度
    double xij = (j - nf - 1) / (double)nf * Xf;
    double xDij = xij / y;
    double xij1 = (j - nf) / (double)nf * Xf;
    double xDij1 = xij1 / y;

    return xDij1 - xDij;
}

double ModelWidget2::integralBesselK0(double XDkv, double YDkv, double yDij,
                                      double fz, double xDij, double xDij1)
{
    return gaussQuadrature(XDkv, YDkv, yDij, fz, xDij, xDij1);
}

double ModelWidget2::besselK0(double x)
{
    // 快速精确的贝塞尔函数K0实现
    if (x <= 0) return 1e10;

    if (x <= 2.0) {
        double t = x / 2.0;
        double t2 = t * t;

        double I0 = 1.0 + t2 * (1.0 + t2 * (0.25 + t2 * (1.0/36.0 +
                                                         t2 * (1.0/576.0 + t2 * 1.0/14400.0))));

        double lnTerm = -log(t) * I0;
        double series = -0.5772156649015329 + t2 * (0.4227843350984671 +
                                                    t2 * (0.2306975606171077 + t2 * 0.0348859015277786));

        return lnTerm + series;
    } else {
        double y = 2.0 / x;
        double y2 = y * y;

        double asymp = 1.0 + y2 * (-0.125 + y2 * (0.0703125 +
                                                  y2 * (-0.0732421875 + y2 * 0.1121520996094)));

        return sqrt(M_PI / (2.0 * x)) * exp(-x) * asymp;
    }
}

double ModelWidget2::gaussQuadrature(double XDkv, double YDkv, double yDij, double fz,
                                     double a, double b)
{
    // 使用8点高斯-勒让德积分
    static const double points[] = {
        -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
        0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363
    };

    static const double weights[] = {
        0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
        0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763
    };

    if (qAbs(b - a) < 1e-15) return 0.0;

    double sum = 0.0;
    double halfWidth = (b - a) / 2.0;
    double center = (a + b) / 2.0;
    double sqrtFz = sqrt(qAbs(fz));

    for (int i = 0; i < 8; ++i) {
        double x = center + halfWidth * points[i];
        double distance = sqrt((XDkv - x) * (XDkv - x) + (YDkv - yDij) * (YDkv - yDij));
        double arg = distance * sqrtFz;

        double value = besselK0(arg);
        sum += weights[i] * value;
    }

    return sum * halfWidth;
}
