#include "plottingwidget.h"
#include <QPaintEvent>
#include <QPainter>
#include <QApplication>
#include <QClipboard>
#include <QPrintPreviewDialog>
#include <QPrinter>
#include <QPixmap>
#include <QDateTime>
#include <QDebug>
#include <QtMath>
#include <QMessageBox>
#include <QDialog>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QListWidgetItem>
#include <QSvgGenerator>
#include <QTableWidget>
#include <QTableWidgetItem>
#include <QMenuBar>
#include <QToolBar>
#include <QStatusBar>
#include <QGroupBox>
#include <QGridLayout>
#include <QFileDialog>

// =======================
// PlotWindow 类实现
// =======================

PlotWindow::PlotWindow(const QString &windowTitle, QWidget *parent)
    : QMainWindow(parent)
    , m_isDragging(false)
    , m_isSelecting(false)
    , m_isPanning(false)
    , m_isDraggingLegend(false)
    , m_legendOffset(0, 0)
{
    setWindowTitle(windowTitle);
    resize(900, 700);
    setupUI();
    setupContextMenu();
    initializePlotSettings();

    setMouseTracking(true);
    m_plotWidget->setMouseTracking(true);
    m_plotWidget->installEventFilter(this);
}

PlotWindow::~PlotWindow()
{
}

bool PlotWindow::eventFilter(QObject *obj, QEvent *event)
{
    if (obj == m_plotWidget && event->type() == QEvent::Paint) {
        QPainter painter(m_plotWidget);
        painter.setRenderHint(QPainter::Antialiasing, true);
        painter.setRenderHint(QPainter::TextAntialiasing, true);
        paintPlotArea(painter);
        return true;
    }
    return QMainWindow::eventFilter(obj, event);
}

void PlotWindow::setupUI()
{
    m_centralWidget = new QWidget(this);
    setCentralWidget(m_centralWidget);

    QVBoxLayout *layout = new QVBoxLayout(m_centralWidget);

    m_plotWidget = new QWidget();
    m_plotWidget->setMinimumSize(800, 600);
    // 美化绘图区域边框
    m_plotWidget->setStyleSheet(
        "QWidget { "
        "   background-color: white; "
        "   border: 3px solid #2196F3; "
        "   border-radius: 8px; "
        "}"
        );
    layout->addWidget(m_plotWidget);

    // 美化工具栏
    QToolBar *toolbar = addToolBar("🛠️ 绘图工具");
    toolbar->setStyleSheet(
        "QToolBar { "
        "   background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1, "
        "                              stop: 0 #E3F2FD, stop: 1 #BBDEFB); "
        "   border: none; "
        "   border-bottom: 2px solid #2196F3; "
        "   spacing: 8px; "
        "   padding: 4px; "
        "}"
        "QToolButton { "
        "   background: white; "
        "   border: 2px solid #2196F3; "
        "   border-radius: 6px; "
        "   padding: 6px 12px; "
        "   color: #1976D2; "
        "   font-weight: bold; "
        "}"
        "QToolButton:hover { "
        "   background: #2196F3; "
        "   color: white; "
        "}"
        "QToolButton:pressed { "
        "   background: #1976D2; "
        "}"
        );

    toolbar->setIconSize(QSize(20, 20));
    toolbar->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);

    toolbar->addAction("💾 导出图像", this, &PlotWindow::onExportPlot);
    toolbar->addSeparator();
    toolbar->addAction("🔲 网格", this, &PlotWindow::onToggleGrid);
    toolbar->addAction("📋 图例", this, &PlotWindow::onToggleLegend);
    toolbar->addSeparator();
    toolbar->addAction("🔍➕ 放大", this, &PlotWindow::onZoomIn);
    toolbar->addAction("🔍➖ 缩小", this, &PlotWindow::onZoomOut);
    toolbar->addAction("🔄 重置", this, &PlotWindow::onResetZoom);

    // 美化状态栏
    QStatusBar *status = statusBar();
    status->setStyleSheet(
        "QStatusBar { "
        "   background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1, "
        "                              stop: 0 #E3F2FD, stop: 1 #BBDEFB); "
        "   color: #1565C0; "
        "   font-weight: bold; "
        "   border-top: 2px solid #2196F3; "
        "}"
        );
    status->showMessage("✅ 准备就绪");
}

void PlotWindow::setupContextMenu()
{
    m_contextMenu = new QMenu(this);
    m_contextMenu->setStyleSheet(
        "QMenu { "
        "   background-color: white; "
        "   border: 2px solid #2196F3; "
        "   border-radius: 6px; "
        "   padding: 8px; "
        "}"
        "QMenu::item { "
        "   padding: 8px 24px; "
        "   color: #212121; "
        "   border-radius: 4px; "
        "   margin: 2px; "
        "}"
        "QMenu::item:selected { "
        "   background-color: #2196F3; "
        "   color: white; "
        "}"
        "QMenu::separator { "
        "   height: 2px; "
        "   background: #BBDEFB; "
        "   margin: 4px 8px; "
        "}"
        );

    m_exportAction = m_contextMenu->addAction("💾 导出图像", this, &PlotWindow::onExportPlot);
    m_contextMenu->addSeparator();
    m_gridAction = m_contextMenu->addAction("🔲 显示网格", this, &PlotWindow::onToggleGrid);
    m_gridAction->setCheckable(true);
    m_gridAction->setChecked(true);
    m_legendAction = m_contextMenu->addAction("📋 显示图例", this, &PlotWindow::onToggleLegend);
    m_legendAction->setCheckable(true);
    m_legendAction->setChecked(true);
    m_contextMenu->addSeparator();
    m_zoomInAction = m_contextMenu->addAction("🔍➕ 放大", this, &PlotWindow::onZoomIn);
    m_zoomOutAction = m_contextMenu->addAction("🔍➖ 缩小", this, &PlotWindow::onZoomOut);
    m_resetZoomAction = m_contextMenu->addAction("🔄 重置缩放", this, &PlotWindow::onResetZoom);
}

void PlotWindow::initializePlotSettings()
{
    m_plotSettings.showGrid = true;
    m_plotSettings.logScaleX = false;
    m_plotSettings.logScaleY = false;
    m_plotSettings.backgroundColor = Qt::white;
    m_plotSettings.gridColor = QColor(224, 224, 224);
    m_plotSettings.textColor = Qt::black;
    m_plotSettings.lineWidth = 2;
    m_plotSettings.pointSize = 4;
    m_plotSettings.xAxisTitle = "X轴";
    m_plotSettings.yAxisTitle = "Y轴";
    m_plotSettings.plotTitle = "数据曲线";
    m_plotSettings.autoScale = true;
    m_plotSettings.xMin = 0;
    m_plotSettings.xMax = 100;
    m_plotSettings.yMin = 0;
    m_plotSettings.yMax = 100;
    m_plotSettings.showLegend = true;
    m_plotSettings.legendPosition = QPointF(0.8, 0.1);
    m_plotSettings.xAxisType = AxisType::Linear;
    m_plotSettings.yAxisType = AxisType::Linear;
}

void PlotWindow::addCurve(const CurveData &curve)
{
    m_curves.append(curve);

    if (!m_curves.isEmpty()) {
        const CurveData &firstCurve = m_curves.first();
        m_plotSettings.logScaleX = (firstCurve.xAxisType == AxisType::Logarithmic);
        m_plotSettings.logScaleY = (firstCurve.yAxisType == AxisType::Logarithmic);
        m_plotSettings.xAxisType = firstCurve.xAxisType;
        m_plotSettings.yAxisType = firstCurve.yAxisType;
    }

    calculateDataBounds();
    updatePlot();
}

void PlotWindow::updatePlot()
{
    if (m_plotWidget) {
        m_plotWidget->update();
    }
}

void PlotWindow::setAxisSettings(bool logX, bool logY, const QString &xTitle, const QString &yTitle)
{
    m_plotSettings.logScaleX = logX;
    m_plotSettings.logScaleY = logY;
    m_plotSettings.xAxisTitle = xTitle;
    m_plotSettings.yAxisTitle = yTitle;
    m_plotSettings.xAxisType = logX ? AxisType::Logarithmic : AxisType::Linear;
    m_plotSettings.yAxisType = logY ? AxisType::Logarithmic : AxisType::Linear;
    updatePlot();
}

void PlotWindow::setPlotTitle(const QString &title)
{
    m_plotSettings.plotTitle = title;
    updatePlot();
}

void PlotWindow::paintEvent(QPaintEvent *event)
{
    QMainWindow::paintEvent(event);
}

void PlotWindow::paintPlotArea(QPainter &painter)
{
    QRect widgetRect = m_plotWidget->rect();
    m_plotArea = QRect(80, 50, widgetRect.width() - 160, widgetRect.height() - 100);

    drawBackground(painter);

    if (m_plotSettings.showGrid) {
        drawGrid(painter);
    }

    drawAxes(painter);

    if (!m_curves.isEmpty()) {
        drawCurves(painter);
    }

    drawMarkers(painter);
    drawAnnotations(painter);

    if (m_isSelecting) {
        drawSelection(painter);
    }

    if (m_plotSettings.showLegend && !m_curves.isEmpty()) {
        drawLegend(painter);
    }
}

// 绘图函数实现
void PlotWindow::drawBackground(QPainter &painter)
{
    painter.fillRect(m_plotArea, m_plotSettings.backgroundColor);
    // 移除外围边框，只保留绘图区域内的边框
    painter.setPen(QPen(Qt::black, 1));
    painter.drawRect(m_plotArea);
}

void PlotWindow::drawGrid(QPainter &painter)
{
    painter.setPen(QPen(m_plotSettings.gridColor, 1, Qt::DotLine));

    QVector<double> xLabels = generateOptimizedAxisLabels(m_plotSettings.xMin, m_plotSettings.xMax, m_plotSettings.xAxisType);
    QVector<double> yLabels = generateOptimizedAxisLabels(m_plotSettings.yMin, m_plotSettings.yMax, m_plotSettings.yAxisType);

    // 绘制X轴网格线 - 使用正确的对数转换
    for (double value : xLabels) {
        double x;
        if (m_plotSettings.xAxisType == AxisType::Logarithmic && m_plotSettings.xMin > 0) {
            double normalizedX = (log10(value) - log10(m_plotSettings.xMin)) /
                                 (log10(m_plotSettings.xMax) - log10(m_plotSettings.xMin));
            x = m_plotArea.left() + normalizedX * m_plotArea.width();
        } else {
            x = m_plotArea.left() + (value - m_plotSettings.xMin) /
                                        (m_plotSettings.xMax - m_plotSettings.xMin) * m_plotArea.width();
        }

        if (x >= m_plotArea.left() && x <= m_plotArea.right()) {
            painter.drawLine(x, m_plotArea.top(), x, m_plotArea.bottom());
        }
    }

    // 绘制Y轴网格线 - 使用正确的对数转换
    for (double value : yLabels) {
        double y;
        if (m_plotSettings.yAxisType == AxisType::Logarithmic && m_plotSettings.yMin > 0) {
            double normalizedY = (log10(value) - log10(m_plotSettings.yMin)) /
                                 (log10(m_plotSettings.yMax) - log10(m_plotSettings.yMin));
            y = m_plotArea.bottom() - normalizedY * m_plotArea.height();
        } else {
            y = m_plotArea.bottom() - (value - m_plotSettings.yMin) /
                                          (m_plotSettings.yMax - m_plotSettings.yMin) * m_plotArea.height();
        }

        if (y >= m_plotArea.top() && y <= m_plotArea.bottom()) {
            painter.drawLine(m_plotArea.left(), y, m_plotArea.right(), y);
        }
    }
}

void PlotWindow::drawAxes(QPainter &painter)
{
    painter.setPen(QPen(Qt::black, 2));
    painter.setFont(QFont("Arial", 9));

    QVector<double> xLabels = generateOptimizedAxisLabels(m_plotSettings.xMin, m_plotSettings.xMax, m_plotSettings.xAxisType);
    QVector<double> yLabels = generateOptimizedAxisLabels(m_plotSettings.yMin, m_plotSettings.yMax, m_plotSettings.yAxisType);

    // 绘制X轴标签 - 参考Saphir软件风格
    for (double value : xLabels) {
        double x;
        if (m_plotSettings.xAxisType == AxisType::Logarithmic && m_plotSettings.xMin > 0) {
            double normalizedX = (log10(value) - log10(m_plotSettings.xMin)) /
                                 (log10(m_plotSettings.xMax) - log10(m_plotSettings.xMin));
            x = m_plotArea.left() + normalizedX * m_plotArea.width();
        } else {
            x = m_plotArea.left() + (value - m_plotSettings.xMin) /
                                        (m_plotSettings.xMax - m_plotSettings.xMin) * m_plotArea.width();
        }

        if (x >= m_plotArea.left() && x <= m_plotArea.right()) {
            painter.drawLine(x, m_plotArea.bottom(), x, m_plotArea.bottom() - 8);
            QString label = formatAxisLabel(value, m_plotSettings.logScaleX);
            QRect textRect(x - 30, m_plotArea.bottom() + 5, 60, 15);
            painter.drawText(textRect, Qt::AlignCenter, label);
        }
    }

    // 绘制Y轴标签 - 参考Saphir软件风格
    for (double value : yLabels) {
        double y;
        if (m_plotSettings.yAxisType == AxisType::Logarithmic && m_plotSettings.yMin > 0) {
            double normalizedY = (log10(value) - log10(m_plotSettings.yMin)) /
                                 (log10(m_plotSettings.yMax) - log10(m_plotSettings.yMin));
            y = m_plotArea.bottom() - normalizedY * m_plotArea.height();
        } else {
            y = m_plotArea.bottom() - (value - m_plotSettings.yMin) /
                                          (m_plotSettings.yMax - m_plotSettings.yMin) * m_plotArea.height();
        }

        if (y >= m_plotArea.top() && y <= m_plotArea.bottom()) {
            painter.drawLine(m_plotArea.left(), y, m_plotArea.left() + 8, y);
            QString label = formatAxisLabel(value, m_plotSettings.logScaleY);
            QRect textRect(m_plotArea.left() - 75, y - 8, 70, 16);
            painter.drawText(textRect, Qt::AlignRight | Qt::AlignVCenter, label);
        }
    }

    painter.setFont(QFont("Arial", 10, QFont::Bold));

    QRect xTitleRect(m_plotArea.left(), m_plotArea.bottom() + 30, m_plotArea.width(), 15);
    painter.drawText(xTitleRect, Qt::AlignCenter, m_plotSettings.xAxisTitle);

    painter.save();
    painter.translate(15, m_plotArea.center().y());
    painter.rotate(-90);
    painter.drawText(-80, -3, 160, 15, Qt::AlignCenter, m_plotSettings.yAxisTitle);
    painter.restore();

    painter.setFont(QFont("Arial", 12, QFont::Bold));
    QRect titleRect(m_plotArea.left(), 5, m_plotArea.width(), 35);
    painter.drawText(titleRect, Qt::AlignCenter, m_plotSettings.plotTitle);
}

void PlotWindow::drawCurves(QPainter &painter)
{
    for (const CurveData &curve : m_curves) {
        if (curve.visible) {
            drawCurve(painter, curve);
        }
    }
}

void PlotWindow::drawCurve(QPainter &painter, const CurveData &curve)
{
    if (curve.drawType == CurveType::Step) {
        drawStepCurve(painter, curve);
        return;
    }

    if (curve.xData.isEmpty() || curve.yData.isEmpty()) {
        return;
    }

    Qt::PenStyle penStyle = lineStyleToQt(curve.lineStyle);
    painter.setPen(QPen(curve.color, curve.lineWidth, penStyle));

    QVector<QPointF> points;
    int maxPoints = qMin(curve.xData.size(), curve.yData.size());

    for (int i = 0; i < maxPoints; ++i) {
        if (!isValidDataPoint(curve.xData[i], curve.yData[i])) {
            continue;
        }

        QPointF dataPoint(curve.xData[i], curve.yData[i]);
        QPointF pixelPoint = dataToPixel(dataPoint);

        if (pixelPoint.x() >= m_plotArea.left() - 50 && pixelPoint.x() <= m_plotArea.right() + 50 &&
            pixelPoint.y() >= m_plotArea.top() - 50 && pixelPoint.y() <= m_plotArea.bottom() + 50) {
            points.append(pixelPoint);
        }
    }

    if (points.size() > 1) {
        painter.setClipRect(m_plotArea.adjusted(-5, -5, 5, 5));
        for (int i = 0; i < points.size() - 1; ++i) {
            painter.drawLine(points[i], points[i + 1]);
        }
        painter.setClipping(false);
    }

    painter.setBrush(curve.color);
    for (const QPointF &point : points) {
        if (m_plotArea.contains(point.toPoint())) {
            painter.drawEllipse(point, curve.pointSize/2, curve.pointSize/2);
        }
    }
}

void PlotWindow::drawStepCurve(QPainter &painter, const CurveData &curve)
{
    if (curve.xData.isEmpty() || curve.yData.isEmpty()) {
        return;
    }

    Qt::PenStyle penStyle = lineStyleToQt(curve.lineStyle);
    painter.setPen(QPen(curve.color, curve.lineWidth, penStyle));

    painter.setClipRect(m_plotArea.adjusted(-5, -5, 5, 5));

    int dataSize = qMin(curve.xData.size(), curve.yData.size());

    for (int i = 0; i < dataSize - 1; i += 2) {
        if (i + 1 >= dataSize) break;

        if (!isValidDataPoint(curve.xData[i], curve.yData[i]) ||
            !isValidDataPoint(curve.xData[i+1], curve.yData[i+1])) {
            continue;
        }

        QPointF segmentStart = dataToPixel(QPointF(curve.xData[i], curve.yData[i]));
        QPointF segmentEnd = dataToPixel(QPointF(curve.xData[i+1], curve.yData[i+1]));

        painter.drawLine(segmentStart, segmentEnd);

        if (i + 2 < dataSize && i + 3 < dataSize) {
            if (isValidDataPoint(curve.xData[i+2], curve.yData[i+2])) {
                QPointF nextSegmentStart = dataToPixel(QPointF(curve.xData[i+2], curve.yData[i+2]));
                QPointF verticalStart = segmentEnd;
                QPointF verticalEnd = QPointF(segmentEnd.x(), nextSegmentStart.y());
                painter.drawLine(verticalStart, verticalEnd);
            }
        }
    }

    painter.setClipping(false);

    painter.setBrush(curve.color);
    for (int i = 0; i < dataSize; i += 2) {
        if (i < dataSize && isValidDataPoint(curve.xData[i], curve.yData[i])) {
            QPointF pixelPoint = dataToPixel(QPointF(curve.xData[i], curve.yData[i]));
            if (m_plotArea.contains(pixelPoint.toPoint())) {
                painter.drawEllipse(pixelPoint, curve.pointSize/2, curve.pointSize/2);
            }
        }
    }
}

void PlotWindow::drawLegend(QPainter &painter)
{
    if (m_curves.isEmpty()) return;

    int visibleCurveCount = 0;
    for (const CurveData &curve : m_curves) {
        if (curve.visible) {
            visibleCurveCount++;
        }
    }

    if (visibleCurveCount == 0) return;

    painter.setFont(QFont("Arial", 9));
    QFontMetrics fm(painter.font());

    int lineHeight = fm.height() + 4;
    int legendWidth = 0;
    int legendHeight = visibleCurveCount * lineHeight + 16;

    for (const CurveData &curve : m_curves) {
        if (curve.visible) {
            int textWidth = fm.horizontalAdvance(curve.name) + 40;
            legendWidth = qMax(legendWidth, textWidth);
        }
    }
    legendWidth += 20;

    int legendX = m_plotArea.right() - legendWidth - 10 + m_legendOffset.x();
    int legendY = m_plotArea.top() + 10 + m_legendOffset.y();

    m_legendArea = QRect(legendX, legendY, legendWidth, legendHeight);

    painter.setPen(QPen(QColor(150, 150, 150, 180), 1));
    painter.drawRect(m_legendArea);

    painter.setPen(QPen(Qt::black, 1));
    painter.setFont(QFont("Arial", 9, QFont::Bold));
    painter.drawText(legendX + 8, legendY + 15, "图例");

    painter.setFont(QFont("Arial", 8));
    int currentY = legendY + 25;

    for (const CurveData &curve : m_curves) {
        if (!curve.visible) continue;

        Qt::PenStyle penStyle = lineStyleToQt(curve.lineStyle);
        painter.setPen(QPen(curve.color, qMax(1, curve.lineWidth - 1), penStyle));
        painter.drawLine(legendX + 8, currentY + lineHeight/2 - 2,
                         legendX + 28, currentY + lineHeight/2 - 2);

        painter.setBrush(curve.color);
        painter.setPen(QPen(curve.color, 1));
        painter.drawEllipse(legendX + 18 - 2, currentY + lineHeight/2 - 4, 4, 4);

        painter.setPen(QPen(Qt::black, 1));
        painter.drawText(legendX + 35, currentY + lineHeight/2 + 3, curve.name);

        currentY += lineHeight;
    }
}

void PlotWindow::drawMarkers(QPainter &painter)
{
    painter.setPen(QPen(Qt::red, 2));
    painter.setBrush(Qt::red);

    for (const QPointF &marker : m_markers) {
        QPointF pixelPoint = dataToPixel(marker);
        painter.drawEllipse(pixelPoint, 6, 6);
        painter.drawLine(pixelPoint.x() - 10, pixelPoint.y(), pixelPoint.x() + 10, pixelPoint.y());
        painter.drawLine(pixelPoint.x(), pixelPoint.y() - 10, pixelPoint.x(), pixelPoint.y() + 10);
    }
}

void PlotWindow::drawAnnotations(QPainter &painter)
{
    painter.setPen(QPen(Qt::darkBlue, 1));
    painter.setFont(QFont("Arial", 8));

    for (const auto &annotation : m_annotations) {
        QPointF pixelPoint = dataToPixel(annotation.first);
        QString text = annotation.second;

        QFontMetrics fm(painter.font());
        QRect textRect = fm.boundingRect(text);
        textRect.moveCenter(pixelPoint.toPoint());
        textRect.adjust(-3, -1, 3, 1);

        painter.fillRect(textRect, QColor(255, 255, 255, 200));
        painter.drawRect(textRect);
        painter.drawText(textRect, Qt::AlignCenter, text);
    }
}

void PlotWindow::drawSelection(QPainter &painter)
{
    painter.setPen(QPen(Qt::blue, 1, Qt::DashLine));
    painter.setBrush(QColor(0, 0, 255, 30));
    painter.drawRect(m_selectionRect);
}

// 坐标转换函数 - 修正对数坐标转换
QPointF PlotWindow::dataToPixel(const QPointF &dataPoint)
{
    double x, y;

    // X轴坐标转换
    if (m_plotSettings.xAxisType == AxisType::Logarithmic && dataPoint.x() > 0 && m_plotSettings.xMin > 0) {
        double normalizedX = (log10(dataPoint.x()) - log10(m_plotSettings.xMin)) /
                             (log10(m_plotSettings.xMax) - log10(m_plotSettings.xMin));
        x = m_plotArea.left() + normalizedX * m_plotArea.width();
    } else {
        x = m_plotArea.left() + (dataPoint.x() - m_plotSettings.xMin) /
                                    (m_plotSettings.xMax - m_plotSettings.xMin) * m_plotArea.width();
    }

    // Y轴坐标转换
    if (m_plotSettings.yAxisType == AxisType::Logarithmic && dataPoint.y() > 0 && m_plotSettings.yMin > 0) {
        double normalizedY = (log10(dataPoint.y()) - log10(m_plotSettings.yMin)) /
                             (log10(m_plotSettings.yMax) - log10(m_plotSettings.yMin));
        y = m_plotArea.bottom() - normalizedY * m_plotArea.height();
    } else {
        y = m_plotArea.bottom() - (dataPoint.y() - m_plotSettings.yMin) /
                                      (m_plotSettings.yMax - m_plotSettings.yMin) * m_plotArea.height();
    }

    return QPointF(x, y);
}

QPointF PlotWindow::pixelToData(const QPointF &pixelPoint)
{
    double x, y;

    // X轴反向转换
    if (m_plotSettings.xAxisType == AxisType::Logarithmic && m_plotSettings.xMin > 0 && m_plotSettings.xMax > 0) {
        double normalizedX = (pixelPoint.x() - m_plotArea.left()) / m_plotArea.width();
        x = pow(10, log10(m_plotSettings.xMin) + normalizedX * (log10(m_plotSettings.xMax) - log10(m_plotSettings.xMin)));
    } else {
        x = m_plotSettings.xMin + (pixelPoint.x() - m_plotArea.left()) / m_plotArea.width() *
                                      (m_plotSettings.xMax - m_plotSettings.xMin);
    }

    // Y轴反向转换
    if (m_plotSettings.yAxisType == AxisType::Logarithmic && m_plotSettings.yMin > 0 && m_plotSettings.yMax > 0) {
        double normalizedY = (m_plotArea.bottom() - pixelPoint.y()) / m_plotArea.height();
        y = pow(10, log10(m_plotSettings.yMin) + normalizedY * (log10(m_plotSettings.yMax) - log10(m_plotSettings.yMin)));
    } else {
        y = m_plotSettings.yMin + (m_plotArea.bottom() - pixelPoint.y()) / m_plotArea.height() *
                                      (m_plotSettings.yMax - m_plotSettings.yMin);
    }

    return QPointF(x, y);
}

// 数据边界计算 - 增强的自适应功能，特别针对产量数据
void PlotWindow::calculateDataBounds()
{
    if (m_curves.isEmpty()) return;

    double minX = 1e10, maxX = -1e10;
    double minY = 1e10, maxY = -1e10;
    bool hasValidData = false;

    for (const CurveData &curve : m_curves) {
        if (!curve.visible || curve.xData.isEmpty() || curve.yData.isEmpty()) {
            continue;
        }

        for (int i = 0; i < qMin(curve.xData.size(), curve.yData.size()); ++i) {
            double x = curve.xData[i];
            double y = curve.yData[i];

            if (!qIsFinite(x) || !qIsFinite(y)) continue;

            // 对数坐标系需要正值
            if (m_plotSettings.xAxisType == AxisType::Logarithmic && x <= 0) continue;
            if (m_plotSettings.yAxisType == AxisType::Logarithmic && y <= 0) continue;

            minX = qMin(minX, x);
            maxX = qMax(maxX, x);
            minY = qMin(minY, y);
            maxY = qMax(maxY, y);
            hasValidData = true;
        }
    }

    if (hasValidData && minX < maxX && minY < maxY) {
        // 计算X轴范围
        QPair<double, double> xRange = calculateOptimalRange(minX, maxX, m_plotSettings.xAxisType == AxisType::Logarithmic);
        m_plotSettings.xMin = xRange.first;
        m_plotSettings.xMax = xRange.second;

        // 计算Y轴范围 - 产量数据最小值固定为0
        QPair<double, double> yRange;
        bool isProductionData = m_plotSettings.yAxisTitle.contains("产量") || m_plotSettings.yAxisTitle.contains("mÂ³");

        if (isProductionData && m_plotSettings.yAxisType == AxisType::Linear) {
            // 产量数据线性坐标，最小值固定为0
            yRange = calculateOptimalRange(0, maxY, false);
        } else {
            yRange = calculateOptimalRange(minY, maxY, m_plotSettings.yAxisType == AxisType::Logarithmic);
        }

        m_plotSettings.yMin = yRange.first;
        m_plotSettings.yMax = yRange.second;
    }
}

// 优化的范围计算函数 - 更智能的自适应
QPair<double, double> PlotWindow::calculateOptimalRange(double min, double max, bool isLog)
{
    if (max <= min) {
        return QPair<double, double>(min, max);
    }

    if (isLog) {
        // 对数坐标系：扩展到最近的10的整数次幂，参考Saphir软件
        if (min <= 0) min = 1e-10;  // 确保正值
        if (max <= 0) max = 1;

        double logMin = log10(min);
        double logMax = log10(max);

        // 向下取整到最近的10的次幂
        double rangeMin = pow(10, floor(logMin));

        // 向上取整到最近的10的次幂
        double rangeMax = pow(10, ceil(logMax));

        // 如果数据范围太小，至少显示一个数量级
        if (rangeMax / rangeMin < 10) {
            rangeMin = rangeMin / 10;
            rangeMax = rangeMax * 10;
        }

        return QPair<double, double>(rangeMin, rangeMax);
    } else {
        // 线性坐标系：智能计算合适的范围
        double range = max - min;

        if (range == 0) {
            // 如果所有数据相同，扩展范围
            return QPair<double, double>(min - 1, max + 1);
        }

        // 计算合适的刻度间隔
        double orderOfMagnitude = pow(10, floor(log10(range)));
        double normalizedRange = range / orderOfMagnitude;

        double tickInterval;
        if (normalizedRange <= 1.5) {
            tickInterval = orderOfMagnitude * 0.2;
        } else if (normalizedRange <= 3) {
            tickInterval = orderOfMagnitude * 0.5;
        } else if (normalizedRange <= 7) {
            tickInterval = orderOfMagnitude;
        } else {
            tickInterval = orderOfMagnitude * 2;
        }

        // 计算范围，确保包含所有数据点并留有余量
        double rangeMin = floor(min / tickInterval) * tickInterval;
        double rangeMax = ceil(max / tickInterval) * tickInterval;

        // 添加5%的边距，使数据点不会贴边
        double margin = (rangeMax - rangeMin) * 0.05;
        rangeMin -= margin;
        rangeMax += margin;

        // 如果范围包含0，确保0点对齐到刻度
        if (rangeMin < 0 && rangeMax > 0) {
            rangeMin = floor(rangeMin / tickInterval) * tickInterval;
            rangeMax = ceil(rangeMax / tickInterval) * tickInterval;
        }

        return QPair<double, double>(rangeMin, rangeMax);
    }
}

// 缩放和平移
void PlotWindow::zoomAtPoint(const QPointF &point, double factor)
{
    QPointF dataPoint = pixelToData(point);

    double xRange = m_plotSettings.xMax - m_plotSettings.xMin;
    double yRange = m_plotSettings.yMax - m_plotSettings.yMin;

    double newXRange = xRange / factor;
    double newYRange = yRange / factor;

    m_plotSettings.xMin = dataPoint.x() - newXRange * (dataPoint.x() - m_plotSettings.xMin) / xRange;
    m_plotSettings.xMax = dataPoint.x() + newXRange * (m_plotSettings.xMax - dataPoint.x()) / xRange;
    m_plotSettings.yMin = dataPoint.y() - newYRange * (dataPoint.y() - m_plotSettings.yMin) / yRange;
    m_plotSettings.yMax = dataPoint.y() + newYRange * (m_plotSettings.yMax - dataPoint.y()) / yRange;

    updatePlot();
}

void PlotWindow::panView(const QPointF &delta)
{
    QPointF dataDelta = pixelToData(delta) - pixelToData(QPointF(0, 0));

    m_plotSettings.xMin -= dataDelta.x();
    m_plotSettings.xMax -= dataDelta.x();
    m_plotSettings.yMin -= dataDelta.y();
    m_plotSettings.yMax -= dataDelta.y();

    updatePlot();
}

// 辅助函数
double PlotWindow::transformToLogScale(double value, bool isLog)
{
    if (isLog && value > 0) {
        return log10(value);
    }
    return value;
}

double PlotWindow::transformFromLogScale(double value, bool isLog)
{
    if (isLog) {
        return pow(10, value);
    }
    return value;
}

// 优化的坐标轴标签生成函数，参考Saphir软件
QVector<double> PlotWindow::generateOptimizedAxisLabels(double min, double max, AxisType axisType)
{
    QVector<double> labels;

    if (max <= min) {
        return labels;
    }

    if (axisType == AxisType::Logarithmic) {
        // 对数坐标系：生成Saphir风格的对数标签
        if (min <= 0) min = 1e-10;
        if (max <= 0) max = 1;

        double logMin = log10(min);
        double logMax = log10(max);

        int startPower = qFloor(logMin);
        int endPower = qCeil(logMax);

        // 添加主要的10的整数次幂
        for (int power = startPower; power <= endPower; ++power) {
            double value = qPow(10, power);
            if (value >= min * 0.999 && value <= max * 1.001) {
                labels.append(value);
            }
        }

        // 根据范围决定是否显示中间刻度
        bool showMinorTicks = (endPower - startPower) <= 3;

        if (showMinorTicks) {
            // 添加次要刻度（2, 3, 5）
            for (int power = startPower; power < endPower; ++power) {
                for (int mult : {2, 3, 5}) {
                    double value = mult * qPow(10, power);
                    if (value > min && value < max) {
                        labels.append(value);
                    }
                }
            }
        }

        // 排序去重
        std::sort(labels.begin(), labels.end());
        labels.erase(std::unique(labels.begin(), labels.end(), [](double a, double b) {
                         return qAbs(a - b) < 1e-10;
                     }), labels.end());

    } else {
        // 线性坐标系：生成合适的刻度
        double range = max - min;
        double orderOfMagnitude = pow(10, floor(log10(range)));
        double normalizedRange = range / orderOfMagnitude;

        double tickInterval;
        int maxTicks;

        if (normalizedRange <= 1.5) {
            tickInterval = orderOfMagnitude * 0.2;
            maxTicks = 8;
        } else if (normalizedRange <= 3) {
            tickInterval = orderOfMagnitude * 0.5;
            maxTicks = 7;
        } else if (normalizedRange <= 7) {
            tickInterval = orderOfMagnitude;
            maxTicks = 8;
        } else {
            tickInterval = orderOfMagnitude * 2;
            maxTicks = 6;
        }

        // 确保第一个刻度对齐
        double firstTick = ceil(min / tickInterval) * tickInterval;

        // 生成刻度
        for (int i = 0; i < maxTicks; ++i) {
            double value = firstTick + i * tickInterval;
            if (value > max) break;
            if (value >= min) {
                labels.append(value);
            }
        }

        // 确保包含0（如果在范围内）
        if (min <= 0 && max >= 0 && !labels.contains(0)) {
            labels.append(0);
            std::sort(labels.begin(), labels.end());
        }
    }

    return labels;
}

// 格式化坐标轴标签，参考Saphir软件
QString PlotWindow::formatAxisLabel(double value, bool isLog)
{
    if (isLog) {
        // 对数坐标：对于10的整数次幂显示简洁格式
        double logValue = log10(qAbs(value));
        if (qAbs(logValue - qRound(logValue)) < 0.01) {
            int power = qRound(logValue);
            if (power == 0) return "1";
            if (power == 1) return "10";
            if (power == 2) return "100";
            if (power == 3) return "1000";
            if (power == -1) return "0.1";
            if (power == -2) return "0.01";
            return QString("10^%1").arg(power);
        } else {
            // 非10的整数次幂，显示实际数值
            if (value >= 1000) {
                return QString::number(value, 'g', 2);
            } else if (value >= 1) {
                return QString::number(value, 'f', 0);
            } else {
                return QString::number(value, 'g', 2);
            }
        }
    } else {
        // 线性坐标：智能格式化
        if (qAbs(value) >= 100000) {
            return QString::number(value, 'e', 1);
        } else if (qAbs(value) >= 1000) {
            return QString::number(value, 'f', 0);
        } else if (qAbs(value) >= 1) {
            if (qAbs(value - qRound(value)) < 0.01) {
                return QString::number(qRound(value));
            } else {
                return QString::number(value, 'f', 1);
            }
        } else if (value == 0) {
            return "0";
        } else {
            return QString::number(value, 'g', 3);
        }
    }
}

QString PlotWindow::formatScientific(double value, int decimals)
{
    if (qIsNaN(value) || !qIsFinite(value)) {
        return "N/A";
    }

    if (qAbs(value) >= 1000 || (qAbs(value) < 0.01 && value != 0)) {
        return QString::number(value, 'e', decimals);
    } else {
        return QString::number(value, 'f', decimals);
    }
}

bool PlotWindow::isValidDataPoint(double x, double y)
{
    return !qIsNaN(x) && !qIsNaN(y) && qIsFinite(x) && qIsFinite(y);
}

// 鼠标事件处理
void PlotWindow::mousePressEvent(QMouseEvent *event)
{
    QPoint plotPos = m_plotWidget->mapFromParent(event->pos());

    if (!m_plotWidget->rect().contains(plotPos)) {
        return;
    }

    m_lastMousePos = plotPos;

    if (m_legendArea.contains(plotPos)) {
        if (event->button() == Qt::LeftButton) {
            m_isDraggingLegend = true;
            return;
        }
    }

    if (!m_plotArea.contains(plotPos)) {
        return;
    }

    if (event->button() == Qt::LeftButton) {
        if (event->modifiers() & Qt::ControlModifier) {
            m_isSelecting = true;
            m_selectionStart = plotPos;
            m_selectionRect = QRect(m_selectionStart, m_selectionStart);
        } else {
            m_isDragging = true;
            m_isPanning = true;
        }
    }
}

void PlotWindow::mouseMoveEvent(QMouseEvent *event)
{
    QPoint plotPos = m_plotWidget->mapFromParent(event->pos());

    if (m_isDraggingLegend) {
        QPoint delta = plotPos - m_lastMousePos;
        m_legendOffset += delta;
        m_lastMousePos = plotPos;
        updatePlot();
        return;
    }

    m_lastMousePos = plotPos;

    if (m_isSelecting) {
        m_selectionRect = QRect(m_selectionStart, plotPos).normalized();
        updatePlot();
    } else if (m_isDragging && m_isPanning) {
        QPointF delta = plotPos - m_lastMousePos;
        panView(delta);
        m_lastMousePos = plotPos;
    }

    updatePlot();
}

void PlotWindow::mouseReleaseEvent(QMouseEvent *event)
{
    Q_UNUSED(event)

    if (m_isDraggingLegend) {
        m_isDraggingLegend = false;
        return;
    }

    if (m_isSelecting && !m_selectionRect.isEmpty()) {
        QPointF topLeft = pixelToData(m_selectionRect.topLeft());
        QPointF bottomRight = pixelToData(m_selectionRect.bottomRight());

        m_plotSettings.xMin = qMin(topLeft.x(), bottomRight.x());
        m_plotSettings.xMax = qMax(topLeft.x(), bottomRight.x());
        m_plotSettings.yMin = qMin(topLeft.y(), bottomRight.y());
        m_plotSettings.yMax = qMax(topLeft.y(), bottomRight.y());
    }

    m_isDragging = false;
    m_isSelecting = false;
    m_isPanning = false;
    updatePlot();
}

void PlotWindow::wheelEvent(QWheelEvent *event)
{
    QPoint plotPos = m_plotWidget->mapFromParent(event->position().toPoint());

    if (m_plotWidget->rect().contains(plotPos) && m_plotArea.contains(plotPos)) {
        double factor = 1.0 + event->angleDelta().y() / 1200.0;
        zoomAtPoint(plotPos, factor);
    }
}

void PlotWindow::contextMenuEvent(QContextMenuEvent *event)
{
    QPoint plotPos = m_plotWidget->mapFromParent(event->pos());

    if (m_plotWidget->rect().contains(plotPos) && m_plotArea.contains(plotPos)) {
        m_lastMousePos = plotPos;
        m_contextMenu->exec(event->globalPos());
    }
}

// 槽函数实现
void PlotWindow::onExportPlot()
{
    QString fileName = QFileDialog::getSaveFileName(this,
                                                    "导出图像",
                                                    QString("曲线图_%1.png").arg(QDateTime::currentDateTime().toString("yyyyMMdd_hhmmss")),
                                                    "PNG图像 (*.png);;JPEG图像 (*.jpg);;PDF文件 (*.pdf)");

    if (!fileName.isEmpty()) {
        QPixmap pixmap(m_plotWidget->size());
        pixmap.fill(Qt::white);
        m_plotWidget->render(&pixmap);

        if (pixmap.save(fileName)) {
            statusBar()->showMessage("图像已导出到: " + fileName, 3000);
        } else {
            QMessageBox::warning(this, "导出失败", "无法保存文件: " + fileName);
        }
    }
}

void PlotWindow::onToggleGrid()
{
    m_plotSettings.showGrid = !m_plotSettings.showGrid;
    m_gridAction->setChecked(m_plotSettings.showGrid);
    updatePlot();
}

void PlotWindow::onToggleLegend()
{
    m_plotSettings.showLegend = !m_plotSettings.showLegend;
    m_legendAction->setChecked(m_plotSettings.showLegend);
    updatePlot();
}

void PlotWindow::onResetZoom()
{
    calculateDataBounds();
    updatePlot();
}

void PlotWindow::onZoomIn()
{
    QPointF center = m_plotArea.center();
    zoomAtPoint(center, 1.25);
}

void PlotWindow::onZoomOut()
{
    QPointF center = m_plotArea.center();
    zoomAtPoint(center, 0.8);
}

// =======================
// DualPlotWindow 类实现 - 完整版本，修复所有问题
// =======================

DualPlotWindow::DualPlotWindow(const QString &windowTitle, QWidget *parent)
    : QMainWindow(parent)
    , m_pressureDragging(false)
    , m_productionDragging(false)
    , m_pressureSelecting(false)
    , m_productionSelecting(false)
    , m_pressurePanning(false)
    , m_productionPanning(false)
    , m_pressureLegendDragging(false)
    , m_productionLegendDragging(false)
    , m_lastPressureMousePos(0, 0)
    , m_lastProductionMousePos(0, 0)
    , m_pressureLegendOffset(0, 0)
    , m_productionLegendOffset(0, 0)
    , m_syncZoom(true)
    , m_syncPan(false)
{
    setWindowTitle(windowTitle);
    resize(1000, 900);
    setupUI();
    initializePlotSettings();
}

DualPlotWindow::~DualPlotWindow()
{
}

void DualPlotWindow::setupUI()
{
    m_centralWidget = new QWidget(this);
    setCentralWidget(m_centralWidget);

    QVBoxLayout *mainLayout = new QVBoxLayout(m_centralWidget);

    // 创建美化的工具栏
    QToolBar *toolbar = addToolBar("🛠️ 工具栏");
    toolbar->setStyleSheet(
        "QToolBar { "
        "   background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1, "
        "                              stop: 0 #E3F2FD, stop: 1 #BBDEFB); "
        "   border: none; "
        "   border-bottom: 2px solid #2196F3; "
        "   spacing: 10px; "
        "   padding: 6px; "
        "}"
        );

    toolbar->addAction("💾 导出图像", [this]() {
        QString fileName = QFileDialog::getSaveFileName(this,
                                                        "导出图像",
                                                        QString("压力产量分析_%1.png").arg(QDateTime::currentDateTime().toString("yyyyMMdd_hhmmss")),
                                                        "PNG图像 (*.png);;JPEG图像 (*.jpg);;PDF文件 (*.pdf)");

        if (!fileName.isEmpty()) {
            QPixmap pixmap(m_centralWidget->size());
            pixmap.fill(Qt::white);
            m_centralWidget->render(&pixmap);

            if (pixmap.save(fileName)) {
                statusBar()->showMessage("✅ 图像已导出到: " + fileName, 3000);
            } else {
                QMessageBox::warning(this, "❌ 导出失败", "无法保存文件: " + fileName);
            }
        }
    });

    toolbar->addSeparator();

    // 创建美化的开关按钮
    auto createToggleButton = [](const QString &textOn, const QString &textOff, bool checked) {
        QPushButton *btn = new QPushButton();
        btn->setCheckable(true);
        btn->setChecked(checked);
        btn->setMinimumSize(120, 36);

        QString checkedStyle =
            "QPushButton { "
            "   background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1, "
            "                              stop: 0 #66BB6A, stop: 1 #4CAF50); "
            "   color: white; "
            "   border: none; "
            "   border-radius: 6px; "
            "   font-weight: bold; "
            "   padding: 8px 16px; "
            "}";

        QString uncheckedStyle =
            "QPushButton { "
            "   background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1, "
            "                              stop: 0 #EEEEEE, stop: 1 #E0E0E0); "
            "   color: #757575; "
            "   border: 2px solid #BDBDBD; "
            "   border-radius: 6px; "
            "   font-weight: bold; "
            "   padding: 8px 16px; "
            "}";

        auto updateStyle = [btn, textOn, textOff, checkedStyle, uncheckedStyle](bool checked) {
            btn->setText(checked ? textOn : textOff);
            btn->setStyleSheet(checked ? checkedStyle : uncheckedStyle);
        };

        updateStyle(checked);
        QObject::connect(btn, &QPushButton::toggled, updateStyle);

        return btn;
    };

    // 网格开关按钮
    QPushButton *gridButton = createToggleButton("🔲 网格:开", "🔲 网格:关", true);
    connect(gridButton, &QPushButton::toggled, [this](bool checked) {
        m_pressureSettings.showGrid = checked;
        m_productionSettings.showGrid = checked;
        updatePlots();
    });
    toolbar->addWidget(gridButton);

    // 同步缩放按钮
    QPushButton *syncZoomButton = createToggleButton("🔗 同步缩放:开", "⛓️ 同步缩放:关", true);
    connect(syncZoomButton, &QPushButton::toggled, [this](bool checked) {
        m_syncZoom = checked;
    });
    toolbar->addWidget(syncZoomButton);

    // 同步平移按钮
    QPushButton *syncPanButton = createToggleButton("🔗 同步平移:开", "⛓️ 同步平移:关", false);
    connect(syncPanButton, &QPushButton::toggled, [this](bool checked) {
        m_syncPan = checked;
    });
    toolbar->addWidget(syncPanButton);

    toolbar->addSeparator();
    toolbar->addAction("🔄 重置视图", [this]() {
        calculatePressureBounds();
        calculateProductionBounds();
        synchronizeXAxis();
        updatePlots();
        statusBar()->showMessage("✅ 视图已重置", 2000);
    });

    // 创建分割器
    m_splitter = new QSplitter(Qt::Vertical);
    m_splitter->setStyleSheet(
        "QSplitter::handle { "
        "   background: qlineargradient(x1: 0, y1: 0, x2: 1, y2: 0, "
        "                              stop: 0 #BBDEFB, stop: 0.5 #2196F3, stop: 1 #BBDEFB); "
        "   height: 4px; "
        "}"
        "QSplitter::handle:hover { "
        "   background: #2196F3; "
        "}"
        );

    // 创建美化的压力绘图区域
    m_pressurePlotWidget = new QWidget();
    m_pressurePlotWidget->setMinimumHeight(350);
    m_pressurePlotWidget->setStyleSheet(
        "QWidget { "
        "   background-color: white; "
        "   border: 3px solid #2196F3; "
        "   border-radius: 8px; "
        "}"
        );
    m_pressurePlotWidget->setMouseTracking(true);

    // 创建美化的产量绘图区域
    m_productionPlotWidget = new QWidget();
    m_productionPlotWidget->setMinimumHeight(350);
    m_productionPlotWidget->setStyleSheet(
        "QWidget { "
        "   background-color: white; "
        "   border: 3px solid #4CAF50; "
        "   border-radius: 8px; "
        "}"
        );
    m_productionPlotWidget->setMouseTracking(true);

    m_splitter->addWidget(m_pressurePlotWidget);
    m_splitter->addWidget(m_productionPlotWidget);
    m_splitter->setStretchFactor(0, 1);
    m_splitter->setStretchFactor(1, 1);

    mainLayout->addWidget(m_splitter);

    // 设置绘图更新
    m_pressurePlotWidget->installEventFilter(this);
    m_productionPlotWidget->installEventFilter(this);

    // 美化状态栏
    statusBar()->setStyleSheet(
        "QStatusBar { "
        "   background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1, "
        "                              stop: 0 #E3F2FD, stop: 1 #BBDEFB); "
        "   color: #1565C0; "
        "   font-weight: bold; "
        "   border-top: 2px solid #2196F3; "
        "}"
        );
    statusBar()->showMessage("✅ 准备就绪");
}

void DualPlotWindow::initializePlotSettings()
{
    // 压力图设置
    m_pressureSettings.showGrid = true;
    m_pressureSettings.backgroundColor = Qt::white;
    m_pressureSettings.gridColor = QColor(224, 224, 224);
    m_pressureSettings.textColor = Qt::black;
    m_pressureSettings.lineWidth = 2;
    m_pressureSettings.pointSize = 4;
    m_pressureSettings.xAxisTitle = "时间 (小时)";
    m_pressureSettings.yAxisTitle = "压力 (MPa)";
    m_pressureSettings.plotTitle = "压力数据";
    m_pressureSettings.autoScale = true;
    m_pressureSettings.xMin = 0.1;
    m_pressureSettings.xMax = 1000;
    m_pressureSettings.yMin = 1;
    m_pressureSettings.yMax = 100;
    m_pressureSettings.showLegend = true;
    m_pressureSettings.xAxisType = AxisType::Linear;
    m_pressureSettings.yAxisType = AxisType::Linear;

    // 产量图设置
    m_productionSettings = m_pressureSettings;
    m_productionSettings.yAxisTitle = "产量 (m³/d)";
    m_productionSettings.plotTitle = "产量数据";
    // 产量数据最小值固定为0
    m_productionSettings.yMin = 0;
}

bool DualPlotWindow::eventFilter(QObject *obj, QEvent *event)
{
    if (event->type() == QEvent::Paint) {
        if (obj == m_pressurePlotWidget) {
            paintPressurePlot();
            return true;
        } else if (obj == m_productionPlotWidget) {
            paintProductionPlot();
            return true;
        }
    } else if (event->type() == QEvent::MouseButtonPress) {
        QMouseEvent *mouseEvent = static_cast<QMouseEvent*>(event);
        if (obj == m_pressurePlotWidget) {
            onPressureMousePress(mouseEvent);
        } else if (obj == m_productionPlotWidget) {
            onProductionMousePress(mouseEvent);
        }
    } else if (event->type() == QEvent::MouseMove) {
        QMouseEvent *mouseEvent = static_cast<QMouseEvent*>(event);
        if (obj == m_pressurePlotWidget) {
            onPressureMouseMove(mouseEvent);
        } else if (obj == m_productionPlotWidget) {
            onProductionMouseMove(mouseEvent);
        }
    } else if (event->type() == QEvent::MouseButtonRelease) {
        QMouseEvent *mouseEvent = static_cast<QMouseEvent*>(event);
        if (obj == m_pressurePlotWidget) {
            onPressureMouseRelease(mouseEvent);
        } else if (obj == m_productionPlotWidget) {
            onProductionMouseRelease(mouseEvent);
        }
    } else if (event->type() == QEvent::Wheel) {
        QWheelEvent *wheelEvent = static_cast<QWheelEvent*>(event);
        if (obj == m_pressurePlotWidget) {
            onPressureWheel(wheelEvent);
        } else if (obj == m_productionPlotWidget) {
            onProductionWheel(wheelEvent);
        }
    }
    return QMainWindow::eventFilter(obj, event);
}

void DualPlotWindow::addPressureCurve(const CurveData &curve)
{
    CurveData modifiedCurve = curve;

    // 根据用户选择设置坐标轴类型
    if (!m_pressureCurves.isEmpty()) {
        modifiedCurve.xAxisType = m_pressureCurves.first().xAxisType;
        modifiedCurve.yAxisType = m_pressureCurves.first().yAxisType;
    }

    m_pressureCurves.append(modifiedCurve);

    // 设置坐标系类型
    m_pressureSettings.xAxisType = modifiedCurve.xAxisType;
    m_pressureSettings.yAxisType = modifiedCurve.yAxisType;
    m_pressureSettings.logScaleX = (modifiedCurve.xAxisType == AxisType::Logarithmic);
    m_pressureSettings.logScaleY = (modifiedCurve.yAxisType == AxisType::Logarithmic);

    calculatePressureBounds();
    synchronizeXAxis();
    updatePlots();
}

void DualPlotWindow::addProductionCurve(const CurveData &curve)
{
    CurveData modifiedCurve = curve;

    // 根据用户选择设置坐标轴类型
    if (!m_productionCurves.isEmpty()) {
        modifiedCurve.xAxisType = m_productionCurves.first().xAxisType;
        modifiedCurve.yAxisType = m_productionCurves.first().yAxisType;
    }

    m_productionCurves.append(modifiedCurve);

    // 设置坐标系类型
    m_productionSettings.xAxisType = modifiedCurve.xAxisType;
    m_productionSettings.yAxisType = modifiedCurve.yAxisType;
    m_productionSettings.logScaleX = (modifiedCurve.xAxisType == AxisType::Logarithmic);
    m_productionSettings.logScaleY = (modifiedCurve.yAxisType == AxisType::Logarithmic);

    calculateProductionBounds();
    synchronizeXAxis();
    updatePlots();
}

void DualPlotWindow::updatePlots()
{
    if (m_pressurePlotWidget) {
        m_pressurePlotWidget->update();
    }
    if (m_productionPlotWidget) {
        m_productionPlotWidget->update();
    }
}

void DualPlotWindow::setAxisSettings(const QString &xTitle, const QString &pressureTitle, const QString &productionTitle)
{
    m_pressureSettings.xAxisTitle = xTitle;
    m_pressureSettings.yAxisTitle = pressureTitle;
    m_productionSettings.xAxisTitle = xTitle;
    m_productionSettings.yAxisTitle = productionTitle;
    updatePlots();
}

void DualPlotWindow::resizeEvent(QResizeEvent *event)
{
    QMainWindow::resizeEvent(event);
    updatePlots();
}

void DualPlotWindow::paintPressurePlot()
{
    QPainter painter(m_pressurePlotWidget);
    painter.setRenderHint(QPainter::Antialiasing, true);
    painter.setRenderHint(QPainter::TextAntialiasing, true);

    QRect widgetRect = m_pressurePlotWidget->rect();
    m_pressurePlotArea = QRect(80, 40, widgetRect.width() - 160, widgetRect.height() - 80);

    // 绘制背景 - 移除外围边框，只保留绘图区域内的边框
    painter.fillRect(m_pressurePlotArea, m_pressureSettings.backgroundColor);
    painter.setPen(QPen(Qt::black, 1));
    painter.drawRect(m_pressurePlotArea);

    // 绘制网格
    if (m_pressureSettings.showGrid) {
        drawGridOnWidget(painter, m_pressurePlotArea, m_pressureSettings);
    }

    // 绘制坐标轴
    drawAxesOnWidget(painter, m_pressurePlotArea, m_pressureSettings);

    // 绘制曲线
    for (const CurveData &curve : m_pressureCurves) {
        if (curve.visible) {
            drawCurveOnWidget(painter, curve, m_pressurePlotArea, m_pressureSettings);
        }
    }

    // 绘制图例
    if (m_pressureSettings.showLegend && !m_pressureCurves.isEmpty()) {
        drawLegendOnWidget(painter, m_pressureCurves, m_pressurePlotArea, m_pressureLegendOffset);
    }

    // 绘制选择框
    if (m_pressureSelecting) {
        painter.setPen(QPen(Qt::blue, 1, Qt::DashLine));
        painter.setBrush(QColor(0, 0, 255, 30));
        painter.drawRect(m_pressureSelectionRect);
    }
}

void DualPlotWindow::paintProductionPlot()
{
    QPainter painter(m_productionPlotWidget);
    painter.setRenderHint(QPainter::Antialiasing, true);
    painter.setRenderHint(QPainter::TextAntialiasing, true);

    QRect widgetRect = m_productionPlotWidget->rect();
    m_productionPlotArea = QRect(80, 40, widgetRect.width() - 160, widgetRect.height() - 80);

    // 绘制背景 - 移除外围边框，只保留绘图区域内的边框
    painter.fillRect(m_productionPlotArea, m_productionSettings.backgroundColor);
    painter.setPen(QPen(Qt::black, 1));
    painter.drawRect(m_productionPlotArea);

    // 绘制网格
    if (m_productionSettings.showGrid) {
        drawGridOnWidget(painter, m_productionPlotArea, m_productionSettings);
    }

    // 绘制坐标轴
    drawAxesOnWidget(painter, m_productionPlotArea, m_productionSettings);

    // 绘制曲线
    for (const CurveData &curve : m_productionCurves) {
        if (curve.visible) {
            drawCurveOnWidget(painter, curve, m_productionPlotArea, m_productionSettings);
        }
    }

    // 绘制图例
    if (m_productionSettings.showLegend && !m_productionCurves.isEmpty()) {
        drawLegendOnWidget(painter, m_productionCurves, m_productionPlotArea, m_productionLegendOffset);
    }

    // 绘制选择框
    if (m_productionSelecting) {
        painter.setPen(QPen(Qt::blue, 1, Qt::DashLine));
        painter.setBrush(QColor(0, 0, 255, 30));
        painter.drawRect(m_productionSelectionRect);
    }
}
void drawBeautifulLegend(QPainter &painter, const QVector<CurveData> &curves,
                         const QRect &plotArea, const QPoint &offset)
{
    if (curves.isEmpty()) return;

    int visibleCurveCount = 0;
    for (const CurveData &curve : curves) {
        if (curve.visible) visibleCurveCount++;
    }
    if (visibleCurveCount == 0) return;

    painter.setFont(QFont("Arial", 9));
    QFontMetrics fm(painter.font());

    int lineHeight = fm.height() + 4;
    int legendWidth = 0;
    int legendHeight = visibleCurveCount * lineHeight + 20;

    for (const CurveData &curve : curves) {
        if (curve.visible) {
            int textWidth = fm.horizontalAdvance(curve.name) + 50;
            legendWidth = qMax(legendWidth, textWidth);
        }
    }
    legendWidth += 20;

    int legendX = plotArea.right() - legendWidth - 10 + offset.x();
    int legendY = plotArea.top() + 10 + offset.y();

    legendX = qMax(plotArea.left() + 10, qMin(legendX, plotArea.right() - legendWidth - 10));
    legendY = qMax(plotArea.top() + 10, qMin(legendY, plotArea.bottom() - legendHeight - 10));

    QRect legendArea(legendX, legendY, legendWidth, legendHeight);

    painter.save();

    // 绘制阴影
    QRect shadowRect = legendArea.adjusted(3, 3, 3, 3);
    painter.setPen(Qt::NoPen);
    painter.setBrush(QColor(0, 0, 0, 50));
    painter.drawRoundedRect(shadowRect, 8, 8);

    // 绘制渐变背景
    QLinearGradient gradient(legendX, legendY, legendX, legendY + legendHeight);
    gradient.setColorAt(0, QColor(255, 255, 255, 245));
    gradient.setColorAt(1, QColor(227, 242, 253, 245));
    painter.setBrush(gradient);
    painter.setPen(QPen(QColor(33, 150, 243), 2));
    painter.drawRoundedRect(legendArea, 8, 8);

    painter.restore();

    // 绘制图例标题
    painter.setPen(QPen(QColor(25, 118, 210), 1));
    painter.setFont(QFont("Arial", 10, QFont::Bold));
    painter.drawText(legendX + 10, legendY + 18, "📋 图例");

    // 绘制分隔线
    painter.setPen(QPen(QColor(187, 222, 251), 2));
    painter.drawLine(legendX + 10, legendY + 25, legendX + legendWidth - 10, legendY + 25);

    // 绘制曲线列表
    painter.setFont(QFont("Arial", 8));
    int currentY = legendY + 32;

    for (const CurveData &curve : curves) {
        if (!curve.visible) continue;

        Qt::PenStyle penStyle = lineStyleToQt(curve.lineStyle);
        painter.setPen(QPen(curve.color, qMax(2, curve.lineWidth), penStyle));
        painter.drawLine(legendX + 10, currentY + lineHeight/2 - 2,
                         legendX + 35, currentY + lineHeight/2 - 2);

        painter.setBrush(curve.color);
        painter.setPen(QPen(curve.color, 1));
        painter.drawEllipse(legendX + 22 - 3, currentY + lineHeight/2 - 5, 6, 6);

        painter.setPen(QPen(QColor(33, 33, 33), 1));
        painter.drawText(legendX + 40, currentY + lineHeight/2 + 4, curve.name);

        currentY += lineHeight;
    }
}
// 优化的网格绘制函数 - 正确处理对数坐标间距，参考Saphir软件
void DualPlotWindow::drawGridOnWidget(QPainter &painter, const QRect &plotArea, const PlotSettings &settings)
{
    painter.setPen(QPen(settings.gridColor, 1, Qt::DotLine));

    // 生成优化的坐标轴标签
    QVector<double> xTickValues;
    QVector<double> yTickValues;

    // X轴网格线 - 参考Saphir软件的对数坐标处理
    if (settings.xAxisType == AxisType::Logarithmic && settings.xMin > 0) {
        // 对数X轴
        double logMin = log10(settings.xMin);
        double logMax = log10(settings.xMax);
        int startPower = qFloor(logMin);
        int endPower = qCeil(logMax);

        for (int power = startPower; power <= endPower; ++power) {
            double value = qPow(10, power);
            if (value >= settings.xMin && value <= settings.xMax) {
                xTickValues.append(value);
            }
            // 添加中间刻度
            if (power < endPower && (endPower - startPower) <= 3) {
                for (int mult : {2, 3, 5, 7}) {
                    double midValue = mult * qPow(10, power);
                    if (midValue > settings.xMin && midValue < settings.xMax) {
                        xTickValues.append(midValue);
                    }
                }
            }
        }
    } else {
        // 线性X轴
        double range = settings.xMax - settings.xMin;
        double orderOfMagnitude = pow(10, floor(log10(range)));
        double normalizedRange = range / orderOfMagnitude;

        double tickInterval;
        if (normalizedRange <= 2) {
            tickInterval = orderOfMagnitude * 0.2;
        } else if (normalizedRange <= 5) {
            tickInterval = orderOfMagnitude * 0.5;
        } else {
            tickInterval = orderOfMagnitude;
        }

        double firstTick = ceil(settings.xMin / tickInterval) * tickInterval;

        for (double value = firstTick; value <= settings.xMax; value += tickInterval) {
            if (value >= settings.xMin) {
                xTickValues.append(value);
            }
        }
    }

    // Y轴网格线 - 参考Saphir软件的对数坐标处理
    if (settings.yAxisType == AxisType::Logarithmic && settings.yMin > 0) {
        // 对数Y轴
        double logMin = log10(settings.yMin);
        double logMax = log10(settings.yMax);
        int startPower = qFloor(logMin);
        int endPower = qCeil(logMax);

        for (int power = startPower; power <= endPower; ++power) {
            double value = qPow(10, power);
            if (value >= settings.yMin && value <= settings.yMax) {
                yTickValues.append(value);
            }
            // 添加中间刻度
            if (power < endPower && (endPower - startPower) <= 3) {
                for (int mult : {2, 3, 5, 7}) {
                    double midValue = mult * qPow(10, power);
                    if (midValue > settings.yMin && midValue < settings.yMax) {
                        yTickValues.append(midValue);
                    }
                }
            }
        }
    } else {
        // 线性Y轴
        double range = settings.yMax - settings.yMin;
        double orderOfMagnitude = pow(10, floor(log10(range)));
        double normalizedRange = range / orderOfMagnitude;

        double tickInterval;
        if (normalizedRange <= 2) {
            tickInterval = orderOfMagnitude * 0.2;
        } else if (normalizedRange <= 5) {
            tickInterval = orderOfMagnitude * 0.5;
        } else {
            tickInterval = orderOfMagnitude;
        }

        double firstTick = ceil(settings.yMin / tickInterval) * tickInterval;

        for (double value = firstTick; value <= settings.yMax; value += tickInterval) {
            if (value >= settings.yMin) {
                yTickValues.append(value);
            }
        }
    }

    // 绘制X轴网格线 - 使用正确的对数转换
    for (double value : xTickValues) {
        double x;
        if (settings.xAxisType == AxisType::Logarithmic && settings.xMin > 0 && settings.xMax > 0) {
            double normalizedX = (log10(value) - log10(settings.xMin)) / (log10(settings.xMax) - log10(settings.xMin));
            x = plotArea.left() + normalizedX * plotArea.width();
        } else {
            x = plotArea.left() + (value - settings.xMin) / (settings.xMax - settings.xMin) * plotArea.width();
        }

        if (x >= plotArea.left() && x <= plotArea.right()) {
            painter.drawLine(x, plotArea.top(), x, plotArea.bottom());
        }
    }

    // 绘制Y轴网格线 - 使用正确的对数转换
    for (double value : yTickValues) {
        double y;
        if (settings.yAxisType == AxisType::Logarithmic && settings.yMin > 0 && settings.yMax > 0) {
            double normalizedY = (log10(value) - log10(settings.yMin)) / (log10(settings.yMax) - log10(settings.yMin));
            y = plotArea.bottom() - normalizedY * plotArea.height();
        } else {
            y = plotArea.bottom() - (value - settings.yMin) / (settings.yMax - settings.yMin) * plotArea.height();
        }

        if (y >= plotArea.top() && y <= plotArea.bottom()) {
            painter.drawLine(plotArea.left(), y, plotArea.right(), y);
        }
    }
}

// 优化的坐标轴绘制函数，参考Saphir软件
void DualPlotWindow::drawAxesOnWidget(QPainter &painter, const QRect &plotArea, const PlotSettings &settings)
{
    painter.setPen(QPen(Qt::black, 1));
    painter.setFont(QFont("Arial", 8));

    // 生成优化的X轴标签
    QVector<double> xTickValues;
    if (settings.xAxisType == AxisType::Logarithmic && settings.xMin > 0) {
        // 对数坐标系 - 只显示主要的10的整数次幂
        double logMin = log10(settings.xMin);
        double logMax = log10(settings.xMax);
        int startPower = qFloor(logMin);
        int endPower = qCeil(logMax);

        for (int power = startPower; power <= endPower; ++power) {
            double value = qPow(10, power);
            if (value >= settings.xMin && value <= settings.xMax) {
                xTickValues.append(value);
            }
        }
    } else {
        // 线性坐标系
        double range = settings.xMax - settings.xMin;
        double orderOfMagnitude = pow(10, floor(log10(range)));
        double normalizedRange = range / orderOfMagnitude;

        double tickInterval;
        if (normalizedRange <= 2) {
            tickInterval = orderOfMagnitude * 0.2;
        } else if (normalizedRange <= 5) {
            tickInterval = orderOfMagnitude * 0.5;
        } else {
            tickInterval = orderOfMagnitude;
        }

        double firstTick = ceil(settings.xMin / tickInterval) * tickInterval;

        for (double value = firstTick; value <= settings.xMax; value += tickInterval) {
            if (value >= settings.xMin) {
                xTickValues.append(value);
            }
        }
    }

    // 绘制X轴标签
    for (double value : xTickValues) {
        double x;
        if (settings.xAxisType == AxisType::Logarithmic && settings.xMin > 0 && settings.xMax > 0) {
            double normalizedX = (log10(value) - log10(settings.xMin)) / (log10(settings.xMax) - log10(settings.xMin));
            x = plotArea.left() + normalizedX * plotArea.width();
        } else {
            x = plotArea.left() + (value - settings.xMin) / (settings.xMax - settings.xMin) * plotArea.width();
        }

        if (x >= plotArea.left() && x <= plotArea.right()) {
            painter.drawLine(x, plotArea.bottom(), x, plotArea.bottom() - 5);

            QString label;
            if (settings.xAxisType == AxisType::Logarithmic) {
                // 对数坐标：显示整数
                int power = qRound(log10(value));
                if (qAbs(log10(value) - power) < 0.01) {
                    if (power == 0) label = "1";
                    else if (power == 1) label = "10";
                    else if (power == 2) label = "100";
                    else if (power == 3) label = "1000";
                    else label = QString("10^%1").arg(power);
                } else {
                    label = QString::number(value, 'g', 2);
                }
            } else {
                // 线性坐标：显示10的倍数
                if (qAbs(value) >= 10000) {
                    label = QString::number(value, 'e', 1);
                } else if (qAbs(value - qRound(value)) < 0.001) {
                    label = QString::number(qRound(value));
                } else {
                    label = QString::number(value, 'f', 1);
                }
            }

            QRect textRect(x - 30, plotArea.bottom() + 5, 60, 15);
            painter.drawText(textRect, Qt::AlignCenter, label);
        }
    }

    // 生成优化的Y轴标签
    QVector<double> yTickValues;
    if (settings.yAxisType == AxisType::Logarithmic && settings.yMin > 0) {
        // 对数坐标系 - 只显示主要的10的整数次幂
        double logMin = log10(settings.yMin);
        double logMax = log10(settings.yMax);
        int startPower = qFloor(logMin);
        int endPower = qCeil(logMax);

        for (int power = startPower; power <= endPower; ++power) {
            double value = qPow(10, power);
            if (value >= settings.yMin && value <= settings.yMax) {
                yTickValues.append(value);
            }
        }
    } else {
        // 线性坐标系
        double range = settings.yMax - settings.yMin;
        double orderOfMagnitude = pow(10, floor(log10(range)));
        double normalizedRange = range / orderOfMagnitude;

        double tickInterval;
        if (normalizedRange <= 2) {
            tickInterval = orderOfMagnitude * 0.2;
        } else if (normalizedRange <= 5) {
            tickInterval = orderOfMagnitude * 0.5;
        } else {
            tickInterval = orderOfMagnitude;
        }

        double firstTick = ceil(settings.yMin / tickInterval) * tickInterval;

        for (double value = firstTick; value <= settings.yMax; value += tickInterval) {
            if (value >= settings.yMin) {
                yTickValues.append(value);
            }
        }
    }

    // 绘制Y轴标签
    for (double value : yTickValues) {
        double y;
        if (settings.yAxisType == AxisType::Logarithmic && settings.yMin > 0 && settings.yMax > 0) {
            double normalizedY = (log10(value) - log10(settings.yMin)) / (log10(settings.yMax) - log10(settings.yMin));
            y = plotArea.bottom() - normalizedY * plotArea.height();
        } else {
            y = plotArea.bottom() - (value - settings.yMin) / (settings.yMax - settings.yMin) * plotArea.height();
        }

        if (y >= plotArea.top() && y <= plotArea.bottom()) {
            painter.drawLine(plotArea.left(), y, plotArea.left() + 5, y);

            QString label;
            if (settings.yAxisType == AxisType::Logarithmic) {
                // 对数坐标：显示整数
                int power = qRound(log10(value));
                if (qAbs(log10(value) - power) < 0.01) {
                    if (power == 0) label = "1";
                    else if (power == 1) label = "10";
                    else if (power == 2) label = "100";
                    else if (power == 3) label = "1000";
                    else label = QString("10^%1").arg(power);
                } else {
                    label = QString::number(value, 'g', 2);
                }
            } else {
                // 线性坐标：显示10的倍数
                if (qAbs(value) >= 10000) {
                    label = QString::number(value, 'e', 1);
                } else if (qAbs(value - qRound(value)) < 0.001) {
                    label = QString::number(qRound(value));
                } else {
                    label = QString::number(value, 'f', 1);
                }
            }

            QRect textRect(plotArea.left() - 60, y - 8, 55, 16);
            painter.drawText(textRect, Qt::AlignRight | Qt::AlignVCenter, label);
        }
    }

    // 绘制轴标题
    painter.setFont(QFont("Arial", 10, QFont::Bold));

    // X轴标题
    QRect xTitleRect(plotArea.left(), plotArea.bottom() + 25, plotArea.width(), 15);
    painter.drawText(xTitleRect, Qt::AlignCenter, settings.xAxisTitle);

    // Y轴标题
    painter.save();
    painter.translate(15, plotArea.center().y());
    painter.rotate(-90);
    painter.drawText(-60, -3, 120, 15, Qt::AlignCenter, settings.yAxisTitle);
    painter.restore();

    // 图标题
    painter.setFont(QFont("Arial", 11, QFont::Bold));
    QRect titleRect(plotArea.left(), 5, plotArea.width(), 30);
    painter.drawText(titleRect, Qt::AlignCenter, settings.plotTitle);
}

void DualPlotWindow::drawLegendOnWidget(QPainter &painter, const QVector<CurveData> &curves,
                                        const QRect &plotArea, const QPoint &offset)
{
    if (curves.isEmpty()) return;

    int visibleCurveCount = 0;
    for (const CurveData &curve : curves) {
        if (curve.visible) {
            visibleCurveCount++;
        }
    }

    if (visibleCurveCount == 0) return;

    painter.setFont(QFont("Arial", 9));
    QFontMetrics fm(painter.font());

    int lineHeight = fm.height() + 4;
    int legendWidth = 0;
    int legendHeight = visibleCurveCount * lineHeight + 16;

    for (const CurveData &curve : curves) {
        if (curve.visible) {
            int textWidth = fm.horizontalAdvance(curve.name) + 40;
            legendWidth = qMax(legendWidth, textWidth);
        }
    }
    legendWidth += 20;

    int legendX = plotArea.right() - legendWidth - 10 + offset.x();
    int legendY = plotArea.top() + 10 + offset.y();

    legendX = qMax(plotArea.left() + 10, qMin(legendX, plotArea.right() - legendWidth - 10));
    legendY = qMax(plotArea.top() + 10, qMin(legendY, plotArea.bottom() - legendHeight - 10));

    QRect legendArea(legendX, legendY, legendWidth, legendHeight);

    painter.fillRect(legendArea, QColor(255, 255, 255, 200));
    painter.setPen(QPen(QColor(150, 150, 150, 180), 1));
    painter.drawRect(legendArea);

    painter.setPen(QPen(Qt::black, 1));
    painter.setFont(QFont("Arial", 9, QFont::Bold));
    painter.drawText(legendX + 8, legendY + 15, "图例");

    painter.setFont(QFont("Arial", 8));
    int currentY = legendY + 25;

    for (const CurveData &curve : curves) {
        if (!curve.visible) continue;

        Qt::PenStyle penStyle = lineStyleToQt(curve.lineStyle);
        painter.setPen(QPen(curve.color, qMax(1, curve.lineWidth - 1), penStyle));
        painter.drawLine(legendX + 8, currentY + lineHeight/2 - 2,
                         legendX + 28, currentY + lineHeight/2 - 2);

        painter.setBrush(curve.color);
        painter.setPen(QPen(curve.color, 1));
        painter.drawEllipse(legendX + 18 - 2, currentY + lineHeight/2 - 4, 4, 4);

        painter.setPen(QPen(Qt::black, 1));
        painter.drawText(legendX + 35, currentY + lineHeight/2 + 3, curve.name);

        currentY += lineHeight;
    }
}

// 鼠标事件处理函数
void DualPlotWindow::onPressureMousePress(QMouseEvent *event)
{
    QPoint pos = event->pos();
    m_lastPressureMousePos = pos;

    if (!m_pressurePlotArea.contains(pos)) {
        return;
    }

    if (event->button() == Qt::LeftButton) {
        if (event->modifiers() & Qt::ControlModifier) {
            m_pressureSelecting = true;
            m_pressureSelectionStart = pos;
            m_pressureSelectionRect = QRect(pos, pos);
        } else {
            m_pressureDragging = true;
            m_pressurePanning = true;
        }
    }
}

void DualPlotWindow::onPressureMouseMove(QMouseEvent *event)
{
    QPoint pos = event->pos();

    if (m_pressureSelecting) {
        m_pressureSelectionRect = QRect(m_pressureSelectionStart, pos).normalized();
        updatePlots();
    } else if (m_pressureDragging && m_pressurePanning) {
        QPointF delta = pos - m_lastPressureMousePos;
        panPressureView(delta);

        if (m_syncPan) {
            panProductionView(delta);
        }

        m_lastPressureMousePos = pos;
    }
}

void DualPlotWindow::onPressureMouseRelease(QMouseEvent *event)
{
    Q_UNUSED(event)

    if (m_pressureSelecting && !m_pressureSelectionRect.isEmpty()) {
        QPointF topLeft = pixelToPressureData(m_pressureSelectionRect.topLeft());
        QPointF bottomRight = pixelToPressureData(m_pressureSelectionRect.bottomRight());

        m_pressureSettings.xMin = qMin(topLeft.x(), bottomRight.x());
        m_pressureSettings.xMax = qMax(topLeft.x(), bottomRight.x());
        m_pressureSettings.yMin = qMin(topLeft.y(), bottomRight.y());
        m_pressureSettings.yMax = qMax(topLeft.y(), bottomRight.y());

        if (m_syncZoom) {
            m_productionSettings.xMin = m_pressureSettings.xMin;
            m_productionSettings.xMax = m_pressureSettings.xMax;
        }
    }

    m_pressureDragging = false;
    m_pressureSelecting = false;
    m_pressurePanning = false;
    updatePlots();
}

void DualPlotWindow::onPressureWheel(QWheelEvent *event)
{
    QPoint pos = event->position().toPoint();

    if (m_pressurePlotArea.contains(pos)) {
        double factor = 1.0 + event->angleDelta().y() / 1200.0;
        zoomPressureAtPoint(pos, factor);

        if (m_syncZoom) {
            // 同步X轴缩放到产量图
            m_productionSettings.xMin = m_pressureSettings.xMin;
            m_productionSettings.xMax = m_pressureSettings.xMax;
            updatePlots();
        }
    }
}

void DualPlotWindow::onProductionMousePress(QMouseEvent *event)
{
    QPoint pos = event->pos();
    m_lastProductionMousePos = pos;

    if (!m_productionPlotArea.contains(pos)) {
        return;
    }

    if (event->button() == Qt::LeftButton) {
        if (event->modifiers() & Qt::ControlModifier) {
            m_productionSelecting = true;
            m_productionSelectionStart = pos;
            m_productionSelectionRect = QRect(pos, pos);
        } else {
            m_productionDragging = true;
            m_productionPanning = true;
        }
    }
}

void DualPlotWindow::onProductionMouseMove(QMouseEvent *event)
{
    QPoint pos = event->pos();

    if (m_productionSelecting) {
        m_productionSelectionRect = QRect(m_productionSelectionStart, pos).normalized();
        updatePlots();
    } else if (m_productionDragging && m_productionPanning) {
        QPointF delta = pos - m_lastProductionMousePos;
        panProductionView(delta);

        if (m_syncPan) {
            panPressureView(delta);
        }

        m_lastProductionMousePos = pos;
    }
}

void DualPlotWindow::onProductionMouseRelease(QMouseEvent *event)
{
    Q_UNUSED(event)

    if (m_productionSelecting && !m_productionSelectionRect.isEmpty()) {
        QPointF topLeft = pixelToProductionData(m_productionSelectionRect.topLeft());
        QPointF bottomRight = pixelToProductionData(m_productionSelectionRect.bottomRight());

        m_productionSettings.xMin = qMin(topLeft.x(), bottomRight.x());
        m_productionSettings.xMax = qMax(topLeft.x(), bottomRight.x());
        m_productionSettings.yMin = qMin(topLeft.y(), bottomRight.y());
        m_productionSettings.yMax = qMax(topLeft.y(), bottomRight.y());

        if (m_syncZoom) {
            m_pressureSettings.xMin = m_productionSettings.xMin;
            m_pressureSettings.xMax = m_productionSettings.xMax;
        }
    }

    m_productionDragging = false;
    m_productionSelecting = false;
    m_productionPanning = false;
    updatePlots();
}

void DualPlotWindow::onProductionWheel(QWheelEvent *event)
{
    QPoint pos = event->position().toPoint();

    if (m_productionPlotArea.contains(pos)) {
        double factor = 1.0 + event->angleDelta().y() / 1200.0;
        zoomProductionAtPoint(pos, factor);

        if (m_syncZoom) {
            // 同步X轴缩放到压力图
            m_pressureSettings.xMin = m_productionSettings.xMin;
            m_pressureSettings.xMax = m_productionSettings.xMax;
            updatePlots();
        }
    }
}

// 坐标转换函数 - 修正对数坐标
QPointF DualPlotWindow::pressureDataToPixel(const QPointF &dataPoint)
{
    double x, y;

    if (m_pressureSettings.xMax > m_pressureSettings.xMin) {
        if (m_pressureSettings.xAxisType == AxisType::Logarithmic && dataPoint.x() > 0 && m_pressureSettings.xMin > 0) {
            double normalizedX = (log10(dataPoint.x()) - log10(m_pressureSettings.xMin)) /
                                 (log10(m_pressureSettings.xMax) - log10(m_pressureSettings.xMin));
            x = m_pressurePlotArea.left() + normalizedX * m_pressurePlotArea.width();
        } else {
            x = m_pressurePlotArea.left() + (dataPoint.x() - m_pressureSettings.xMin) /
                                                (m_pressureSettings.xMax - m_pressureSettings.xMin) * m_pressurePlotArea.width();
        }
    } else {
        x = m_pressurePlotArea.left();
    }

    if (m_pressureSettings.yMax > m_pressureSettings.yMin) {
        if (m_pressureSettings.yAxisType == AxisType::Logarithmic && dataPoint.y() > 0 && m_pressureSettings.yMin > 0) {
            double normalizedY = (log10(dataPoint.y()) - log10(m_pressureSettings.yMin)) /
                                 (log10(m_pressureSettings.yMax) - log10(m_pressureSettings.yMin));
            y = m_pressurePlotArea.bottom() - normalizedY * m_pressurePlotArea.height();
        } else {
            y = m_pressurePlotArea.bottom() - (dataPoint.y() - m_pressureSettings.yMin) /
                                                  (m_pressureSettings.yMax - m_pressureSettings.yMin) * m_pressurePlotArea.height();
        }
    } else {
        y = m_pressurePlotArea.bottom();
    }

    return QPointF(x, y);
}

QPointF DualPlotWindow::pixelToPressureData(const QPointF &pixelPoint)
{
    double x, y;

    if (m_pressureSettings.xAxisType == AxisType::Logarithmic && m_pressureSettings.xMin > 0 && m_pressureSettings.xMax > 0) {
        double normalizedX = (pixelPoint.x() - m_pressurePlotArea.left()) / m_pressurePlotArea.width();
        x = pow(10, log10(m_pressureSettings.xMin) + normalizedX * (log10(m_pressureSettings.xMax) - log10(m_pressureSettings.xMin)));
    } else {
        x = m_pressureSettings.xMin + (pixelPoint.x() - m_pressurePlotArea.left()) / m_pressurePlotArea.width() *
                                          (m_pressureSettings.xMax - m_pressureSettings.xMin);
    }

    if (m_pressureSettings.yAxisType == AxisType::Logarithmic && m_pressureSettings.yMin > 0 && m_pressureSettings.yMax > 0) {
        double normalizedY = (m_pressurePlotArea.bottom() - pixelPoint.y()) / m_pressurePlotArea.height();
        y = pow(10, log10(m_pressureSettings.yMin) + normalizedY * (log10(m_pressureSettings.yMax) - log10(m_pressureSettings.yMin)));
    } else {
        y = m_pressureSettings.yMin + (m_pressurePlotArea.bottom() - pixelPoint.y()) / m_pressurePlotArea.height() *
                                          (m_pressureSettings.yMax - m_pressureSettings.yMin);
    }

    return QPointF(x, y);
}

QPointF DualPlotWindow::productionDataToPixel(const QPointF &dataPoint)
{
    double x, y;

    if (m_productionSettings.xMax > m_productionSettings.xMin) {
        if (m_productionSettings.xAxisType == AxisType::Logarithmic && dataPoint.x() > 0 && m_productionSettings.xMin > 0) {
            double normalizedX = (log10(dataPoint.x()) - log10(m_productionSettings.xMin)) /
                                 (log10(m_productionSettings.xMax) - log10(m_productionSettings.xMin));
            x = m_productionPlotArea.left() + normalizedX * m_productionPlotArea.width();
        } else {
            x = m_productionPlotArea.left() + (dataPoint.x() - m_productionSettings.xMin) /
                                                  (m_productionSettings.xMax - m_productionSettings.xMin) * m_productionPlotArea.width();
        }
    } else {
        x = m_productionPlotArea.left();
    }

    if (m_productionSettings.yMax > m_productionSettings.yMin) {
        if (m_productionSettings.yAxisType == AxisType::Logarithmic && dataPoint.y() > 0 && m_productionSettings.yMin > 0) {
            double normalizedY = (log10(dataPoint.y()) - log10(m_productionSettings.yMin)) /
                                 (log10(m_productionSettings.yMax) - log10(m_productionSettings.yMin));
            y = m_productionPlotArea.bottom() - normalizedY * m_productionPlotArea.height();
        } else {
            y = m_productionPlotArea.bottom() - (dataPoint.y() - m_productionSettings.yMin) /
                                                    (m_productionSettings.yMax - m_productionSettings.yMin) * m_productionPlotArea.height();
        }
    } else {
        y = m_productionPlotArea.bottom();
    }

    return QPointF(x, y);
}

QPointF DualPlotWindow::pixelToProductionData(const QPointF &pixelPoint)
{
    double x, y;

    if (m_productionSettings.xAxisType == AxisType::Logarithmic && m_productionSettings.xMin > 0 && m_productionSettings.xMax > 0) {
        double normalizedX = (pixelPoint.x() - m_productionPlotArea.left()) / m_productionPlotArea.width();
        x = pow(10, log10(m_productionSettings.xMin) + normalizedX * (log10(m_productionSettings.xMax) - log10(m_productionSettings.xMin)));
    } else {
        x = m_productionSettings.xMin + (pixelPoint.x() - m_productionPlotArea.left()) / m_productionPlotArea.width() *
                                            (m_productionSettings.xMax - m_productionSettings.xMin);
    }

    if (m_productionSettings.yAxisType == AxisType::Logarithmic && m_productionSettings.yMin > 0 && m_productionSettings.yMax > 0) {
        double normalizedY = (m_productionPlotArea.bottom() - pixelPoint.y()) / m_productionPlotArea.height();
        y = pow(10, log10(m_productionSettings.yMin) + normalizedY * (log10(m_productionSettings.yMax) - log10(m_productionSettings.yMin)));
    } else {
        y = m_productionSettings.yMin + (m_productionPlotArea.bottom() - pixelPoint.y()) / m_productionPlotArea.height() *
                                            (m_productionSettings.yMax - m_productionSettings.yMin);
    }

    return QPointF(x, y);
}

// 缩放和平移函数
void DualPlotWindow::zoomPressureAtPoint(const QPointF &point, double factor)
{
    QPointF dataPoint = pixelToPressureData(point);

    double xRange = m_pressureSettings.xMax - m_pressureSettings.xMin;
    double yRange = m_pressureSettings.yMax - m_pressureSettings.yMin;

    double newXRange = xRange / factor;
    double newYRange = yRange / factor;

    m_pressureSettings.xMin = dataPoint.x() - newXRange * (dataPoint.x() - m_pressureSettings.xMin) / xRange;
    m_pressureSettings.xMax = dataPoint.x() + newXRange * (m_pressureSettings.xMax - dataPoint.x()) / xRange;
    m_pressureSettings.yMin = dataPoint.y() - newYRange * (dataPoint.y() - m_pressureSettings.yMin) / yRange;
    m_pressureSettings.yMax = dataPoint.y() + newYRange * (m_pressureSettings.yMax - dataPoint.y()) / yRange;

    updatePlots();
}

void DualPlotWindow::zoomProductionAtPoint(const QPointF &point, double factor)
{
    QPointF dataPoint = pixelToProductionData(point);

    double xRange = m_productionSettings.xMax - m_productionSettings.xMin;
    double yRange = m_productionSettings.yMax - m_productionSettings.yMin;

    double newXRange = xRange / factor;
    double newYRange = yRange / factor;

    m_productionSettings.xMin = dataPoint.x() - newXRange * (dataPoint.x() - m_productionSettings.xMin) / xRange;
    m_productionSettings.xMax = dataPoint.x() + newXRange * (m_productionSettings.xMax - dataPoint.x()) / xRange;
    m_productionSettings.yMin = dataPoint.y() - newYRange * (dataPoint.y() - m_productionSettings.yMin) / yRange;
    m_productionSettings.yMax = dataPoint.y() + newYRange * (m_productionSettings.yMax - dataPoint.y()) / yRange;

    updatePlots();
}

void DualPlotWindow::panPressureView(const QPointF &delta)
{
    QPointF dataDelta = pixelToPressureData(delta) - pixelToPressureData(QPointF(0, 0));

    m_pressureSettings.xMin -= dataDelta.x();
    m_pressureSettings.xMax -= dataDelta.x();
    m_pressureSettings.yMin -= dataDelta.y();
    m_pressureSettings.yMax -= dataDelta.y();

    updatePlots();
}

void DualPlotWindow::panProductionView(const QPointF &delta)
{
    QPointF dataDelta = pixelToProductionData(delta) - pixelToProductionData(QPointF(0, 0));

    m_productionSettings.xMin -= dataDelta.x();
    m_productionSettings.xMax -= dataDelta.x();
    m_productionSettings.yMin -= dataDelta.y();
    m_productionSettings.yMax -= dataDelta.y();

    updatePlots();
}

// 增强的数据自适应功能 - 特别针对产量数据
void DualPlotWindow::synchronizeXAxis()
{
    // 同步两个图的X轴范围
    double minX = 1e10, maxX = -1e10;
    bool hasValidData = false;

    // 检查压力数据
    for (const CurveData &curve : m_pressureCurves) {
        for (int i = 0; i < curve.xData.size(); ++i) {
            double x = curve.xData[i];
            if (qIsFinite(x) && x > 0) {  // 对数坐标需要正值
                minX = qMin(minX, x);
                maxX = qMax(maxX, x);
                hasValidData = true;
            }
        }
    }

    // 检查产量数据
    for (const CurveData &curve : m_productionCurves) {
        for (int i = 0; i < curve.xData.size(); ++i) {
            double x = curve.xData[i];
            if (qIsFinite(x) && x > 0) {  // 对数坐标需要正值
                minX = qMin(minX, x);
                maxX = qMax(maxX, x);
                hasValidData = true;
            }
        }
    }

    if (hasValidData && minX < maxX) {
        // 根据坐标类型计算合适的范围
        if (m_pressureSettings.xAxisType == AxisType::Logarithmic) {
            // 对数坐标系
            double logMin = log10(minX);
            double logMax = log10(maxX);
            m_pressureSettings.xMin = pow(10, floor(logMin));
            m_pressureSettings.xMax = pow(10, ceil(logMax));

            // 如果范围太小，扩展一个数量级
            if (m_pressureSettings.xMax / m_pressureSettings.xMin < 10) {
                m_pressureSettings.xMin /= 10;
                m_pressureSettings.xMax *= 10;
            }
        } else {
            // 线性坐标系
            double range = maxX - minX;
            double orderOfMagnitude = pow(10, floor(log10(range)));

            double minTick = floor(minX / orderOfMagnitude) * orderOfMagnitude;
            double maxTick = ceil(maxX / orderOfMagnitude) * orderOfMagnitude;

            // 添加5%的边距
            double margin = (maxTick - minTick) * 0.05;
            m_pressureSettings.xMin = minTick - margin;
            m_pressureSettings.xMax = maxTick + margin;
        }

        m_productionSettings.xMin = m_pressureSettings.xMin;
        m_productionSettings.xMax = m_pressureSettings.xMax;
    }
}

void DualPlotWindow::calculatePressureBounds()
{
    if (m_pressureCurves.isEmpty()) return;

    double minY = 1e10, maxY = -1e10;
    bool hasValidData = false;

    for (const CurveData &curve : m_pressureCurves) {
        for (int i = 0; i < curve.yData.size(); ++i) {
            double y = curve.yData[i];
            if (qIsFinite(y) && y > 0) {  // 对数坐标需要正值
                minY = qMin(minY, y);
                maxY = qMax(maxY, y);
                hasValidData = true;
            }
        }
    }

    if (hasValidData && minY < maxY) {
        if (m_pressureSettings.yAxisType == AxisType::Logarithmic) {
            // 对数坐标系
            double logMin = log10(minY);
            double logMax = log10(maxY);
            m_pressureSettings.yMin = pow(10, floor(logMin));
            m_pressureSettings.yMax = pow(10, ceil(logMax));

            // 如果范围太小，扩展一个数量级
            if (m_pressureSettings.yMax / m_pressureSettings.yMin < 10) {
                m_pressureSettings.yMin /= 10;
                m_pressureSettings.yMax *= 10;
            }
        } else {
            // 线性坐标系
            double range = maxY - minY;
            double orderOfMagnitude = pow(10, floor(log10(range)));

            double minTick = floor(minY / orderOfMagnitude) * orderOfMagnitude;
            double maxTick = ceil(maxY / orderOfMagnitude) * orderOfMagnitude;

            // 添加5%的边距
            double margin = (maxTick - minTick) * 0.05;
            m_pressureSettings.yMin = minTick - margin;
            m_pressureSettings.yMax = maxTick + margin;

            // 如果包含0，确保0对齐
            if (m_pressureSettings.yMin < 0 && m_pressureSettings.yMax > 0) {
                m_pressureSettings.yMin = floor(m_pressureSettings.yMin / orderOfMagnitude) * orderOfMagnitude;
                m_pressureSettings.yMax = ceil(m_pressureSettings.yMax / orderOfMagnitude) * orderOfMagnitude;
            }
        }
    }
}

void DualPlotWindow::calculateProductionBounds()
{
    if (m_productionCurves.isEmpty()) return;

    double minY = 1e10, maxY = -1e10;
    bool hasValidData = false;

    for (const CurveData &curve : m_productionCurves) {
        for (int i = 0; i < curve.yData.size(); ++i) {
            double y = curve.yData[i];
            if (qIsFinite(y) && y > 0) {  // 对数坐标需要正值
                minY = qMin(minY, y);
                maxY = qMax(maxY, y);
                hasValidData = true;
            }
        }
    }

    if (hasValidData && minY < maxY) {
        if (m_productionSettings.yAxisType == AxisType::Logarithmic) {
            // 对数坐标系
            double logMin = log10(minY);
            double logMax = log10(maxY);
            m_productionSettings.yMin = pow(10, floor(logMin));
            m_productionSettings.yMax = pow(10, ceil(logMax));

            // 如果范围太小，扩展一个数量级
            if (m_productionSettings.yMax / m_productionSettings.yMin < 10) {
                m_productionSettings.yMin /= 10;
                m_productionSettings.yMax *= 10;
            }
        } else {
            // 线性坐标系 - 产量数据最小值固定为0
            double range = maxY - 0;  // 产量最小值为0
            double orderOfMagnitude = pow(10, floor(log10(range)));

            double maxTick = ceil(maxY / orderOfMagnitude) * orderOfMagnitude;

            // 添加5%的边距
            double margin = maxTick * 0.05;
            m_productionSettings.yMin = 0;  // 产量最小值固定为0
            m_productionSettings.yMax = maxTick + margin;
        }
    }
}

// 优化的曲线绘制函数 - 正确处理对数坐标
void DualPlotWindow::drawCurveOnWidget(QPainter &painter, const CurveData &curve,
                                       const QRect &plotArea, const PlotSettings &settings)
{
    if (curve.drawType == CurveType::Step) {
        // 阶梯状曲线绘制
        Qt::PenStyle penStyle = lineStyleToQt(curve.lineStyle);
        painter.setPen(QPen(curve.color, curve.lineWidth, penStyle));

        painter.setClipRect(plotArea.adjusted(-5, -5, 5, 5));

        int dataSize = qMin(curve.xData.size(), curve.yData.size());

        for (int i = 0; i < dataSize - 1; i += 2) {
            if (i + 1 >= dataSize) break;

            // 数据有效性检查
            if (!qIsFinite(curve.xData[i]) || !qIsFinite(curve.yData[i]) ||
                !qIsFinite(curve.xData[i+1]) || !qIsFinite(curve.yData[i+1])) {
                continue;
            }

            // 对数坐标需要正值
            if (settings.xAxisType == AxisType::Logarithmic &&
                (curve.xData[i] <= 0 || curve.xData[i+1] <= 0)) continue;
            if (settings.yAxisType == AxisType::Logarithmic &&
                (curve.yData[i] <= 0 || curve.yData[i+1] <= 0)) continue;

            // 坐标转换
            double x1, y1, x2, y2;

            if (settings.xAxisType == AxisType::Logarithmic) {
                x1 = (log10(curve.xData[i]) - log10(settings.xMin)) /
                     (log10(settings.xMax) - log10(settings.xMin));
                x2 = (log10(curve.xData[i+1]) - log10(settings.xMin)) /
                     (log10(settings.xMax) - log10(settings.xMin));
            } else {
                x1 = (curve.xData[i] - settings.xMin) / (settings.xMax - settings.xMin);
                x2 = (curve.xData[i+1] - settings.xMin) / (settings.xMax - settings.xMin);
            }

            if (settings.yAxisType == AxisType::Logarithmic) {
                y1 = (log10(curve.yData[i]) - log10(settings.yMin)) /
                     (log10(settings.yMax) - log10(settings.yMin));
                y2 = (log10(curve.yData[i+1]) - log10(settings.yMin)) /
                     (log10(settings.yMax) - log10(settings.yMin));
            } else {
                y1 = (curve.yData[i] - settings.yMin) / (settings.yMax - settings.yMin);
                y2 = (curve.yData[i+1] - settings.yMin) / (settings.yMax - settings.yMin);
            }

            double pixelX1 = plotArea.left() + x1 * plotArea.width();
            double pixelY1 = plotArea.bottom() - y1 * plotArea.height();
            double pixelX2 = plotArea.left() + x2 * plotArea.width();
            double pixelY2 = plotArea.bottom() - y2 * plotArea.height();

            QPointF segmentStart(pixelX1, pixelY1);
            QPointF segmentEnd(pixelX2, pixelY2);

            painter.drawLine(segmentStart, segmentEnd);

            // 绘制垂直连接线
            if (i + 2 < dataSize && i + 3 < dataSize) {
                if (qIsFinite(curve.xData[i+2]) && qIsFinite(curve.yData[i+2])) {
                    // 对数坐标检查
                    if (settings.xAxisType == AxisType::Logarithmic && curve.xData[i+2] <= 0) continue;
                    if (settings.yAxisType == AxisType::Logarithmic && curve.yData[i+2] <= 0) continue;

                    double x3, y3;
                    if (settings.xAxisType == AxisType::Logarithmic) {
                        x3 = (log10(curve.xData[i+2]) - log10(settings.xMin)) /
                             (log10(settings.xMax) - log10(settings.xMin));
                    } else {
                        x3 = (curve.xData[i+2] - settings.xMin) / (settings.xMax - settings.xMin);
                    }

                    if (settings.yAxisType == AxisType::Logarithmic) {
                        y3 = (log10(curve.yData[i+2]) - log10(settings.yMin)) /
                             (log10(settings.yMax) - log10(settings.yMin));
                    } else {
                        y3 = (curve.yData[i+2] - settings.yMin) / (settings.yMax - settings.yMin);
                    }

                    double pixelX3 = plotArea.left() + x3 * plotArea.width();
                    double pixelY3 = plotArea.bottom() - y3 * plotArea.height();
                    QPointF nextSegmentStart(pixelX3, pixelY3);
                    QPointF verticalStart = segmentEnd;
                    QPointF verticalEnd = QPointF(segmentEnd.x(), nextSegmentStart.y());
                    painter.drawLine(verticalStart, verticalEnd);
                }
            }
        }

        painter.setClipping(false);

        // 绘制数据点
        painter.setBrush(curve.color);
        for (int i = 0; i < dataSize; i += 2) {
            if (i < dataSize && qIsFinite(curve.xData[i]) && qIsFinite(curve.yData[i])) {
                // 对数坐标检查
                if (settings.xAxisType == AxisType::Logarithmic && curve.xData[i] <= 0) continue;
                if (settings.yAxisType == AxisType::Logarithmic && curve.yData[i] <= 0) continue;

                double x, y;
                if (settings.xAxisType == AxisType::Logarithmic) {
                    x = (log10(curve.xData[i]) - log10(settings.xMin)) /
                        (log10(settings.xMax) - log10(settings.xMin));
                } else {
                    x = (curve.xData[i] - settings.xMin) / (settings.xMax - settings.xMin);
                }

                if (settings.yAxisType == AxisType::Logarithmic) {
                    y = (log10(curve.yData[i]) - log10(settings.yMin)) /
                        (log10(settings.yMax) - log10(settings.yMin));
                } else {
                    y = (curve.yData[i] - settings.yMin) / (settings.yMax - settings.yMin);
                }

                double pixelX = plotArea.left() + x * plotArea.width();
                double pixelY = plotArea.bottom() - y * plotArea.height();
                QPointF pixelPoint(pixelX, pixelY);
                if (plotArea.contains(pixelPoint.toPoint())) {
                    painter.drawEllipse(pixelPoint, curve.pointSize/2, curve.pointSize/2);
                }
            }
        }
        return;
    }

    // 普通曲线绘制
    Qt::PenStyle penStyle = lineStyleToQt(curve.lineStyle);
    painter.setPen(QPen(curve.color, curve.lineWidth, penStyle));

    QVector<QPointF> points;
    int maxPoints = qMin(curve.xData.size(), curve.yData.size());

    for (int i = 0; i < maxPoints; ++i) {
        if (qIsFinite(curve.xData[i]) && qIsFinite(curve.yData[i])) {
            // 对数坐标需要正值
            if (settings.xAxisType == AxisType::Logarithmic && curve.xData[i] <= 0) continue;
            if (settings.yAxisType == AxisType::Logarithmic && curve.yData[i] <= 0) continue;

            double x, y;

            // X轴坐标转换
            if (settings.xAxisType == AxisType::Logarithmic) {
                x = (log10(curve.xData[i]) - log10(settings.xMin)) /
                    (log10(settings.xMax) - log10(settings.xMin));
            } else {
                x = (curve.xData[i] - settings.xMin) / (settings.xMax - settings.xMin);
            }

            // Y轴坐标转换
            if (settings.yAxisType == AxisType::Logarithmic) {
                y = (log10(curve.yData[i]) - log10(settings.yMin)) /
                    (log10(settings.yMax) - log10(settings.yMin));
            } else {
                y = (curve.yData[i] - settings.yMin) / (settings.yMax - settings.yMin);
            }

            double pixelX = plotArea.left() + x * plotArea.width();
            double pixelY = plotArea.bottom() - y * plotArea.height();

            points.append(QPointF(pixelX, pixelY));
        }
    }

    // 绘制曲线
    if (points.size() > 1) {
        painter.setClipRect(plotArea.adjusted(-5, -5, 5, 5));
        for (int i = 0; i < points.size() - 1; ++i) {
            painter.drawLine(points[i], points[i + 1]);
        }
        painter.setClipping(false);
    }

    // 绘制数据点
    painter.setBrush(curve.color);
    for (const QPointF &point : points) {
        if (plotArea.contains(point.toPoint())) {
            painter.drawEllipse(point, curve.pointSize/2, curve.pointSize/2);
        }
    }
}
