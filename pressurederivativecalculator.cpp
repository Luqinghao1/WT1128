#include "PressureDerivativeCalculator.h"
#include <QStandardItem>
#include <QRegularExpression>
#include <QDebug>
#include <cmath>

PressureDerivativeCalculator::PressureDerivativeCalculator(QObject *parent)
    : QObject(parent)
{
}

PressureDerivativeCalculator::~PressureDerivativeCalculator()
{
}

PressureDerivativeResult PressureDerivativeCalculator::calculatePressureDerivative(
    QStandardItemModel* model, const PressureDerivativeConfig& config)
{
    PressureDerivativeResult result;
    result.success = false;
    result.addedColumnIndex = -1;
    result.processedRows = 0;

    // 检查数据模型
    if (!model) {
        result.errorMessage = "数据模型不存在";
        return result;
    }

    int rowCount = model->rowCount();
    if (rowCount < 3) {
        result.errorMessage = "数据行数不足（至少需要3行）";
        return result;
    }

    // 检查列索引
    if (config.pressureColumnIndex < 0 || config.pressureColumnIndex >= model->columnCount()) {
        result.errorMessage = "压力列索引无效";
        return result;
    }

    if (config.timeColumnIndex < 0 || config.timeColumnIndex >= model->columnCount()) {
        result.errorMessage = "时间列索引无效";
        return result;
    }

    // 检查L-Spacing参数
    if (config.lSpacing <= 0) {
        result.errorMessage = "L-Spacing参数必须大于0";
        return result;
    }

    emit progressUpdated(10, "正在读取数据...");

    // 读取时间和压力数据 (使用 QVector 提高性能)
    QVector<double> timeData;
    QVector<double> pressureData;
    timeData.reserve(rowCount);
    pressureData.reserve(rowCount);

    for (int row = 0; row < rowCount; ++row) {
        QStandardItem* timeItem = model->item(row, config.timeColumnIndex);
        QStandardItem* pressureItem = model->item(row, config.pressureColumnIndex);

        double timeValue = 0.0;
        double pressureValue = 0.0;

        if (timeItem) {
            timeValue = parseNumericValue(timeItem->text());
        }
        if (pressureItem) {
            pressureValue = parseNumericValue(pressureItem->text());
        }

        // 检查时间值有效性（允许从0开始）
        if (timeValue < 0) {
            result.errorMessage = QString("检测到无效时间值（行 %1），时间不能为负数").arg(row + 1);
            return result;
        }

        timeData.append(timeValue);
        pressureData.append(pressureValue);
    }

    // 检查是否需要添加时间偏移（处理t=0的情况）
    double actualTimeOffset = 0.0;
    if (config.autoTimeOffset) {
        // 查找最小非零时间值
        double minPositiveTime = -1;
        bool hasZeroTime = false;

        for (double t : timeData) {
            if (t <= 0) {
                hasZeroTime = true;
            } else {
                if (minPositiveTime < 0 || t < minPositiveTime) {
                    minPositiveTime = t;
                }
            }
        }

        // 如果有零或负时间值，添加偏移
        if (hasZeroTime) {
            if (minPositiveTime > 0) {
                // 使用最小正时间值的1/10作为偏移
                actualTimeOffset = minPositiveTime * 0.1;
            } else {
                // 如果所有时间都<=0，使用配置的偏移量
                actualTimeOffset = config.timeOffset;
            }
            qDebug() << "检测到时间从0开始，自动添加时间偏移：" << actualTimeOffset;
        }
    } else {
        actualTimeOffset = config.timeOffset;
    }

    // 应用时间偏移
    QVector<double> adjustedTimeData;
    adjustedTimeData.reserve(rowCount);
    for (double t : timeData) {
        adjustedTimeData.append(t + actualTimeOffset);
    }

    emit progressUpdated(30, "正在计算压降...");

    // 计算压降 (初始压力 - 当前压力，假定是压降测试)
    // 注意：如果是恢复测试(Buildup)，通常需要输入Delta t和Delta P。
    // 这里保持原逻辑：Delta P = P_initial - P_current
    QVector<double> pressureDropData;
    pressureDropData.reserve(rowCount);
    double initialPressure = pressureData.isEmpty() ? 0.0 : pressureData[0];

    for (int i = 0; i < pressureData.size(); ++i) {
        double pressureDrop = initialPressure - pressureData[i];
        pressureDropData.append(pressureDrop);
    }

    emit progressUpdated(50, "正在计算Bourdet导数（L-Spacing平滑）...");

    // 调用静态统一算法
    QVector<double> derivativeData = calculateBourdetDerivative(adjustedTimeData, pressureDropData, config.lSpacing);

    if (derivativeData.size() != rowCount) {
        result.errorMessage = "导数计算结果数量不匹配";
        return result;
    }

    emit progressUpdated(80, "正在写入结果...");

    // 在压力列后面插入新列
    int newColumnIndex = config.pressureColumnIndex + 1;
    model->insertColumn(newColumnIndex);

    // 设置列标题
    QString columnName = QString("压力导数\\%1").arg(config.pressureUnit);
    QStandardItem* headerItem = new QStandardItem(columnName);
    model->setHorizontalHeaderItem(newColumnIndex, headerItem);

    // 写入导数数据
    for (int row = 0; row < rowCount; ++row) {
        QString valueText = formatValue(derivativeData[row], 6);
        QStandardItem* item = new QStandardItem(valueText);
        item->setForeground(QBrush(QColor("#1565C0"))); // 蓝色文字
        model->setItem(row, newColumnIndex, item);
        result.processedRows++;
    }

    emit progressUpdated(100, "计算完成");

    // 设置返回结果
    result.success = true;
    result.addedColumnIndex = newColumnIndex;
    result.columnName = columnName;

    emit calculationCompleted(result);

    return result;
}

// 静态方法实现：Bourdet 导数核心算法 (Saphir 方法)
QVector<double> PressureDerivativeCalculator::calculateBourdetDerivative(
    const QVector<double>& timeData,
    const QVector<double>& pressureDropData,
    double lSpacing)
{
    QVector<double> derivativeData;
    int n = timeData.size();
    derivativeData.reserve(n);

    if (n == 0) return derivativeData;

    for (int i = 0; i < n; ++i) {
        double derivative = 0.0;
        double ti = timeData[i];
        double pi = pressureDropData[i];

        // 寻找左侧点j：ln(ti) - ln(tj) ≥ L
        int leftIndex = findLeftPoint(timeData, i, lSpacing);

        // 寻找右侧点k：ln(tk) - ln(ti) ≥ L
        int rightIndex = findRightPoint(timeData, i, lSpacing);

        // 1. 如果找到左右两个点，使用加权平均法 (Bourdet Standard)
        if (leftIndex >= 0 && rightIndex >= 0) {
            double tj = timeData[leftIndex];
            double pj = pressureDropData[leftIndex];
            double tk = timeData[rightIndex];
            double pk = pressureDropData[rightIndex];

            // 计算对数差值
            double deltaXL = std::log(ti) - std::log(tj);  // ΔXL = ln(ti) - ln(tj)
            double deltaXR = std::log(tk) - std::log(ti);  // ΔXR = ln(tk) - ln(ti)

            // 计算左导数和右导数
            double mL = calculateDerivativeValue(ti, tj, pi, pj);  // 左导数 slope
            double mR = calculateDerivativeValue(tk, ti, pk, pi);  // 右导数 slope

            // 加权平均公式：P' = (mL * ΔXR + mR * ΔXL) / (ΔXL + ΔXR)
            if (deltaXL + deltaXR > 1e-12) {
                derivative = (mL * deltaXR + mR * deltaXL) / (deltaXL + deltaXR);
            } else {
                derivative = 0.0;
            }
        }
        // 2. 边界情况：只找到左侧点 (曲线末端)
        else if (leftIndex >= 0 && rightIndex < 0) {
            double tj = timeData[leftIndex];
            double pj = pressureDropData[leftIndex];
            derivative = calculateDerivativeValue(ti, tj, pi, pj);
        }
        // 3. 边界情况：只找到右侧点 (曲线开端)
        else if (leftIndex < 0 && rightIndex >= 0) {
            double tk = timeData[rightIndex];
            double pk = pressureDropData[rightIndex];
            derivative = calculateDerivativeValue(tk, ti, pk, pi);
        }
        // 4. L-Spacing 范围内点不足 (通常是数据极少或 L 设置过大)
        else {
            // 使用简单的相邻点差分作为保底
            if (i > 0) {
                double t_prev = timeData[i-1];
                double p_prev = pressureDropData[i-1];
                derivative = calculateDerivativeValue(ti, t_prev, pi, p_prev);
            } else if (i < n - 1) {
                double t_next = timeData[i+1];
                double p_next = pressureDropData[i+1];
                derivative = calculateDerivativeValue(t_next, ti, p_next, pi);
            } else {
                derivative = 0.0;
            }
        }

        derivativeData.append(derivative);
    }

    return derivativeData;
}

int PressureDerivativeCalculator::findLeftPoint(const QVector<double>& timeData, int currentIndex, double lSpacing)
{
    if (currentIndex <= 0 || timeData.isEmpty()) return -1;

    double ti = timeData[currentIndex];
    if (ti <= 0) return -1;
    double lnTi = std::log(ti);

    // 从当前点向左搜索，找到第一个满足距离 >= L 的点
    for (int j = currentIndex - 1; j >= 0; --j) {
        double tj = timeData[j];
        if (tj <= 0) continue;

        double lnTj = std::log(tj);
        if ((lnTi - lnTj) >= lSpacing) {
            return j;
        }
    }
    return -1;
}

int PressureDerivativeCalculator::findRightPoint(const QVector<double>& timeData, int currentIndex, double lSpacing)
{
    int n = timeData.size();
    if (currentIndex >= n - 1 || timeData.isEmpty()) return -1;

    double ti = timeData[currentIndex];
    if (ti <= 0) return -1;
    double lnTi = std::log(ti);

    // 从当前点向右搜索，找到第一个满足距离 >= L 的点
    for (int k = currentIndex + 1; k < n; ++k) {
        double tk = timeData[k];
        if (tk <= 0) continue;

        double lnTk = std::log(tk);
        if ((lnTk - lnTi) >= lSpacing) {
            return k;
        }
    }
    return -1;
}

double PressureDerivativeCalculator::calculateDerivativeValue(double t1, double t2, double p1, double p2)
{
    // 计算单边导数：dP/d(ln t) = (p1 - p2) / (ln(t1) - ln(t2))
    if (t1 <= 0 || t2 <= 0) return 0.0;

    double lnT1 = std::log(t1);
    double lnT2 = std::log(t2);
    double deltaLnT = lnT1 - lnT2;

    if (std::abs(deltaLnT) < 1e-10) {
        return 0.0;
    }

    return (p1 - p2) / deltaLnT;
}

PressureDerivativeConfig PressureDerivativeCalculator::autoDetectColumns(QStandardItemModel* model)
{
    PressureDerivativeConfig config;
    if (!model) return config;

    config.pressureColumnIndex = findPressureColumn(model);
    config.timeColumnIndex = findTimeColumn(model);
    return config;
}

int PressureDerivativeCalculator::findPressureColumn(QStandardItemModel* model)
{
    if (!model) return -1;
    QStringList pressureKeywords = {"压力", "pressure", "pres", "P\\", "压力\\"};

    for (int col = 0; col < model->columnCount(); ++col) {
        QStandardItem* headerItem = model->horizontalHeaderItem(col);
        if (headerItem) {
            QString headerText = headerItem->text();
            for (const QString& keyword : pressureKeywords) {
                if (headerText.contains(keyword, Qt::CaseInsensitive)) {
                    if (!headerText.contains("压降") && !headerText.contains("导数")) {
                        return col;
                    }
                }
            }
        }
    }
    return -1;
}

int PressureDerivativeCalculator::findTimeColumn(QStandardItemModel* model)
{
    if (!model) return -1;
    QStringList timeKeywords = {"时间", "time", "t\\", "小时", "hour", "min", "sec"};

    for (int col = 0; col < model->columnCount(); ++col) {
        QStandardItem* headerItem = model->horizontalHeaderItem(col);
        if (headerItem) {
            QString headerText = headerItem->text();
            for (const QString& keyword : timeKeywords) {
                if (headerText.contains(keyword, Qt::CaseInsensitive)) {
                    return col;
                }
            }
        }
    }
    return -1;
}

double PressureDerivativeCalculator::parseNumericValue(const QString& str)
{
    if (str.isEmpty()) return 0.0;
    QString cleanStr = str.trimmed();
    bool ok;
    double value = cleanStr.toDouble(&ok);
    if (ok) return value;

    cleanStr.remove(QRegularExpression("[a-zA-Z%\\s]+$"));
    value = cleanStr.toDouble(&ok);
    return ok ? value : 0.0;
}

QString PressureDerivativeCalculator::formatValue(double value, int precision)
{
    if (std::isnan(value) || std::isinf(value)) return "0";
    return QString::number(value, 'g', precision);
}
