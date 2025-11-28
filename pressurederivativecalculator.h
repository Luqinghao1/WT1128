#ifndef PRESSUREDERIVATIVECALCULATOR_H
#define PRESSUREDERIVATIVECALCULATOR_H

#include <QObject>
#include <QString>
#include <QVector>
#include <QStandardItemModel>

// 压力导数计算结果结构
struct PressureDerivativeResult {
    bool success;
    QString errorMessage;
    int addedColumnIndex;
    QString columnName;
    int processedRows;

    PressureDerivativeResult() :
        success(false),
        addedColumnIndex(-1),
        processedRows(0) {}
};

// 压力导数计算配置
struct PressureDerivativeConfig {
    int timeColumnIndex;      // 时间列索引
    int pressureColumnIndex;  // 压力列索引
    QString timeUnit;         // 时间单位 ("s", "min", "h")
    QString pressureUnit;     // 压力单位
    double lSpacing;          // L-Spacing平滑参数（对数周期，通常0.1-0.5）
    double timeOffset;        // 时间偏移量（用于处理t=0的情况）
    bool autoTimeOffset;      // 是否自动添加时间偏移

    PressureDerivativeConfig() :
        timeColumnIndex(-1),
        pressureColumnIndex(-1),
        timeUnit("h"),
        pressureUnit("MPa"),
        lSpacing(0.15),        // 默认值0.15个对数周期 (Saphir标准默认值附近)
        timeOffset(0.0001),    // 默认偏移量0.0001
        autoTimeOffset(true) {}// 默认自动添加偏移
};

/**
 * @brief 压力导数计算器类
 *
 * 使用Bourdet导数算法（L-Spacing平滑算法）计算压力导数：
 * P' = dP/d(ln t) = t * dP/dt
 *
 * 核心算法已提取为静态方法，供 FittingWidget 等模块复用。
 */
class PressureDerivativeCalculator : public QObject
{
    Q_OBJECT

public:
    explicit PressureDerivativeCalculator(QObject *parent = nullptr);
    ~PressureDerivativeCalculator();

    /**
     * @brief 计算压力导数（针对表格模型的封装）
     * @param model 数据模型
     * @param config 计算配置
     * @return 计算结果
     */
    PressureDerivativeResult calculatePressureDerivative(QStandardItemModel* model,
                                                         const PressureDerivativeConfig& config);

    /**
     * @brief 自动检测压力列和时间列
     * @param model 数据模型
     * @return 配置对象，包含检测到的列索引
     */
    PressureDerivativeConfig autoDetectColumns(QStandardItemModel* model);

    // =========================================================================
    // 静态核心算法接口 (Saphir 风格 Bourdet 导数)
    // =========================================================================

    /**
     * @brief 使用Bourdet导数算法计算导数
     * @param timeData 时间数据
     * @param pressureDropData 压降数据 (Delta P)
     * @param lSpacing L-Spacing参数 (通常0.1-0.5)
     * @return 导数数据向量
     */
    static QVector<double> calculateBourdetDerivative(const QVector<double>& timeData,
                                                      const QVector<double>& pressureDropData,
                                                      double lSpacing);

signals:
    void progressUpdated(int progress, const QString& message);
    void calculationCompleted(const PressureDerivativeResult& result);

private:
    // 内部静态辅助函数
    static int findLeftPoint(const QVector<double>& timeData, int currentIndex, double lSpacing);
    static int findRightPoint(const QVector<double>& timeData, int currentIndex, double lSpacing);
    static double calculateDerivativeValue(double t1, double t2, double p1, double p2);

    int findPressureColumn(QStandardItemModel* model);
    int findTimeColumn(QStandardItemModel* model);
    double parseNumericValue(const QString& str);
    QString formatValue(double value, int precision = 6);
};

#endif // PRESSUREDERIVATIVECALCULATOR_H
