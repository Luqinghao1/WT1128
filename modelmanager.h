#ifndef MODELMANAGER_H
#define MODELMANAGER_H

#include <QObject>
#include <QWidget>
#include <QStackedWidget>
#include <QComboBox>
#include <QMap>
#include <QVector>
#include <tuple>

// 定义模型曲线数据类型: <时间, 压力, 导数>
typedef std::tuple<QVector<double>, QVector<double>, QVector<double>> ModelCurveData;

class ModelWidget1;
class ModelWidget2;
class ModelWidget3;

class ModelManager : public QObject
{
    Q_OBJECT

public:
    enum ModelType {
        InfiniteConductive = 0,    // 无限导流
        FiniteConductive = 1,      // 有限导流
        SegmentedMultiCluster = 2  // 分段多簇
    };
    Q_ENUM(ModelType)

    explicit ModelManager(QWidget* parent = nullptr);
    ~ModelManager();

    // UI初始化
    void initializeModels(QWidget* parentWidget);
    QWidget* getMainWidget() const { return m_mainWidget; }

    // 模型切换
    void switchToModel(ModelType modelType);
    ModelType getCurrentModelType() const { return m_currentModelType; }

    // 静态辅助信息
    static QString getModelTypeName(ModelType type);
    static QStringList getAvailableModelTypes();

    // =========================================================================
    // 核心计算接口
    // =========================================================================

    // 获取指定模型的默认参数
    QMap<QString, double> getDefaultParameters(ModelType type);

    // 计算理论曲线 (返回: 时间, 压力, 导数)
    ModelCurveData calculateTheoreticalCurve(ModelType type, const QMap<QString, double>& params);

signals:
    void modelSwitched(ModelType newType, ModelType oldType);
    void calculationCompleted(const QString& title, const QMap<QString, double>& results);

private slots:
    void onModelTypeSelectionChanged(int index);
    void onModel1CalculationCompleted(const QString& t, const QMap<QString, double>& r);
    void onModel2CalculationCompleted(const QString& t, const QMap<QString, double>& r);
    void onModel3CalculationCompleted(const QString& t, const QMap<QString, double>& r);

private:
    void createMainWidget();
    void setupModelSelection();
    void connectModelSignals();

    // 具体模型计算实现
    ModelCurveData calculateModel1(const QMap<QString, double>& params);
    ModelCurveData calculateModel2(const QMap<QString, double>& params);
    ModelCurveData calculateModel3(const QMap<QString, double>& params);

    // 使用 PressureDerivativeCalculator 替代原来的手动计算
    void computePressureDerivative(const QVector<double>& tD, const QVector<double>& pd,
                                   double cD, QVector<double>& dpd, QVector<double>& td_dpd);

    // 数学辅助函数
    double flaplace1(double z, const QMap<QString, double>& p);
    double flaplace2(double z, const QMap<QString, double>& p);
    double e_function(double z, int i, int j, int k, int v, int mf, int nf, double omega, double lambda, double Xf, double yy, double y);
    double f_function(int j, int nf, double Xf, double y);
    double integralBesselK0(double XDkv, double YDkv, double yDij, double fz, double a, double b);
    double gaussQuadrature(double XDkv, double YDkv, double yDij, double fz, double a, double b);
    double besselK0(double x);
    QVector<double> solveLinearSystem(const QVector<QVector<double>>& A, const QVector<double>& b);
    double stefestCoefficient(int i, int N);
    double factorial(int n);

    // UI 组件
    QWidget* m_mainWidget;
    QComboBox* m_modelTypeCombo;
    QStackedWidget* m_modelStack;

    ModelWidget1* m_modelWidget1;
    ModelWidget2* m_modelWidget2;
    ModelWidget3* m_modelWidget3;

    ModelType m_currentModelType;
};

#endif // MODELMANAGER_H
