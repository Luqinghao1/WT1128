#ifndef MODELMANAGER_H
#define MODELMANAGER_H

#include <QObject>
#include <QWidget>
#include <QComboBox>
#include <QStackedWidget>
#include <QMap>
#include <tuple>
#include <QVector>

// Forward declarations
class ModelWidget1;
class ModelWidget2;
class ModelWidget3;

// Define the return type for curve data: (Time, Pressure, Derivative)
typedef std::tuple<QVector<double>, QVector<double>, QVector<double>> ModelCurveData;

class ModelManager : public QObject
{
    Q_OBJECT
public:
    enum ModelType {
        InfiniteConductive = 0,
        FiniteConductive,
        SegmentedMultiCluster
    };
    Q_ENUM(ModelType)

    explicit ModelManager(QWidget* parent = nullptr);
    ~ModelManager();

    void initializeModels(QWidget* parentWidget);

    // Get current model type
    ModelType getCurrentModelType() const { return m_currentModelType; }

    // Get available model names
    static QStringList getAvailableModelTypes();
    static QString getModelTypeName(ModelType type);

    // Get default parameters for fitting
    QMap<QString, double> getDefaultParameters(ModelType type);

    // Calculate theoretical curve based on parameters
    ModelCurveData calculateTheoreticalCurve(ModelType type, const QMap<QString, double>& params);

signals:
    void modelSwitched(ModelType newType, ModelType oldType);
    void calculationCompleted(const QString& analysisType, const QMap<QString, double>& results);

private slots:
    void onModelTypeSelectionChanged(int index);

    // Slots to receive results from individual model widgets
    void onModel1CalculationCompleted(const QString& type, const QMap<QString, double>& results);
    void onModel2CalculationCompleted(const QString& type, const QMap<QString, double>& results);
    void onModel3CalculationCompleted(const QString& type, const QMap<QString, double>& results);

private:
    QWidget* m_mainWidget;
    QComboBox* m_modelTypeCombo;
    QStackedWidget* m_modelStack;

    ModelWidget1* m_modelWidget1;
    ModelWidget2* m_modelWidget2;
    ModelWidget3* m_modelWidget3;

    ModelType m_currentModelType;

    void createMainWidget();
    void setupModelSelection();
    void connectModelSignals();
    void switchToModel(ModelType modelType);

    // --- Mathematical Calculation Core ---

    // Calculate curves for specific models
    ModelCurveData calculateModel1(const QMap<QString, double>& params);
    ModelCurveData calculateModel2(const QMap<QString, double>& params);
    ModelCurveData calculateModel3(const QMap<QString, double>& params);

    // Generic Derivative Calculation
    void computePressureDerivative(const QVector<double>& tD, const QVector<double>& pd,
                                   double cD, QVector<double>& dpd, QVector<double>& td_dpd);

    // Numerical Inversion (Stehfest) helper
    double stefestCoefficient(int i, int N);
    double factorial(int n);

    // Laplace space functions (Model Specific)
    // [FIX] Added declarations for specific flaplace functions
    double flaplace(double z, const QMap<QString, double>& params); // Keep generic if needed, or remove
    double flaplace1(double z, const QMap<QString, double>& p); // Infinite Conductivity
    double flaplace2(double z, const QMap<QString, double>& p); // Finite Conductivity
    double flaplace3(double z, const QMap<QString, double>& p); // Segmented Multi-Cluster

    // Bessel functions and Integration
    double e_function(double z, int i, int j, int k, int v, int mf, int nf, double omega, double lambda, double Xf, double yy, double y);
    double f_function(int j, int nf, double Xf, double y);

    double integralBesselK0(double XDkv, double YDkv, double yDij, double fz, double a, double b);
    double gaussQuadrature(double XDkv, double YDkv, double yDij, double fz, double a, double b);
    double besselK0(double x);

    // Linear Algebra
    QVector<double> solveLinearSystem(const QVector<QVector<double>>& A, const QVector<double>& b);
};

#endif // MODELMANAGER_H
