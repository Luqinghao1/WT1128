#ifndef FITTINGWIDGET_H
#define FITTINGWIDGET_H

#include <QWidget>
#include <QMap>
#include <QVector>
#include <QFuture>
#include <QFutureWatcher>
#include <QDialog>
#include <QComboBox>
#include <QTableWidget>
#include <QLabel>
#include <QPushButton>
#include "modelmanager.h"
#include "qcustomplot.h"

namespace Ui { class FittingWidget; }

// ===========================================================================
// 数据加载对话框类
// ===========================================================================
class FittingDataLoadDialog : public QDialog
{
    Q_OBJECT
public:
    explicit FittingDataLoadDialog(const QList<QStringList>& previewData, QWidget *parent = nullptr);
    int getTimeColumnIndex() const;
    int getPressureColumnIndex() const;
    int getDerivativeColumnIndex() const;
    int getSkipRows() const;
private slots:
    void validateSelection();
private:
    QTableWidget* m_previewTable;
    QComboBox* m_comboTime, *m_comboPressure, *m_comboDeriv, *m_comboSkipRows;
};

// ===========================================================================
// 拟合参数结构体
// ===========================================================================
struct FitParameter {
    QString name;       // 内部参数名 (如 "omega")
    QString displayName;// 显示名称 (如 "储容比 (ω)")
    double value;       // 当前值
    double min;         // 最小值下限
    double max;         // 最大值上限
    bool isFit;         // 是否参与拟合
};

// ===========================================================================
// 拟合主窗口类
// ===========================================================================
class FittingWidget : public QWidget
{
    Q_OBJECT
public:
    explicit FittingWidget(QWidget *parent = nullptr);
    ~FittingWidget();

    void setModelManager(ModelManager* manager);
    void setObservedData(const QVector<double>& time, const QVector<double>& pressure, const QVector<double>& derivative);

signals:
    void fittingCompleted(ModelManager::ModelType modelType, const QMap<QString, double> &parameters);
    void sigIterationUpdated(double currentError, const QMap<QString, double>& currentParams);
    void sigProgress(int value);

private slots:
    void on_comboModelSelect_currentIndexChanged(int index);
    void on_btnResetParams_clicked();
    void on_btnRunFit_clicked();
    void on_btnStop_clicked();
    void on_btnUpdateModel_clicked();
    void on_btnLoadData_clicked();

    void onIterationUpdate(double currentError, const QMap<QString, double>& currentParams);
    void onFitFinished();

private:
    Ui::FittingWidget *ui;
    QCustomPlot *m_plot;
    ModelManager* m_modelManager;

    QVector<double> m_obsTime, m_obsPressure, m_obsDerivative;
    QList<FitParameter> m_parameters;

    bool m_isFitting;
    bool m_stopRequested;
    QFutureWatcher<void> m_watcher;

    void setupPlot();
    void initModelCombo();
    void loadParamsToTable();
    void updateParamsFromTable();

    void plotCurves(const QVector<double>& t, const QVector<double>& p, const QVector<double>& dp, bool isModel);
    QString getParamDisplayName(const QString& key);

    // --- 核心优化逻辑 (修复：传递参数副本以实现线程安全) ---
    void runOptimizationTask(int algIndex, ModelManager::ModelType modelType, QList<FitParameter> fitParams);

    // 具体算法实现 (修复：接收参数列表)
    void runLevenbergMarquardtOptimization(ModelManager::ModelType modelType, QList<FitParameter> params);
    void runNelderMeadOptimization(ModelManager::ModelType modelType, QList<FitParameter> params);

    QVector<double> calculateResiduals(const QMap<QString, double>& params, ModelManager::ModelType modelType);
    double calculateSumSquaredError(const QVector<double>& residuals);

    QVector<QVector<double>> computeJacobian(const QMap<QString, double>& params, const QVector<double>& currentResiduals, const QVector<int>& fitIndices, ModelManager::ModelType modelType, const QList<FitParameter>& currentFitParams);
    QVector<double> solveLinearSystem(const QVector<QVector<double>>& A, const QVector<double>& b);

    QStringList parseLine(const QString& line);
    double calculateError(const QMap<QString, double>& trialParams, ModelManager::ModelType modelType);
};
#endif // FITTINGWIDGET_H
