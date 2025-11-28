#ifndef FITTINGWIDGET_H
#define FITTINGWIDGET_H

#include <QWidget>
#include <QFutureWatcher>
#include "modelmanager.h"
#include "qcustomplot.h"

namespace Ui {
class FittingWidget;
}

// 拟合参数结构
struct FitParameter {
    QString name;
    QString displayName;
    double value;
    double min;
    double max;
    bool isFit;
};

class FittingWidget : public QWidget
{
    Q_OBJECT

public:
    explicit FittingWidget(QWidget *parent = nullptr);
    ~FittingWidget();

    void setModelManager(ModelManager* m);
    void setObservedData(const QVector<double>& t, const QVector<double>& p, const QVector<double>& d);

signals:
    void sigIterationUpdated(double error, const QMap<QString,double>& currentParams);
    void sigProgress(int percent);
    void fittingCompleted(ModelManager::ModelType type, const QMap<QString,double>& finalParams);

private slots:
    void on_btnLoadData_clicked();
    void on_btnRunFit_clicked();
    void on_btnStop_clicked();
    void on_btnResetParams_clicked();
    void on_btnUpdateModel_clicked();
    void on_comboModelSelect_currentIndexChanged(int index);

    void onIterationUpdate(double err, const QMap<QString,double>& p);
    void onFitFinished();

private:
    Ui::FittingWidget *ui;
    ModelManager* m_modelManager;
    QCustomPlot* m_plot;
    bool m_isFitting;
    bool m_stopRequested;

    // 观测数据
    QVector<double> m_obsTime;
    QVector<double> m_obsPressure;
    QVector<double> m_obsDerivative;

    // 拟合参数
    QList<FitParameter> m_parameters;

    // 异步任务
    QFutureWatcher<void> m_watcher;

    void setupPlot();
    void initModelCombo();
    void loadParamsToTable();
    void updateParamsFromTable();
    QString getParamDisplayName(const QString& key);
    QStringList parseLine(const QString& line);
    void plotCurves(const QVector<double>& t, const QVector<double>& p, const QVector<double>& d, bool isModel);

    // 拟合核心逻辑
    void runOptimizationTask(ModelManager::ModelType modelType, QList<FitParameter> fitParams); // 移除算法索引参数

    // 仅保留 Levenberg-Marquardt (仿 lsqcurvefit)
    void runLevenbergMarquardtOptimization(ModelManager::ModelType modelType, QList<FitParameter> params);

    // 数学辅助
    QVector<double> calculateResiduals(const QMap<QString,double>& params, ModelManager::ModelType modelType);
    double calculateSumSquaredError(const QVector<double>& residuals);
    QVector<QVector<double>> computeJacobian(const QMap<QString,double>& params,
                                             const QVector<double>& baseResiduals,
                                             const QVector<int>& fitIndices,
                                             ModelManager::ModelType modelType,
                                             const QList<FitParameter>& currentFitParams);

    QVector<double> solveLinearSystem(const QVector<QVector<double>>& A, const QVector<double>& b);
    double calculateError(const QMap<QString,double>& trialParams, ModelManager::ModelType modelType);
};

// 数据加载对话框 (保持不变)
class QTableWidget;
class QComboBox;

class FittingDataLoadDialog : public QDialog {
    Q_OBJECT
public:
    FittingDataLoadDialog(const QList<QStringList>& previewData, QWidget* parent=nullptr);
    int getTimeColumnIndex() const;
    int getPressureColumnIndex() const;
    int getDerivativeColumnIndex() const;
    int getSkipRows() const;
private slots:
    void validateSelection();
private:
    QTableWidget* m_previewTable;
    QComboBox* m_comboTime;
    QComboBox* m_comboPressure;
    QComboBox* m_comboDeriv;
    QComboBox* m_comboSkipRows;
};

#endif // FITTINGWIDGET_H
