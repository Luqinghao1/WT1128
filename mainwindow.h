#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QMap>
#include <QTimer>
#include <QStandardItemModel>
#include "modelmanager.h" // [修复] 必须包含，因为用到了 ModelType

class NavBtn;
class MonitorWidget;
class DataEditorWidget;
class PlottingWidget;
class FittingWidget;
class SettingsWidget;

struct WellTestData;

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    void init();

    void initMonitorForm();
    void initDataEditorForm();
    void initModelForm();
    void initPlottingForm();
    void initFittingForm();

private slots:
    void onProjectCreated();
    void onFileLoaded(const QString& filePath, const QString& fileType);
    void onPlotAnalysisCompleted(const QString &analysisType, const QMap<QString, double> &results);
    void onDataReadyForPlotting();
    void onTransferDataToPlotting();
    void onDataEditorDataChanged();
    void onSystemSettingsChanged();
    void onAutoSaveIntervalChanged(int interval);
    void onBackupSettingsChanged(bool enabled);
    void onPerformanceSettingsChanged();
    void onModelCalculationCompleted(const QString &analysisType, const QMap<QString, double> &results);

    // [修复] 参数类型 ModelManager::ModelType
    void onFittingCompleted(ModelManager::ModelType modelType, const QMap<QString, double> &parameters);

    void onFittingProgressChanged(int progress);

private:
    Ui::MainWindow *ui;
    MonitorWidget* m_MonitorWidget;
    DataEditorWidget* m_DataEditorWidget;
    ModelManager* m_ModelManager;
    PlottingWidget* m_PlottingWidget;
    FittingWidget* m_FittingWidget;
    SettingsWidget* m_SettingsWidget;
    QMap<QString,NavBtn*> m_NavBtnMap;
    QTimer m_timer;
    bool m_hasValidData = false;

    void transferDataFromEditorToPlotting();
    void updateNavigationState();
    void transferDataToFitting();

    QStandardItemModel* getDataEditorModel() const;
    QString getCurrentFileName() const;
    bool hasDataLoaded();
    WellTestData createDemoWellTestData();
};

#endif // MAINWINDOW_H
