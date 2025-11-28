#include "fittingwidget.h"
#include "ui_fittingwidget.h"
#include "PressureDerivativeCalculator.h"
#include <QtConcurrent>
#include <QMessageBox>
#include <QDebug>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <QFileDialog>
#include <QFile>
#include <QTextStream>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGroupBox>
#include <QHeaderView>


// ===========================================================================
// FittingDataLoadDialog 实现
// ===========================================================================

FittingDataLoadDialog::FittingDataLoadDialog(const QList<QStringList>& previewData, QWidget *parent)
    : QDialog(parent)
{
    setWindowTitle("数据列映射配置");
    resize(800, 500);
    setStyleSheet("QDialog { background-color: #f0f0f0; } QLabel, QComboBox, QPushButton, QTableWidget, QGroupBox { color: black; }");

    QVBoxLayout* layout = new QVBoxLayout(this);
    layout->addWidget(new QLabel("请指定数据列含义 (时间必选):", this));

    m_previewTable = new QTableWidget(this);
    if(!previewData.isEmpty()) {
        int rows = qMin(previewData.size(), 50);
        int cols = previewData[0].size();
        m_previewTable->setRowCount(rows); m_previewTable->setColumnCount(cols);
        QStringList headers; for(int i=0;i<cols;++i) headers<<QString("Col %1").arg(i+1);
        m_previewTable->setHorizontalHeaderLabels(headers);
        for(int i=0;i<rows;++i) for(int j=0;j<cols && j<previewData[i].size();++j)
                m_previewTable->setItem(i,j,new QTableWidgetItem(previewData[i][j]));
    }
    m_previewTable->setAlternatingRowColors(true);
    layout->addWidget(m_previewTable);

    QGroupBox* grp = new QGroupBox("列映射", this);
    QGridLayout* grid = new QGridLayout(grp);
    QStringList opts; for(int i=0;i<m_previewTable->columnCount();++i) opts<<QString("Col %1").arg(i+1);

    grid->addWidget(new QLabel("时间 *:",this),0,0);
    m_comboTime = new QComboBox(this); m_comboTime->addItems(opts);
    grid->addWidget(m_comboTime,0,1);

    grid->addWidget(new QLabel("压力:",this),0,2);
    m_comboPressure = new QComboBox(this); m_comboPressure->addItem("不导入",-1); m_comboPressure->addItems(opts);
    if(opts.size()>1) m_comboPressure->setCurrentIndex(2);
    grid->addWidget(m_comboPressure,0,3);

    grid->addWidget(new QLabel("导数:",this),1,0);
    m_comboDeriv = new QComboBox(this); m_comboDeriv->addItem("自动计算 (Bourdet)",-1); m_comboDeriv->addItems(opts);
    grid->addWidget(m_comboDeriv,1,1);

    grid->addWidget(new QLabel("跳过行:",this),1,2);
    m_comboSkipRows = new QComboBox(this);
    for(int i=0;i<=20;++i) m_comboSkipRows->addItem(QString::number(i),i);
    m_comboSkipRows->setCurrentIndex(1);
    grid->addWidget(m_comboSkipRows,1,3);

    layout->addWidget(grp);

    QHBoxLayout* btns = new QHBoxLayout;
    QPushButton* ok = new QPushButton("确定",this);
    QPushButton* cancel = new QPushButton("取消",this);
    connect(ok, &QPushButton::clicked, this, &FittingDataLoadDialog::validateSelection);
    connect(cancel, &QPushButton::clicked, this, &QDialog::reject);
    btns->addStretch(); btns->addWidget(ok); btns->addWidget(cancel);
    layout->addLayout(btns);
}

void FittingDataLoadDialog::validateSelection() { if(m_comboTime->currentIndex()<0) return; accept(); }
int FittingDataLoadDialog::getTimeColumnIndex() const { return m_comboTime->currentIndex(); }
int FittingDataLoadDialog::getPressureColumnIndex() const { return m_comboPressure->currentIndex()-1; }
int FittingDataLoadDialog::getDerivativeColumnIndex() const { return m_comboDeriv->currentIndex()-1; }
int FittingDataLoadDialog::getSkipRows() const { return m_comboSkipRows->currentData().toInt(); }


// ===========================================================================
// FittingWidget 实现
// ===========================================================================

FittingWidget::FittingWidget(QWidget *parent) : QWidget(parent), ui(new Ui::FittingWidget), m_modelManager(nullptr), m_isFitting(false)
{
    ui->setupUi(this);
    // ... UI setup (same as before) ...
    this->setStyleSheet("QWidget { color: black; font-family: Arial; } "
                        "QGroupBox { font-weight: bold; border: 1px solid gray; margin-top: 10px; } "
                        "QGroupBox::title { subcontrol-origin: margin; left: 10px; padding: 0 3px; }");

    ui->splitter->setSizes(QList<int>{300, 800});
    ui->splitter->setCollapsible(0, false);

    ui->tableParams->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);

    m_plot = new QCustomPlot(this);
    ui->plotContainer->layout()->addWidget(m_plot);
    setupPlot();

    qRegisterMetaType<QMap<QString,double>>("QMap<QString,double>");
    qRegisterMetaType<ModelManager::ModelType>("ModelManager::ModelType");

    connect(this, &FittingWidget::sigIterationUpdated, this, &FittingWidget::onIterationUpdate, Qt::QueuedConnection);
    connect(this, &FittingWidget::sigProgress, ui->progressBar, &QProgressBar::setValue);
    connect(&m_watcher, &QFutureWatcher<void>::finished, this, &FittingWidget::onFitFinished);
}

FittingWidget::~FittingWidget() { delete ui; }

void FittingWidget::setModelManager(ModelManager *m) { m_modelManager = m; initModelCombo(); }

void FittingWidget::initModelCombo() {
    if(!m_modelManager) return;
    ui->comboModelSelect->clear();
    ui->comboModelSelect->addItems(ModelManager::getAvailableModelTypes());
    ui->comboModelSelect->setCurrentIndex((int)m_modelManager->getCurrentModelType());
}

void FittingWidget::setupPlot() {
    // ... (Same as before) ...
    m_plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom);
    m_plot->setBackground(Qt::white);
    m_plot->xAxis->setScaleType(QCPAxis::stLogarithmic);
    m_plot->yAxis->setScaleType(QCPAxis::stLogarithmic);
    QSharedPointer<QCPAxisTickerLog> logTicker(new QCPAxisTickerLog);
    m_plot->xAxis->setTicker(logTicker);
    m_plot->yAxis->setTicker(logTicker);

    m_plot->xAxis->setNumberFormat("eb"); m_plot->xAxis->setNumberPrecision(0);
    m_plot->xAxis->setLabel("tD / CD");
    m_plot->yAxis->setLabel("PD & dPD");

    m_plot->xAxis->setBasePen(QPen(Qt::black)); m_plot->xAxis->setTickPen(QPen(Qt::black)); m_plot->xAxis->setSubTickPen(QPen(Qt::black));
    m_plot->yAxis->setBasePen(QPen(Qt::black)); m_plot->yAxis->setTickPen(QPen(Qt::black)); m_plot->yAxis->setSubTickPen(QPen(Qt::black));
    m_plot->xAxis->setTickLabelColor(Qt::black); m_plot->xAxis->setLabelColor(Qt::black);
    m_plot->yAxis->setTickLabelColor(Qt::black); m_plot->yAxis->setLabelColor(Qt::black);

    m_plot->addGraph(); m_plot->graph(0)->setPen(Qt::NoPen); m_plot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, Qt::blue, 5)); m_plot->graph(0)->setName("实测压力");
    m_plot->addGraph(); m_plot->graph(1)->setPen(Qt::NoPen); m_plot->graph(1)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssTriangle, Qt::red, 5)); m_plot->graph(1)->setName("实测导数");
    m_plot->addGraph(); m_plot->graph(2)->setPen(QPen(Qt::black, 2)); m_plot->graph(2)->setName("理论压力");
    m_plot->addGraph(); m_plot->graph(3)->setPen(QPen(Qt::green, 2)); m_plot->graph(3)->setName("理论导数");

    m_plot->legend->setVisible(true);
    m_plot->legend->setBrush(QColor(255,255,255,200));
    m_plot->legend->setTextColor(Qt::black);
}

void FittingWidget::setObservedData(const QVector<double>& t, const QVector<double>& p, const QVector<double>& d) {
    // ... (Same as before) ...
    m_obsTime = t; m_obsPressure = p; m_obsDerivative = d;
    QVector<double> vt, vp, vd;
    for(int i=0; i<t.size(); ++i) {
        if(t[i]>1e-6 && p[i]>1e-6) {
            vt<<t[i]; vp<<p[i];
            if(i<d.size() && d[i]>1e-6) vd<<d[i]; else vd<<1e-10;
        }
    }
    m_plot->graph(0)->setData(vt, vp);
    m_plot->graph(1)->setData(vt, vd);
    m_plot->rescaleAxes();
    if(m_plot->xAxis->range().lower<=0) m_plot->xAxis->setRangeLower(1e-3);
    if(m_plot->yAxis->range().lower<=0) m_plot->yAxis->setRangeLower(1e-3);
    m_plot->replot();
}

void FittingWidget::on_comboModelSelect_currentIndexChanged(int) { on_btnResetParams_clicked(); }

QString FittingWidget::getParamDisplayName(const QString& key) {
    // ... (Same as before) ...
    if (key == "omega") return "储容比 (ω)";
    if (key == "omega1") return "内区储容比 (ω1)";
    if (key == "omega2") return "外区储容比 (ω2)";
    if (key == "lambda") return "窜流系数 (λ)";
    if (key == "lambda1") return "内区窜流系数 (λ1)";
    if (key == "lambda2") return "外区窜流系数 (λ2)";
    if (key == "S") return "表皮系数 (S)";
    if (key == "cD") return "井筒储存 (Cd)";
    if (key == "k") return "渗透率 (k)";
    if (key == "mf") return "裂缝条数 (Mf)";
    if (key == "mf1") return "内区裂缝数 (Mf1)";
    if (key == "mf2") return "外区裂缝数 (Mf2)";
    if (key == "nf") return "离散段数 (Nf)";
    if (key == "Xf") return "裂缝半长 (Xf)";
    if (key == "Xf1") return "内区半长 (Xf1)";
    if (key == "Xf2") return "外区半长 (Xf2)";
    if (key == "yy") return "裂缝间距 (Ye)";
    if (key == "yy1") return "内区间距 (Ye1)";
    if (key == "yy2") return "外区间距 (Ye2)";
    if (key == "y") return "水平井长 (L)";
    if (key == "CFD") return "导流能力 (Fcd)";
    if (key == "CFD1") return "内区导流 (Fcd1)";
    if (key == "CFD2") return "外区导流 (Fcd2)";
    if (key == "kpd") return "应力敏感系数 (γ)";
    if (key == "N") return "Stehfest N";
    return key;
}

void FittingWidget::on_btnResetParams_clicked() {
    if(!m_modelManager) return;

    QMap<QString,double> defs = m_modelManager->getDefaultParameters((ModelManager::ModelType)ui->comboModelSelect->currentIndex());

    m_parameters.clear();
    QMapIterator<QString,double> i(defs);
    while(i.hasNext()) {
        i.next();
        FitParameter p;
        p.name = i.key();
        p.displayName = getParamDisplayName(p.name);
        p.value = i.value();

        bool isReservoirParam = (
            p.name == "k" || p.name == "S" || p.name == "cD" ||
            p.name == "omega" || p.name == "omega1" || p.name == "omega2" ||
            p.name == "lambda" || p.name == "lambda1" || p.name == "lambda2" ||
            p.name == "CFD" || p.name == "CFD1" || p.name == "CFD2" ||
            p.name == "kpd"
            );
        p.isFit = isReservoirParam;

        if(p.value > 0) { p.min = p.value * 0.01; p.max = p.value * 100.0; }
        else { p.min = -10.0; p.max = 10.0; }

        m_parameters.append(p);
    }

    loadParamsToTable();

    QMap<QString,double> map;
    for(auto &p : m_parameters) map[p.name] = p.value;
    double err = calculateError(map, (ModelManager::ModelType)ui->comboModelSelect->currentIndex());
    onIterationUpdate(err, map);
}

void FittingWidget::loadParamsToTable() {
    // ... (Same as before) ...
    ui->tableParams->setRowCount(0);
    ui->tableParams->blockSignals(true);

    for(int i=0; i<m_parameters.size(); ++i) {
        ui->tableParams->insertRow(i);

        QTableWidgetItem* nameItem = new QTableWidgetItem(m_parameters[i].displayName);
        nameItem->setData(Qt::UserRole, m_parameters[i].name);
        nameItem->setFlags(nameItem->flags() ^ Qt::ItemIsEditable);
        ui->tableParams->setItem(i,0,nameItem);

        ui->tableParams->setItem(i,1,new QTableWidgetItem(QString::number(m_parameters[i].value)));

        QTableWidgetItem* chk = new QTableWidgetItem();
        chk->setFlags(Qt::ItemIsUserCheckable | Qt::ItemIsEnabled | Qt::ItemIsSelectable);
        chk->setCheckState(m_parameters[i].isFit ? Qt::Checked : Qt::Unchecked);
        chk->setTextAlignment(Qt::AlignCenter);
        ui->tableParams->setItem(i,2,chk);

        ui->tableParams->setItem(i,3,new QTableWidgetItem(QString("[%1,%2]").arg(m_parameters[i].min).arg(m_parameters[i].max)));
    }
    ui->tableParams->blockSignals(false);
}

void FittingWidget::updateParamsFromTable() {
    // ... (Same as before) ...
    for(int i=0; i<ui->tableParams->rowCount(); ++i) {
        if(i < m_parameters.size()) {
            QString key = ui->tableParams->item(i,0)->data(Qt::UserRole).toString();

            if(m_parameters[i].name == key) {
                m_parameters[i].value = ui->tableParams->item(i,1)->text().toDouble();
                m_parameters[i].isFit = (ui->tableParams->item(i,2)->checkState() == Qt::Checked);
            }
        }
    }
}

QStringList FittingWidget::parseLine(const QString& line) { return line.split(QRegularExpression("[,\\s\\t]+"), Qt::SkipEmptyParts); }

void FittingWidget::on_btnLoadData_clicked() {
    // ... (Same as before) ...
    QString path = QFileDialog::getOpenFileName(this, "加载试井数据", "", "文本文件 (*.txt *.csv)");
    if(path.isEmpty()) return;
    QFile f(path); if(!f.open(QIODevice::ReadOnly)) return;
    QTextStream in(&f); QList<QStringList> data;
    while(!in.atEnd()) { QString l=in.readLine().trimmed(); if(!l.isEmpty()) data<<parseLine(l); }
    f.close();

    FittingDataLoadDialog dlg(data, this);
    if(dlg.exec()!=QDialog::Accepted) return;

    int tCol=dlg.getTimeColumnIndex(), pCol=dlg.getPressureColumnIndex(), dCol=dlg.getDerivativeColumnIndex();
    QVector<double> t, p, d;
    double p_init = 0;

    if(pCol>=0) {
        for(int i=dlg.getSkipRows(); i<data.size(); ++i) {
            if(pCol<data[i].size()) { p_init=data[i][pCol].toDouble(); break; }
        }
    }

    for(int i=dlg.getSkipRows(); i<data.size(); ++i) {
        if(tCol<data[i].size()) {
            double tv = data[i][tCol].toDouble();
            double pv = (pCol>=0 && pCol<data[i].size()) ? std::abs(data[i][pCol].toDouble()-p_init) : 0;
            if(tv>0) { t<<tv; p<<pv; }
        }
    }

    if (dCol >= 0) {
        for(int i=dlg.getSkipRows(); i<data.size(); ++i) {
            if(tCol<data[i].size() && data[i][tCol].toDouble() > 0) {
                if(dCol<data[i].size()) d << data[i][dCol].toDouble();
                else d << 0;
            }
        }
    } else {
        double lSpacing = 0.15;
        d = PressureDerivativeCalculator::calculateBourdetDerivative(t, p, lSpacing);
    }

    setObservedData(t, p, d);
}

void FittingWidget::on_btnRunFit_clicked() {
    if(m_isFitting) return;
    if(m_obsTime.isEmpty()) { QMessageBox::warning(this,"错误","请先加载观测数据。"); return; }

    updateParamsFromTable();

    m_isFitting = true;
    m_stopRequested = false;
    ui->btnRunFit->setEnabled(false);

    int algIndex = ui->comboAlgorithm->currentIndex();
    ModelManager::ModelType modelType = (ModelManager::ModelType)ui->comboModelSelect->currentIndex();

    // [CRITICAL FIX] 制作参数副本传递给子线程，避免多线程同时访问 m_parameters 导致崩溃
    QList<FitParameter> paramsCopy = m_parameters;

    (void)QtConcurrent::run([this, algIndex, modelType, paramsCopy](){
        // 参数副本传入
        runOptimizationTask(algIndex, modelType, paramsCopy);
    });
}

void FittingWidget::runOptimizationTask(int algIndex, ModelManager::ModelType modelType, QList<FitParameter> fitParams) {
    if (algIndex == 0) {
        runLevenbergMarquardtOptimization(modelType, fitParams);
    } else {
        runNelderMeadOptimization(modelType, fitParams);
    }
}

void FittingWidget::on_btnStop_clicked() { m_stopRequested=true; }

void FittingWidget::on_btnUpdateModel_clicked() {
    if(!m_modelManager) return;
    QMap<QString,double> finalParams;
    for(const auto& p : m_parameters) finalParams.insert(p.name, p.value);
    emit fittingCompleted((ModelManager::ModelType)ui->comboModelSelect->currentIndex(), finalParams);
}

// ===========================================================================
// 优化算法 1: Levenberg-Marquardt (LM)
// ===========================================================================
void FittingWidget::runLevenbergMarquardtOptimization(ModelManager::ModelType modelType, QList<FitParameter> params) {
    // 整个函数内部使用局部变量 params，不触碰 m_parameters
    QVector<int> fitIndices;
    for(int i=0; i<params.size(); ++i) {
        if(params[i].isFit) fitIndices.append(i);
    }

    int nParams = fitIndices.size();
    if(nParams == 0) { QMetaObject::invokeMethod(this, "onFitFinished"); return; }

    double lambda = 0.01;
    double lambdaFactor = 10.0;
    int maxIter = 100;

    QMap<QString, double> currentParamMap;
    for(const auto& p : params) currentParamMap.insert(p.name, p.value);

    QVector<double> residuals = calculateResiduals(currentParamMap, modelType);
    double currentSSE = calculateSumSquaredError(residuals);

    emit sigIterationUpdated(currentSSE / (residuals.isEmpty()?1:residuals.size()), currentParamMap);

    for(int iter = 0; iter < maxIter; ++iter) {
        if(m_stopRequested) break;
        emit sigProgress(iter * 100 / maxIter);

        // 使用 params
        QVector<QVector<double>> J = computeJacobian(currentParamMap, residuals, fitIndices, modelType, params);
        int nResiduals = residuals.size();
        if(nResiduals == 0) break;

        // ... (Hessian 计算保持不变) ...
        QVector<QVector<double>> Hessian(nParams, QVector<double>(nParams, 0.0));
        QVector<double> gradient(nParams, 0.0);

        for(int k = 0; k < nResiduals; ++k) {
            double r_k = residuals[k];
            for(int i = 0; i < nParams; ++i) {
                double J_ki = J[k][i];
                gradient[i] += J_ki * r_k;
                for(int j = 0; j <= i; ++j) {
                    Hessian[i][j] += J_ki * J[k][j];
                }
            }
        }
        for(int i = 0; i < nParams; ++i) {
            for(int j = i + 1; j < nParams; ++j) Hessian[i][j] = Hessian[j][i];
        }

        bool stepAccepted = false;
        for(int innerLoop = 0; innerLoop < 5; ++innerLoop) {
            QVector<QVector<double>> A = Hessian;
            for(int i = 0; i < nParams; ++i) {
                A[i][i] *= (1.0 + lambda);
            }

            QVector<double> negGradient(nParams);
            for(int i=0; i<nParams; ++i) negGradient[i] = -gradient[i];

            QVector<double> delta = solveLinearSystem(A, negGradient);

            QMap<QString, double> trialParamMap = currentParamMap;
            for(int i = 0; i < nParams; ++i) {
                int paramIdx = fitIndices[i];
                // 使用 params (本地副本) 进行边界检查
                double newVal = params[paramIdx].value + delta[i];
                if(newVal < params[paramIdx].min) newVal = params[paramIdx].min;
                if(newVal > params[paramIdx].max) newVal = params[paramIdx].max;
                trialParamMap[params[paramIdx].name] = newVal;
            }

            QVector<double> newResiduals = calculateResiduals(trialParamMap, modelType);
            double newSSE = calculateSumSquaredError(newResiduals);

            if(newSSE < currentSSE) {
                lambda /= lambdaFactor;
                currentSSE = newSSE;
                residuals = newResiduals;
                currentParamMap = trialParamMap;
                // 更新本地 params，以便下一次迭代计算 Jacobian 时使用最新值（虽雅可比计算基于 Map，但为了一致性）
                for(int i = 0; i < nParams; ++i) {
                    params[fitIndices[i]].value = currentParamMap[params[fitIndices[i]].name];
                }
                stepAccepted = true;
                emit sigIterationUpdated(currentSSE / nResiduals, currentParamMap);
                break;
            } else {
                lambda *= lambdaFactor;
            }
        }

        if(!stepAccepted && lambda > 1e6) break;
        if(currentSSE < 1e-8) break;
    }
    // 拟合结束后，主线程的 onFitFinished 可以选择是否同步 m_parameters，这里主要依靠 sigIterationUpdated 实时更新 UI
    QMetaObject::invokeMethod(this, "onFitFinished");
}

// ===========================================================================
// 优化算法 2: Nelder-Mead (单纯形法)
// ===========================================================================
void FittingWidget::runNelderMeadOptimization(ModelManager::ModelType modelType, QList<FitParameter> params) {
    QVector<int> fitIndices;
    for(int i=0; i<params.size(); ++i) if(params[i].isFit) fitIndices.append(i);
    int n = fitIndices.size();
    if(n == 0) { QMetaObject::invokeMethod(this, "onFitFinished"); return; }

    double alpha = 1.0, beta = 0.5, gamma = 2.0;
    QVector<QVector<double>> simplex(n + 1);
    QVector<double> errors(n + 1);

    simplex[0].resize(n);
    for(int i=0; i<n; ++i) simplex[0][i] = params[fitIndices[i]].value;

    for(int i=1; i<=n; ++i) {
        simplex[i] = simplex[0];
        double val = simplex[i][i-1];
        simplex[i][i-1] = (std::abs(val) < 1e-9) ? 0.001 : val * 1.1;
    }

    auto getMap = [&](const QVector<double>& pVals) {
        QMap<QString,double> map;
        for(auto& p : params) map[p.name] = p.value;
        for(int k=0; k<n; ++k) map[params[fitIndices[k]].name] = pVals[k];
        return map;
    };

    for(int i=0; i<=n; ++i) {
        if(m_stopRequested) break;
        errors[i] = calculateError(getMap(simplex[i]), modelType);
    }

    int maxIter = 200;
    for(int iter=0; iter<maxIter; ++iter) {
        if(m_stopRequested) break;

        QVector<int> idx(n+1); std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](int a, int b){ return errors[a] < errors[b]; });
        int best = idx[0], worst = idx[n], secondWorst = idx[n-1];

        emit sigIterationUpdated(errors[best], getMap(simplex[best]));
        emit sigProgress(iter * 100 / maxIter);

        QVector<double> centroid(n, 0.0);
        for(int i=0; i<n; ++i) {
            int ii = idx[i];
            for(int j=0; j<n; ++j) centroid[j] += simplex[ii][j];
        }
        for(int j=0; j<n; ++j) centroid[j] /= n;

        QVector<double> xr(n);
        for(int j=0; j<n; ++j) xr[j] = centroid[j] + alpha * (centroid[j] - simplex[worst][j]);
        double errR = calculateError(getMap(xr), modelType);

        if(errR < errors[secondWorst] && errR >= errors[best]) {
            simplex[worst] = xr; errors[worst] = errR;
        } else if(errR < errors[best]) {
            QVector<double> xe(n);
            for(int j=0; j<n; ++j) xe[j] = centroid[j] + gamma * (xr[j] - centroid[j]);
            double errE = calculateError(getMap(xe), modelType);
            if(errE < errR) { simplex[worst] = xe; errors[worst] = errE; }
            else { simplex[worst] = xr; errors[worst] = errR; }
        } else {
            QVector<double> xc(n); bool accepted = false;
            if(errR < errors[worst]) {
                for(int j=0; j<n; ++j) xc[j] = centroid[j] + beta * (xr[j] - centroid[j]);
                double errC = calculateError(getMap(xc), modelType);
                if(errC <= errR) { simplex[worst] = xc; errors[worst] = errC; accepted = true;}
            } else {
                for(int j=0; j<n; ++j) xc[j] = centroid[j] - beta * (centroid[j] - simplex[worst][j]);
                double errC = calculateError(getMap(xc), modelType);
                if(errC < errors[worst]) { simplex[worst] = xc; errors[worst] = errC; accepted = true;}
            }
            if(!accepted) {
                for(int i=1; i<=n; ++i) {
                    int ii = idx[i];
                    for(int j=0; j<n; ++j) simplex[ii][j] = simplex[best][j] + 0.5 * (simplex[ii][j] - simplex[best][j]);
                    errors[ii] = calculateError(getMap(simplex[ii]), modelType);
                }
            }
        }
        // 更新 params 本地副本
        QMap<QString,double> bestMap = getMap(simplex[best]);
        for(int i=0; i<params.size(); ++i) {
            if(bestMap.contains(params[i].name)) params[i].value = bestMap[params[i].name];
        }
    }
    QMetaObject::invokeMethod(this, "onFitFinished");
}

QVector<double> FittingWidget::calculateResiduals(const QMap<QString, double>& params, ModelManager::ModelType modelType) {
    if(!m_modelManager) return QVector<double>();
    ModelCurveData res = m_modelManager->calculateTheoreticalCurve(modelType, params);

    const QVector<double>& pCal = std::get<1>(res);
    const QVector<double>& dpCal = std::get<2>(res);

    QVector<double> r;
    // 增加对空数据的检查，防止崩溃
    if (m_obsPressure.isEmpty() || pCal.isEmpty()) return r;

    int count = qMin(m_obsPressure.size(), pCal.size());
    r.reserve(count * 2);

    for(int i=0; i<count; ++i) {
        if(m_obsPressure[i] > 1e-6 && pCal[i] > 1e-6) {
            double diff = log(m_obsPressure[i]) - log(pCal[i]);
            r.append(diff);
        }
    }

    double derivWeight = sqrt(3.0);
    // 注意：obsDerivative 和 dpCal 可能长度不一致，或者比 pressure 短
    int dCount = qMin(m_obsDerivative.size(), dpCal.size());
    dCount = qMin(dCount, count);

    for(int i=0; i<dCount; ++i) {
        if(m_obsDerivative[i] > 1e-6 && dpCal[i] > 1e-6) {
            double diff = log(m_obsDerivative[i]) - log(dpCal[i]);
            r.append(diff * derivWeight);
        }
    }
    return r;
}

double FittingWidget::calculateSumSquaredError(const QVector<double>& residuals) {
    double sse = 0.0;
    for(double v : residuals) sse += v*v;
    return sse;
}

QVector<QVector<double>> FittingWidget::computeJacobian(const QMap<QString, double>& params, const QVector<double>& baseResiduals, const QVector<int>& fitIndices, ModelManager::ModelType modelType, const QList<FitParameter>& currentFitParams) {
    int nResiduals = baseResiduals.size();
    int nParams = fitIndices.size();

    QVector<QVector<double>> J(nResiduals, QVector<double>(nParams));
    double stepSize = 1e-4;

    for(int j = 0; j < nParams; ++j) {
        int paramIdx = fitIndices[j];
        // 使用 currentFitParams 获取正确的参数名
        QString pName = currentFitParams[paramIdx].name;
        double originalVal = params.value(pName); // 安全获取
        double h = std::abs(originalVal) * stepSize;
        if(h < 1e-6) h = 1e-6;

        QMap<QString, double> pPlus = params;
        pPlus[pName] = originalVal + h;
        QVector<double> rPlus = calculateResiduals(pPlus, modelType);

        QMap<QString, double> pMinus = params;
        pMinus[pName] = originalVal - h;
        QVector<double> rMinus = calculateResiduals(pMinus, modelType);

        if(rPlus.size() != nResiduals || rMinus.size() != nResiduals) {
            for(int i=0; i<nResiduals; ++i) J[i][j] = 0;
            continue;
        }

        for(int i = 0; i < nResiduals; ++i) {
            J[i][j] = (rPlus[i] - rMinus[i]) / (2.0 * h);
        }
    }
    return J;
}

// ... (solveLinearSystem, calculateError, onIterationUpdate, onFitFinished, plotCurves 保持不变) ...

QVector<double> FittingWidget::solveLinearSystem(const QVector<QVector<double>>& A, const QVector<double>& b) {
    int n = b.size();
    QVector<QVector<double>> M = A;
    QVector<double> x = b;

    for (int k = 0; k < n - 1; ++k) {
        int maxRow = k;
        double maxVal = std::abs(M[k][k]);
        for (int i = k + 1; i < n; ++i) {
            if (std::abs(M[i][k]) > maxVal) {
                maxVal = std::abs(M[i][k]);
                maxRow = i;
            }
        }

        if(maxRow != k) {
            std::swap(M[k], M[maxRow]);
            std::swap(x[k], x[maxRow]);
        }

        if (std::abs(M[k][k]) < 1e-12) continue;

        for (int i = k + 1; i < n; ++i) {
            double factor = M[i][k] / M[k][k];
            for (int j = k; j < n; ++j) {
                M[i][j] -= factor * M[k][j];
            }
            x[i] -= factor * x[k];
        }
    }

    QVector<double> res(n);
    for (int i = n - 1; i >= 0; --i) {
        if (std::abs(M[i][i]) < 1e-12) {
            res[i] = 0;
        } else {
            double sum = 0.0;
            for (int j = i + 1; j < n; ++j) {
                sum += M[i][j] * res[j];
            }
            res[i] = (x[i] - sum) / M[i][i];
        }
    }
    return res;
}

double FittingWidget::calculateError(const QMap<QString,double>& trialParams, ModelManager::ModelType modelType) {
    QVector<double> r = calculateResiduals(trialParams, modelType);
    if(r.isEmpty()) return 1e9;
    return calculateSumSquaredError(r) / r.size();
}

void FittingWidget::onIterationUpdate(double err, const QMap<QString,double>& p) {
    ui->label_Error->setText(QString("当前误差(MSE): %1").arg(err, 0, 'f', 6));

    // 这里是唯一更新 UI 和 m_parameters 的地方（主线程）
    ui->tableParams->blockSignals(true);
    for(int i=0; i<ui->tableParams->rowCount(); ++i) {
        QString key = ui->tableParams->item(i, 0)->data(Qt::UserRole).toString();
        if(p.contains(key)) {
            double val = p[key];
            ui->tableParams->item(i, 1)->setText(QString::number(val, 'g', 5));
            // 同步更新主线程的参数存储，以便下次迭代或者停止时保留结果
            if(i < m_parameters.size() && m_parameters[i].name == key) {
                m_parameters[i].value = val;
            }
        }
    }
    ui->tableParams->blockSignals(false);

    if(m_modelManager) {
        auto res = m_modelManager->calculateTheoreticalCurve((ModelManager::ModelType)ui->comboModelSelect->currentIndex(), p);
        plotCurves(std::get<0>(res), std::get<1>(res), std::get<2>(res), true);
    }
}

void FittingWidget::onFitFinished() {
    m_isFitting = false; ui->btnRunFit->setEnabled(true);
    QMessageBox::information(this, "完成", "拟合优化已完成。");
}

void FittingWidget::plotCurves(const QVector<double>& t, const QVector<double>& p, const QVector<double>& d, bool isModel) {
    QVector<double> vt, vp, vd;

    double cD = 1.0;
    // 使用 m_parameters 查找 cD，这在主线程是安全的
    for(const auto& param : m_parameters) {
        if(param.name == "cD") {
            cD = param.value;
            break;
        }
    }
    if(cD <= 0) cD = 1.0;

    for(int i=0; i<t.size(); ++i) {
        double t_val = isModel ? (t[i] / cD) : t[i];

        if(t_val>1e-6 && p[i]>1e-6) {
            vt<<t_val; vp<<p[i];
            if(i<d.size() && d[i]>1e-6) vd<<d[i]; else vd<<1e-10;
        }
    }

    if(isModel) {
        m_plot->graph(2)->setData(vt, vp);
        m_plot->graph(3)->setData(vt, vd);
        m_plot->replot();
    }
}
