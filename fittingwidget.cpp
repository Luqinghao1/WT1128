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

    // 隐藏优化算法选择框，因为现在只使用一种改进的算法
    if(ui->comboAlgorithm) {
        ui->comboAlgorithm->setVisible(false);
        // 如果有对应的标签也建议隐藏，或者将界面改为显示 "算法: Levenberg-Marquardt"
    }

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

    // 移除算法选择逻辑，直接使用改进后的 LM 算法
    ModelManager::ModelType modelType = (ModelManager::ModelType)ui->comboModelSelect->currentIndex();

    QList<FitParameter> paramsCopy = m_parameters;

    (void)QtConcurrent::run([this, modelType, paramsCopy](){
        // 参数副本传入
        runOptimizationTask(modelType, paramsCopy);
    });
}

void FittingWidget::runOptimizationTask(ModelManager::ModelType modelType, QList<FitParameter> fitParams) {
    // 统一调用改进后的 Levenberg-Marquardt
    runLevenbergMarquardtOptimization(modelType, fitParams);
}

void FittingWidget::on_btnStop_clicked() { m_stopRequested=true; }

void FittingWidget::on_btnUpdateModel_clicked() {
    if(!m_modelManager) return;
    QMap<QString,double> finalParams;
    for(const auto& p : m_parameters) finalParams.insert(p.name, p.value);
    emit fittingCompleted((ModelManager::ModelType)ui->comboModelSelect->currentIndex(), finalParams);
}

// ===========================================================================
// 改进的 Levenberg-Marquardt 优化算法 (模仿 lsqcurvefit)
// ===========================================================================
void FittingWidget::runLevenbergMarquardtOptimization(ModelManager::ModelType modelType, QList<FitParameter> params) {
    QVector<int> fitIndices;
    for(int i=0; i<params.size(); ++i) {
        if(params[i].isFit) fitIndices.append(i);
    }

    int nParams = fitIndices.size();
    if(nParams == 0) { QMetaObject::invokeMethod(this, "onFitFinished"); return; }

    // 初始化参数
    double lambda = 0.01;      // 初始阻尼因子
    double lambdaUp = 10.0;    // 阻尼增加倍数
    double lambdaDown = 0.1;   // 阻尼减少倍数
    int maxIter = 100;         // 最大迭代次数
    double epsilon = 1e-6;     // 收敛阈值

    // 初始状态
    QMap<QString, double> currentParamMap;
    for(const auto& p : params) currentParamMap.insert(p.name, p.value);

    QVector<double> residuals = calculateResiduals(currentParamMap, modelType);
    double currentSSE = calculateSumSquaredError(residuals);

    // 发送初始进度
    emit sigIterationUpdated(currentSSE / (residuals.isEmpty() ? 1 : residuals.size()), currentParamMap);

    bool updateJ = true; // 标记是否需要更新雅可比矩阵
    QVector<QVector<double>> J;
    QVector<QVector<double>> H;
    QVector<double> g;

    for(int iter = 0; iter < maxIter; ++iter) {
        if(m_stopRequested) break;
        emit sigProgress(iter * 100 / maxIter);

        int nResiduals = residuals.size();
        if(nResiduals == 0) break;

        // 1. 计算雅可比矩阵 J, Hessian 近似 H = J^T * J, 梯度 g = J^T * r
        if (updateJ) {
            J = computeJacobian(currentParamMap, residuals, fitIndices, modelType, params);
            H = QVector<QVector<double>>(nParams, QVector<double>(nParams, 0.0));
            g = QVector<double>(nParams, 0.0);

            for(int k = 0; k < nResiduals; ++k) {
                double r_k = residuals[k];
                for(int i = 0; i < nParams; ++i) {
                    double J_ki = J[k][i];
                    g[i] += J_ki * r_k; // g = J^T * r
                    for(int j = 0; j <= i; ++j) {
                        H[i][j] += J_ki * J[k][j];
                    }
                }
            }
            // 对称化 H
            for(int i = 0; i < nParams; ++i) {
                for(int j = i + 1; j < nParams; ++j) H[i][j] = H[j][i];
            }
        }

        // 2. 应用阻尼因子 (Marquardt)
        // 求解 (H + lambda * I) * delta = -g
        QVector<QVector<double>> H_aug = H;
        for(int i = 0; i < nParams; ++i) {
            // 使用 Marquardt 缩放: H[i][i] * (1 + lambda) 或 H[i][i] + lambda
            // MATLAB 风格通常使用对角线缩放
            if (std::abs(H_aug[i][i]) > 1e-9)
                H_aug[i][i] *= (1.0 + lambda);
            else
                H_aug[i][i] += lambda;
        }

        QVector<double> negGradient(nParams);
        for(int i=0; i<nParams; ++i) negGradient[i] = -g[i];

        QVector<double> delta = solveLinearSystem(H_aug, negGradient);

        // 3. 计算新参数 (带边界投影)
        QMap<QString, double> trialParamMap = currentParamMap;
        bool changed = false;

        for(int i = 0; i < nParams; ++i) {
            int paramIdx = fitIndices[i];
            double oldVal = params[paramIdx].value;
            double newVal = oldVal + delta[i];

            // 边界约束 (Projection)
            if(newVal < params[paramIdx].min) newVal = params[paramIdx].min;
            if(newVal > params[paramIdx].max) newVal = params[paramIdx].max;

            trialParamMap[params[paramIdx].name] = newVal;
        }

        // 4. 评估新位置
        QVector<double> newResiduals = calculateResiduals(trialParamMap, modelType);
        double newSSE = calculateSumSquaredError(newResiduals);

        // 5. 接受或拒绝步骤
        if(newSSE < currentSSE) {
            // 成功：减少阻尼，移动到新点
            lambda *= lambdaDown;
            currentSSE = newSSE;
            residuals = newResiduals;
            currentParamMap = trialParamMap;

            // 更新本地 params 副本值
            for(int i = 0; i < nParams; ++i) {
                params[fitIndices[i]].value = currentParamMap[params[fitIndices[i]].name];
            }

            updateJ = true; // 下次迭代需要重新计算 J

            emit sigIterationUpdated(currentSSE / nResiduals, currentParamMap);

            // 检查收敛 (误差变化极小 或 梯度极小)
            if (std::abs(currentSSE - newSSE) < epsilon) break;
        } else {
            // 失败：增加阻尼，不移动
            lambda *= lambdaUp;
            updateJ = false; // J 和 H 保持不变，只改变 lambda 重新求解
        }

        // 防止 lambda 过大溢出
        if (lambda > 1e9) break;
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
