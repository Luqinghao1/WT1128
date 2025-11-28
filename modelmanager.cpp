#include "modelmanager.h"
#include "modelwidget1.h"
#include "modelwidget2.h"
#include "modelwidget3.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QGroupBox>
#include <QDebug>
#include <cmath>
#include <QtMath>
#include <algorithm> // Added for std::max_element

// ==================================================================================
//  ModelManager 构造与初始化
// ==================================================================================

ModelManager::ModelManager(QWidget* parent)
    : QObject(parent)
    , m_mainWidget(nullptr)
    , m_modelTypeCombo(nullptr)
    , m_modelStack(nullptr)
    , m_modelWidget1(nullptr)
    , m_modelWidget2(nullptr)
    , m_modelWidget3(nullptr)
    , m_currentModelType(InfiniteConductive)
{
}

ModelManager::~ModelManager() {}

void ModelManager::initializeModels(QWidget* parentWidget)
{
    if (!parentWidget) return;
    createMainWidget();
    setupModelSelection();

    m_modelStack = new QStackedWidget(m_mainWidget);
    m_modelWidget1 = new ModelWidget1(m_modelStack);
    m_modelWidget2 = new ModelWidget2(m_modelStack);
    m_modelWidget3 = new ModelWidget3(m_modelStack);

    m_modelStack->addWidget(m_modelWidget1);
    m_modelStack->addWidget(m_modelWidget2);
    m_modelStack->addWidget(m_modelWidget3);

    m_mainWidget->layout()->addWidget(m_modelStack);
    connectModelSignals();
    switchToModel(InfiniteConductive);

    if (parentWidget->layout()) parentWidget->layout()->addWidget(m_mainWidget);
    else {
        QVBoxLayout* layout = new QVBoxLayout(parentWidget);
        layout->addWidget(m_mainWidget);
        parentWidget->setLayout(layout);
    }
}

void ModelManager::createMainWidget()
{
    m_mainWidget = new QWidget();
    QVBoxLayout* mainLayout = new QVBoxLayout(m_mainWidget);
    mainLayout->setContentsMargins(10, 5, 10, 10);
    mainLayout->setSpacing(0);
    m_mainWidget->setLayout(mainLayout);
}

void ModelManager::setupModelSelection()
{
    if (!m_mainWidget) return;
    QGroupBox* selectionGroup = new QGroupBox("模型类型选择", m_mainWidget);
    selectionGroup->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);
    QHBoxLayout* selectionLayout = new QHBoxLayout(selectionGroup);
    selectionLayout->setContentsMargins(9, 9, 9, 9);
    selectionLayout->setSpacing(6);

    QLabel* typeLabel = new QLabel("模型类型:", selectionGroup);
    typeLabel->setMinimumWidth(100);
    m_modelTypeCombo = new QComboBox(selectionGroup);
    m_modelTypeCombo->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
    m_modelTypeCombo->setMinimumWidth(200);

    m_modelTypeCombo->addItem(getModelTypeName(InfiniteConductive));
    m_modelTypeCombo->addItem(getModelTypeName(FiniteConductive));
    m_modelTypeCombo->addItem(getModelTypeName(SegmentedMultiCluster));

    m_modelTypeCombo->setStyleSheet("color: black;");
    typeLabel->setStyleSheet("color: black;");
    selectionGroup->setStyleSheet("QGroupBox { color: black; font-weight: bold; }");

    connect(m_modelTypeCombo, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &ModelManager::onModelTypeSelectionChanged);
    selectionLayout->addWidget(typeLabel);
    selectionLayout->addWidget(m_modelTypeCombo);

    QVBoxLayout* mainLayout = qobject_cast<QVBoxLayout*>(m_mainWidget->layout());
    if (mainLayout) {
        mainLayout->addWidget(selectionGroup);
        mainLayout->setStretchFactor(selectionGroup, 0);
    }
}

void ModelManager::connectModelSignals()
{
    if (m_modelWidget1) connect(m_modelWidget1, &ModelWidget1::calculationCompleted, this, &ModelManager::onModel1CalculationCompleted);
    if (m_modelWidget2) connect(m_modelWidget2, &ModelWidget2::calculationCompleted, this, &ModelManager::onModel2CalculationCompleted);
    if (m_modelWidget3) connect(m_modelWidget3, &ModelWidget3::calculationCompleted, this, &ModelManager::onModel3CalculationCompleted);
}

void ModelManager::switchToModel(ModelType modelType)
{
    if (!m_modelStack) return;
    if (modelType == m_currentModelType) return;
    ModelType old = m_currentModelType;
    m_currentModelType = modelType;
    m_modelStack->setCurrentIndex((int)modelType);
    if (m_modelTypeCombo) m_modelTypeCombo->setCurrentIndex((int)modelType);
    emit modelSwitched(modelType, old);
}

void ModelManager::onModelTypeSelectionChanged(int index) { switchToModel((ModelType)index); }

QString ModelManager::getModelTypeName(ModelType type)
{
    switch (type) {
    case InfiniteConductive: return "无限导流双重孔隙介质页岩油藏渗流模型";
    case FiniteConductive: return "有限导流双重孔隙介质页岩油藏渗流模型";
    case SegmentedMultiCluster: return "分段多簇压裂水平井双重孔隙介质页岩油藏渗流模型";
    default: return "未知模型";
    }
}

QStringList ModelManager::getAvailableModelTypes()
{
    return { getModelTypeName(InfiniteConductive), getModelTypeName(FiniteConductive), getModelTypeName(SegmentedMultiCluster) };
}

void ModelManager::onModel1CalculationCompleted(const QString &t, const QMap<QString, double> &r) { emit calculationCompleted(t, r); }
void ModelManager::onModel2CalculationCompleted(const QString &t, const QMap<QString, double> &r) { emit calculationCompleted(t, r); }
void ModelManager::onModel3CalculationCompleted(const QString &t, const QMap<QString, double> &r) { emit calculationCompleted(t, r); }


// ==================================================================================
//  计算引擎核心实现
// ==================================================================================

QMap<QString, double> ModelManager::getDefaultParameters(ModelType type)
{
    QMap<QString, double> p;

    // 公共基础参数
    p.insert("k", 1.0);
    p.insert("S", 0.1);
    p.insert("cD", 1e-8);
    p.insert("N", 4.0);

    if (type == InfiniteConductive) {
        // 模型1：无限导流
        p.insert("omega", 0.05);
        p.insert("lambda", 1e-2);
        p.insert("mf", 3.0);
        p.insert("nf", 5.0);
        p.insert("Xf", 40.0);
        p.insert("yy", 70.0);
        p.insert("y", 1000.0);
    }
    else if (type == FiniteConductive) {
        // 模型2：有限导流
        p.insert("omega", 0.0155);
        p.insert("lambda", 0.083);
        p.insert("mf", 3.0);
        p.insert("nf", 5.0);
        p.insert("Xf", 193.0);
        p.insert("yy", 295.0);
        p.insert("y", 2758.0);
        p.insert("CFD", 0.9);
        p.insert("kpd", 0.04);
        p.insert("S", 0.81);
        p.insert("cD", 8.08e-8);
    }
    else if (type == SegmentedMultiCluster) {
        // 模型3：分段多簇
        p.insert("omega1", 0.05);
        p.insert("omega2", 0.05);
        p.insert("lambda1", 1e-1);
        p.insert("lambda2", 1e-1);
        p.insert("mf1", 2.0);
        p.insert("mf2", 2.0);
        p.insert("nf", 5.0);
        p.insert("Xf1", 40.0);
        p.insert("Xf2", 40.0);
        p.insert("yy1", 70.0);
        p.insert("yy2", 70.0);
        p.insert("y", 800.0);
        p.insert("CFD1", 0.4);
        p.insert("CFD2", 0.4);
        p.insert("kpd", 0.045);
    }
    return p;
}

ModelCurveData ModelManager::calculateTheoreticalCurve(ModelType type, const QMap<QString, double>& params)
{
    switch(type) {
    case InfiniteConductive:
        return calculateModel1(params);
    case FiniteConductive:
        return calculateModel2(params);
    case SegmentedMultiCluster:
        return calculateModel3(params);
    default:
        return calculateModel1(params);
    }
}

// ==================================================================================
//  Model 1: 无限导流
// ==================================================================================
ModelCurveData ModelManager::calculateModel1(const QMap<QString, double>& params)
{
    // 使用 100 个点进行拟合计算，范围 1e-7 到 1e6
    int numPoints = 100;
    int N = (int)params.value("N", 4);
    if (N % 2 != 0) N = 4;

    QVector<double> tD;
    tD.reserve(numPoints);
    for (int i = 0; i < numPoints; ++i) {
        double exponent = -7.0 + 13.0 * i / (numPoints - 1);
        tD.append(pow(10.0, exponent));
    }

    QVector<double> pd(numPoints, 0.0);
    double ln2 = log(2.0);

    for (int k = 0; k < numPoints; ++k) {
        if (tD[k] <= 0) continue;
        for (int m = 1; m <= N; ++m) {
            double s = m * ln2 / tD[k];
            double L = flaplace1(s, params);
            double coeff = stefestCoefficient(m, N);
            pd[k] += coeff * ln2 * L / tD[k];
        }
    }

    QVector<double> dpVec, tdpVec;
    double cD = params.value("cD", 1e-8);
    computePressureDerivative(tD, pd, cD, dpVec, tdpVec);

    return std::make_tuple(tD, pd, dpVec);
}

// ==================================================================================
//  Model 2: 有限导流
// ==================================================================================
ModelCurveData ModelManager::calculateModel2(const QMap<QString, double>& params)
{
    // 保持原有逻辑，但注意使用安全的 value() 方法
    int numPoints = 100;
    int N = (int)params.value("N", 4);
    double kpd = params.value("kpd", 0.04);

    QVector<double> tD;
    for (int i = 0; i < numPoints; ++i) {
        double exponent = -11.0 + 16.0 * i / (numPoints - 1);
        tD.append(pow(10.0, exponent));
    }

    QVector<double> pd(numPoints, 0.0);
    double ln2 = log(2.0);

    for (int k = 0; k < numPoints; ++k) {
        for (int m = 1; m <= N; ++m) {
            double s = m * ln2 / tD[k];
            double L = flaplace2(s, params);
            pd[k] += stefestCoefficient(m, N) * ln2 * L / tD[k];
        }
    }

    QVector<double> pd1(numPoints);
    for(int i=0; i<numPoints; ++i) {
        pd1[i] = -1.0 / kpd * log(1.0 - kpd * pd[i]);
    }

    QVector<double> dpVec, tdpVec;
    double cD = params.value("cD", 1e-8);
    computePressureDerivative(tD, pd1, cD, dpVec, tdpVec);

    return std::make_tuple(tD, pd1, dpVec);
}

// ==================================================================================
//  Model 3: 分段多簇
// ==================================================================================
ModelCurveData ModelManager::calculateModel3(const QMap<QString, double>& params)
{
    // 目前回退到 Model1 接口，防止崩溃，等待完整实现
    return calculateModel1(params);
}

// ==================================================================================
//  通用导数计算 (加权中心差分)
// ==================================================================================
void ModelManager::computePressureDerivative(const QVector<double>& tD, const QVector<double>& pd,
                                             double cD, QVector<double>& dpd, QVector<double>& td_dpd)
{
    // 使用与 ModelWidget1 一致的更健壮的导数算法
    dpd.clear(); td_dpd.clear();
    if (tD.size() < 3) return;

    QVector<double> td;
    td.reserve(tD.size());
    for (double tD_val : tD) td.append(tD_val / cD);

    QVector<double> log_td, log_pd;
    for (int i = 0; i < td.size(); ++i) {
        if (td[i] > 0 && pd[i] > 0) {
            log_td.append(log(td[i]));
            log_pd.append(log(pd[i]));
        }
    }

    if (log_td.size() < 3) return;

    QVector<double> dlogP_dlogt;
    dlogP_dlogt.reserve(log_td.size());

    if (log_td.size() >= 2) {
        dlogP_dlogt.append((log_pd[1]-log_pd[0])/(log_td[1]-log_td[0]));
    }

    for (int i = 1; i < log_td.size() - 1; ++i) {
        double h1 = log_td[i] - log_td[i-1];
        double h2 = log_td[i+1] - log_td[i];
        double weight1 = h2 / (h1 + h2);
        double weight2 = h1 / (h1 + h2);
        double derivative = weight1 * ((log_pd[i]-log_pd[i-1])/h1) + weight2 * ((log_pd[i+1]-log_pd[i])/h2);
        dlogP_dlogt.append(derivative);
    }

    if (log_td.size() >= 2) {
        int n = log_td.size();
        dlogP_dlogt.append((log_pd[n-1]-log_pd[n-2])/(log_td[n-1]-log_td[n-2]));
    }

    // 简单平滑
    if (dlogP_dlogt.size() > 5) {
        QVector<double> smoothed = dlogP_dlogt;
        for (int i = 2; i < dlogP_dlogt.size() - 2; ++i) {
            smoothed[i] = dlogP_dlogt[i-2]*0.1 + dlogP_dlogt[i-1]*0.2 + dlogP_dlogt[i]*0.4 + dlogP_dlogt[i+1]*0.2 + dlogP_dlogt[i+2]*0.1;
        }
        dlogP_dlogt = smoothed;
    }

    for (int i = 1; i < td.size(); ++i) {
        if (i-1 < dlogP_dlogt.size()) {
            double dpd_val = pd[i] * dlogP_dlogt[i-1];
            if (dpd_val > 0 && std::isfinite(dpd_val)) {
                dpd.append(dpd_val);
                td_dpd.append(td[i]);
            }
        }
    }
}

// ==================================================================================
//  数学核心: 辅助函数与积分 (移植自 ModelWidget1)
// ==================================================================================

// 移植 ModelWidget1 的 8 点高斯积分，避免 0 点奇异性崩溃
double ModelManager::gaussQuadrature(double XDkv, double YDkv, double yDij, double fz, double a, double b)
{
    static const double points[] = {
        -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
        0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363
    };

    static const double weights[] = {
        0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
        0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763
    };

    if (std::abs(b - a) < 1e-15) return 0.0;

    double sum = 0.0;
    double halfWidth = (b - a) / 2.0;
    double center = (a + b) / 2.0;
    double sqrtFz = sqrt(std::abs(fz));

    for (int i = 0; i < 8; ++i) {
        double x = center + halfWidth * points[i];
        double distance = sqrt(pow(XDkv - x, 2) + pow(YDkv - yDij, 2));
        double arg = distance * sqrtFz;
        sum += weights[i] * besselK0(arg);
    }
    return sum * halfWidth;
}

double ModelManager::integralBesselK0(double XDkv, double YDkv, double yDij, double fz, double a, double b)
{
    // 直接使用高斯积分，ModelWidget1 的实现更直接
    return gaussQuadrature(XDkv, YDkv, yDij, fz, a, b);
}

// 移植 ModelWidget1 的 BesselK0 实现，处理更稳健
double ModelManager::besselK0(double x)
{
    if (x <= 0) return 1e10; // 避免 Inf
    if (x <= 2.0) {
        double t = x / 2.0;
        double t2 = t * t;
        double I0 = 1.0 + t2 * (1.0 + t2 * (0.25 + t2 * (1.0/36.0 + t2 * (1.0/576.0 + t2 * 1.0/14400.0))));
        double lnTerm = -log(t) * I0;
        double series = -0.5772156649015329 + t2 * (0.4227843350984671 + t2 * (0.2306975606171077 + t2 * 0.0348859015277786));
        return lnTerm + series;
    } else {
        double y = 2.0 / x;
        double y2 = y * y;
        double asymp = 1.0 + y2 * (-0.125 + y2 * (0.0703125 + y2 * (-0.0732421875 + y2 * 0.1121520996094)));
        return sqrt(M_PI / (2.0 * x)) * exp(-x) * asymp;
    }
}

QVector<double> ModelManager::solveLinearSystem(const QVector<QVector<double>>& A, const QVector<double>& b)
{
    // 增加对空矩阵或尺寸不匹配的检查
    if(A.isEmpty() || b.isEmpty() || A.size() != b.size()) return QVector<double>();

    int n = b.size();
    QVector<QVector<double>> Ab(n);
    for (int i = 0; i < n; ++i) {
        if(A[i].size() < n) return QVector<double>(); // 错误保护
        Ab[i] = A[i];
        Ab[i].append(b[i]);
    }

    for (int k = 0; k < n - 1; ++k) {
        int maxRow = k;
        double maxVal = std::abs(Ab[k][k]);
        for (int i = k + 1; i < n; ++i) {
            if (std::abs(Ab[i][k]) > maxVal) {
                maxVal = std::abs(Ab[i][k]);
                maxRow = i;
            }
        }
        if (maxRow != k) Ab[k].swap(Ab[maxRow]);
        if (std::abs(Ab[k][k]) < 1e-12) continue;

        for (int i = k + 1; i < n; ++i) {
            double f = Ab[i][k] / Ab[k][k];
            for (int j = k; j <= n; ++j) Ab[i][j] -= f * Ab[k][j];
        }
    }
    QVector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) sum += Ab[i][j] * x[j];
        x[i] = (std::abs(Ab[i][i]) < 1e-12) ? 0 : (Ab[i][n] - sum) / Ab[i][i];
    }
    return x;
}

double ModelManager::stefestCoefficient(int i, int N)
{
    // 移植 ModelWidget1 的实现
    double sum = 0.0;
    int start = (i + 1) / 2;
    int end = std::min(i, N / 2);
    for (int k = start; k <= end; ++k) {
        double num = pow(k, N/2.0) * factorial(2*k);
        double den = factorial(N/2-k) * factorial(k) * factorial(k-1) * factorial(i-k) * factorial(2*k-i);
        if(den > 0) sum += num/den;
    }
    return ((N/2+i)%2==0?1:-1)*sum;
}

double ModelManager::factorial(int n) {
    if (n<=1) return 1.0;
    static QVector<double> c={1.0,1.0};
    if (n < c.size()) return c[n];
    for(int i=c.size();i<=n;++i) c.append(c.last()*i);
    return c[n];
}

// ==================================================================================
//  Model 1 专用函数 (修复 QMap 访问崩溃问题)
// ==================================================================================
double ModelManager::flaplace1(double z, const QMap<QString, double>& p) {
    // 使用 value() 替代 operator[]，防止崩溃和编译错误
    int mf = (int)p.value("mf", 3);
    int nf = (int)p.value("nf", 5);
    double omega = p.value("omega", 0.05);
    double lambda = p.value("lambda", 1e-2);
    double Xf = p.value("Xf", 40.0);
    double yy = p.value("yy", 70.0);
    double y = p.value("y", 1000.0);
    double S = p.value("S", 0.1);
    double cD = p.value("cD", 1e-8);

    if (mf <= 0 || nf <= 0) return 0.0;

    int sz = mf * 2 * nf;
    QVector<QVector<double>> E(sz + 1, QVector<double>(sz + 1, 0.0));
    QVector<double> F(sz + 1, 0.0);

    for(int i=1; i<=mf; ++i) {
        for(int j=1; j<=2*nf; ++j) {
            for(int k=1; k<=mf; ++k) {
                for(int v=1; v<=2*nf; ++v) {
                    E[(i-1)*2*nf+j-1][(k-1)*2*nf+v-1] = e_function(z,i,j,k,v,mf,nf,omega,lambda,Xf,yy,y);
                }
            }
            E[(i-1)*2*nf+j-1][sz] = -1.0;
            E[sz][(i-1)*2*nf+j-1] = f_function(j,nf,Xf,y);
        }
    }
    F[sz] = 1.0/z;

    QVector<double> res = solveLinearSystem(E, F);
    if (res.size() <= sz) return 0.0; // 求解失败保护

    double pwd = z * res[sz] + S;
    return pwd / (z * (1.0 + z * cD * pwd));
}

double ModelManager::flaplace2(double z, const QMap<QString, double>& p) {
    // 安全读取参数
    int mf = (int)p.value("mf", 3);
    int nf = (int)p.value("nf", 5);
    double omega = p.value("omega", 0.0155);
    double lambda = p.value("lambda", 0.083);
    double Xf = p.value("Xf", 193.0);
    double yy = p.value("yy", 295.0);
    double y = p.value("y", 2758.0);
    double S = p.value("S", 0.81);
    double cD = p.value("cD", 8.08e-8);
    double CFD = p.value("CFD", 0.9);

    if(mf<=0 || nf<=0) return 0.0;

    int sz = mf * nf;
    double deltaL = Xf / (nf * y);

    QVector<QVector<double>> I(sz+mf+1, QVector<double>(sz+mf+1));
    QVector<double> F(sz+mf+1, 0);

    for(int i=1; i<=mf; ++i) {
        for(int j=1; j<=nf; ++j) {
            int r = (i-1)*nf + j - 1;
            for(int k=1; k<=mf; ++k) {
                for(int v=1; v<=nf; ++v) {
                    double val = integralBesselK0(
                        ((2*v-2*nf-1)/(2.0*nf)*Xf)/y, (yy+(y-2*yy)/(mf-1)*(k-1))/y, (yy+(y-2*yy)/(mf-1)*(i-1))/y,
                        z*(omega*z*(1-omega)+lambda)/(lambda+(1-omega)*z), ((j-nf-1)/(double)nf*Xf)/y, ((j-nf)/(double)nf*Xf)/y
                        );
                    if(i==k) val -= (j==v ? M_PI/(CFD*8)*deltaL*deltaL : M_PI/CFD*std::abs(j-v)*deltaL*deltaL);
                    I[r][(k-1)*nf + v - 1] = val;
                }
            }
            I[r][sz] = (2*M_PI/CFD) * (((j-nf)/(double)nf*Xf/y) - ((j-nf-1)/(double)nf*Xf/y));
            I[r][sz+mf] = -1.0;
        }
    }
    for(int i=1; i<=mf; ++i) {
        for(int j=(i-1)*nf; j<i*nf; ++j) I[sz+i-1][j] = deltaL;
        I[sz+i-1][sz+i-1] = -1.0;
    }
    for(int i=sz; i<sz+mf; ++i) I[sz+mf][i] = 1.0;

    F[sz+mf] = 1.0/z;
    QVector<double> res = solveLinearSystem(I, F);
    if(res.size() <= sz+mf) return 0.0;
    double pd1 = res[sz+mf];
    return (z * pd1 + S) / (z + z*z*cD * (z * pd1 + S));
}

double ModelManager::e_function(double z, int i, int j, int k, int v, int mf, int nf, double omega, double lambda, double Xf, double yy, double y)
{
    // 逻辑与 ModelWidget1 保持一致
    double fz = (z * (omega * z * (1 - omega) + lambda)) / (lambda + (1 - omega) * z);
    double xij = (j - nf - 1) / (double)nf * Xf;
    double xDij = xij / y;
    double xij1 = (j - nf) / (double)nf * Xf;
    double xDij1 = xij1 / y;
    double Xkv = (2*v - 2*nf - 1) / (double)(2*nf) * Xf;
    double XDkv = Xkv / y;
    double yij = yy + (y - 2*yy) / (mf - 1) * (i - 1);
    double yDij = yij / y;
    double Ykv = yy + (y - 2*yy) / (mf - 1) * (k - 1);
    double YDkv = Ykv / y;
    return integralBesselK0(XDkv, YDkv, yDij, fz, xDij, xDij1);
}

double ModelManager::f_function(int j, int nf, double Xf, double y)
{
    double xij = (j - nf - 1) / (double)nf * Xf;
    double xDij = xij / y;
    double xij1 = (j - nf) / (double)nf * Xf;
    double xDij1 = xij1 / y;
    return xDij1 - xDij;
}
