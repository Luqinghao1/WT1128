#include "settingswidget.h"
#include "ui_settingswidget.h"
#include <QDebug>

// 常量定义
const int SettingsWidget::DEFAULT_AUTO_SAVE_INTERVAL = 10;
const int SettingsWidget::DEFAULT_MAX_BACKUPS = 10;

SettingsWidget::SettingsWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::SettingsWidget),
    m_settings(nullptr),
    m_settingsChanged(false),
    m_validationTimer(new QTimer(this))
{
    ui->setupUi(this);

    // 创建设置对象
    m_settings = new QSettings("WellTestPro", "WellTestAnalysis", this);

    // 设置验证定时器
    m_validationTimer->setSingleShot(true);
    m_validationTimer->setInterval(500);
    connect(m_validationTimer, &QTimer::timeout, this, &SettingsWidget::validateSettings);

    // 初始化界面
    initializeInterface();

    // 设置信号槽连接
    setupConnections();

    // 加载设置
    loadSettings();

    // 设置默认选中第一个导航项
    ui->navigationList->setCurrentRow(0);
    ui->contentStackedWidget->setCurrentIndex(0);
}

SettingsWidget::~SettingsWidget()
{
    delete ui;
}

void SettingsWidget::initializeInterface()
{
    // 设置窗口属性
    setWindowTitle("系统设置");
    setMinimumSize(800, 600);

    // 设置默认路径
    QString documentsPath = QStandardPaths::writableLocation(QStandardPaths::DocumentsLocation);
    ui->dataPathLineEdit->setText(documentsPath + "/WellTestPro/Data");
    ui->reportPathLineEdit->setText(documentsPath + "/WellTestPro/Reports");
    ui->backupPathLineEdit->setText(documentsPath + "/WellTestPro/Backups");
}

void SettingsWidget::setupConnections()
{
    // 导航相关连接
    connect(ui->navigationList, &QListWidget::currentRowChanged,
            this, &SettingsWidget::on_navigationList_currentRowChanged);

    // 路径设置连接
    connect(ui->browseDataPathButton, &QPushButton::clicked,
            this, &SettingsWidget::on_browseDataPathButton_clicked);
    connect(ui->browseReportPathButton, &QPushButton::clicked,
            this, &SettingsWidget::on_browseReportPathButton_clicked);
    connect(ui->browseBackupPathButton, &QPushButton::clicked,
            this, &SettingsWidget::on_browseBackupPathButton_clicked);

    // 系统设置连接
    connect(ui->autoSaveSpinBox, QOverload<int>::of(&QSpinBox::valueChanged),
            this, &SettingsWidget::on_autoSaveSpinBox_valueChanged);
    connect(ui->backupEnabledCheckBox, &QCheckBox::toggled,
            this, &SettingsWidget::on_backupEnabledCheckBox_toggled);
    connect(ui->maxBackupsSpinBox, QOverload<int>::of(&QSpinBox::valueChanged),
            this, &SettingsWidget::on_maxBackupsSpinBox_valueChanged);
    connect(ui->cleanupOldLogsCheckBox, &QCheckBox::toggled,
            this, &SettingsWidget::on_cleanupOldLogsCheckBox_toggled);

    // 应用取消按钮连接
    connect(ui->applyButton, &QPushButton::clicked,
            this, &SettingsWidget::on_applyButton_clicked);
    connect(ui->cancelButton, &QPushButton::clicked,
            this, &SettingsWidget::on_cancelButton_clicked);
    // connect(ui->closeButton, &QPushButton::clicked,
    //         this, &SettingsWidget::on_closeButton_clicked);

    // 监听设置变化
    QList<QLineEdit*> lineEdits = this->findChildren<QLineEdit*>();
    for (QLineEdit* lineEdit : lineEdits) {
        connect(lineEdit, &QLineEdit::textChanged, this, [this]() {
            m_settingsChanged = true;
        });
    }

    QList<QSpinBox*> spinBoxes = this->findChildren<QSpinBox*>();
    for (QSpinBox* spinBox : spinBoxes) {
        connect(spinBox, QOverload<int>::of(&QSpinBox::valueChanged), this, [this]() {
            m_settingsChanged = true;
        });
    }

    QList<QCheckBox*> checkBoxes = this->findChildren<QCheckBox*>();
    for (QCheckBox* checkBox : checkBoxes) {
        connect(checkBox, &QCheckBox::toggled, this, [this]() {
            m_settingsChanged = true;
        });
    }

    QList<QComboBox*> comboBoxes = this->findChildren<QComboBox*>();
    for (QComboBox* comboBox : comboBoxes) {
        connect(comboBox, QOverload<int>::of(&QComboBox::currentIndexChanged), this, [this]() {
            m_settingsChanged = true;
        });
    }
}

void SettingsWidget::loadSettings()
{
    // 加载路径设置
    QString documentsPath = QStandardPaths::writableLocation(QStandardPaths::DocumentsLocation);
    ui->dataPathLineEdit->setText(
        m_settings->value("paths/dataPath", documentsPath + "/WellTestPro/Data").toString());
    ui->reportPathLineEdit->setText(
        m_settings->value("paths/reportPath", documentsPath + "/WellTestPro/Reports").toString());
    ui->backupPathLineEdit->setText(
        m_settings->value("paths/backupPath", documentsPath + "/WellTestPro/Backups").toString());

    // 加载系统设置
    ui->autoSaveSpinBox->setValue(
        m_settings->value("system/autoSaveInterval", DEFAULT_AUTO_SAVE_INTERVAL).toInt());
    ui->backupEnabledCheckBox->setChecked(
        m_settings->value("system/backupEnabled", true).toBool());
    ui->maxBackupsSpinBox->setValue(
        m_settings->value("system/maxBackups", DEFAULT_MAX_BACKUPS).toInt());
    ui->cleanupOldLogsCheckBox->setChecked(
        m_settings->value("system/cleanupOldLogs", true).toBool());
    ui->logRetentionSpinBox->setValue(
        m_settings->value("system/logRetentionDays", 30).toInt());
    ui->logLevelComboBox->setCurrentIndex(
        m_settings->value("system/logLevel", 2).toInt());

    // 重置设置变化标志
    m_settingsChanged = false;
}

void SettingsWidget::applySettings()
{
    // 保存路径设置
    m_settings->setValue("paths/dataPath", ui->dataPathLineEdit->text());
    m_settings->setValue("paths/reportPath", ui->reportPathLineEdit->text());
    m_settings->setValue("paths/backupPath", ui->backupPathLineEdit->text());

    // 保存系统设置
    m_settings->setValue("system/autoSaveInterval", ui->autoSaveSpinBox->value());
    m_settings->setValue("system/backupEnabled", ui->backupEnabledCheckBox->isChecked());
    m_settings->setValue("system/maxBackups", ui->maxBackupsSpinBox->value());
    m_settings->setValue("system/cleanupOldLogs", ui->cleanupOldLogsCheckBox->isChecked());
    m_settings->setValue("system/logRetentionDays", ui->logRetentionSpinBox->value());
    m_settings->setValue("system/logLevel", ui->logLevelComboBox->currentIndex());

    // 同步设置到磁盘
    m_settings->sync();

    // 发送设置变更信号
    emit systemSettingsChanged();
    emit autoSaveIntervalChanged(ui->autoSaveSpinBox->value());
    emit backupSettingsChanged(ui->backupEnabledCheckBox->isChecked());

    // 重置设置变化标志
    m_settingsChanged = false;

    qDebug() << "系统设置已保存并应用";
}

bool SettingsWidget::validatePaths()
{
    QStringList paths = {
        ui->dataPathLineEdit->text(),
        ui->reportPathLineEdit->text(),
        ui->backupPathLineEdit->text()
    };

    for (const QString &path : paths) {
        if (path.trimmed().isEmpty()) {
            QMessageBox msgBox(this);
            msgBox.setWindowTitle("路径验证");
            msgBox.setText("所有路径都必须填写！");
            msgBox.setIcon(QMessageBox::Warning);
            setupMessageBoxStyle(&msgBox);
            msgBox.exec();
            return false;
        }

        QDir dir(path);
        if (!dir.exists()) {
            // 尝试创建目录
            if (!dir.mkpath(path)) {
                QMessageBox msgBox(this);
                msgBox.setWindowTitle("路径验证");
                msgBox.setText(QString("无法创建目录：%1\n请检查路径权限。").arg(path));
                msgBox.setIcon(QMessageBox::Warning);
                setupMessageBoxStyle(&msgBox);
                msgBox.exec();
                return false;
            }
        }
    }
    return true;
}

bool SettingsWidget::validateSettings()
{
    // 验证路径
    if (!validatePaths()) {
        return false;
    }

    return true;
}

void SettingsWidget::createDirectoryIfNotExists(const QString &path)
{
    QDir dir(path);
    if (!dir.exists()) {
        dir.mkpath(path);
        qDebug() << "创建目录：" << path;
    }
}

QString SettingsWidget::getDefaultPath(const QString &pathType)
{
    QString documentsPath = QStandardPaths::writableLocation(QStandardPaths::DocumentsLocation);
    return documentsPath + "/WellTestPro/" + pathType;
}

void SettingsWidget::setupMessageBoxStyle(QMessageBox *msgBox)
{
    if (!msgBox) return;

    // 设置消息框的样式，确保按钮文字清晰可见
    QString msgBoxStyle = R"(
        QMessageBox {
            background-color: #f8fafc;
            color: #1e293b;
            font-family: "Microsoft YaHei", "微软雅黑", sans-serif;
            font-size: 13px;
        }

        QMessageBox QLabel {
            color: #1e293b;
            font-size: 13px;
        }

        QMessageBox QPushButton {
            background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                        stop:0 #2563eb, stop:1 #1d4ed8);
            color: white;
            border: none;
            border-radius: 6px;
            padding: 8px 16px;
            font-weight: bold;
            font-size: 13px;
            min-width: 80px;
            min-height: 24px;
        }

        QMessageBox QPushButton:hover {
            background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                        stop:0 #3b82f6, stop:1 #2563eb);
        }

        QMessageBox QPushButton:pressed {
            background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                        stop:0 #1d4ed8, stop:1 #1e40af);
        }
    )";

    msgBox->setStyleSheet(msgBoxStyle);
}

// 槽函数实现
void SettingsWidget::on_navigationList_currentRowChanged(int currentRow)
{
    ui->contentStackedWidget->setCurrentIndex(currentRow);

    // 更新状态标签
    QStringList statusTexts = {
        "文件路径 - 配置数据和报告的存储位置",
        "系统设置 - 配置自动保存和日志管理"
    };

    if (currentRow >= 0 && currentRow < statusTexts.size()) {
        ui->statusLabel->setText(statusTexts[currentRow]);
    }
}

void SettingsWidget::on_browseDataPathButton_clicked()
{
    QString path = QFileDialog::getExistingDirectory(
        this, "选择数据文件存储路径", ui->dataPathLineEdit->text());
    if (!path.isEmpty()) {
        ui->dataPathLineEdit->setText(path);
    }
}

void SettingsWidget::on_browseReportPathButton_clicked()
{
    QString path = QFileDialog::getExistingDirectory(
        this, "选择报告输出路径", ui->reportPathLineEdit->text());
    if (!path.isEmpty()) {
        ui->reportPathLineEdit->setText(path);
    }
}

void SettingsWidget::on_browseBackupPathButton_clicked()
{
    QString path = QFileDialog::getExistingDirectory(
        this, "选择备份文件路径", ui->backupPathLineEdit->text());
    if (!path.isEmpty()) {
        ui->backupPathLineEdit->setText(path);
    }
}

void SettingsWidget::on_autoSaveSpinBox_valueChanged(int value)
{
    Q_UNUSED(value)
    m_validationTimer->start();
}

void SettingsWidget::on_backupEnabledCheckBox_toggled(bool checked)
{
    ui->maxBackupsSpinBox->setEnabled(checked);
}

void SettingsWidget::on_maxBackupsSpinBox_valueChanged(int value)
{
    Q_UNUSED(value)
    // 可以在这里添加额外的验证逻辑
}

void SettingsWidget::on_cleanupOldLogsCheckBox_toggled(bool checked)
{
    ui->logRetentionSpinBox->setEnabled(checked);
}

void SettingsWidget::on_applyButton_clicked()
{
    if (validateSettings()) {
        applySettings();

        // 创建自定义样式的消息框
        QMessageBox msgBox(this);
        msgBox.setWindowTitle("设置保存");
        msgBox.setText("设置已成功保存并应用！");
        msgBox.setIcon(QMessageBox::Information);
        setupMessageBoxStyle(&msgBox);
        msgBox.exec();

        // 不要关闭设置界面，让用户可以继续修改其他设置
        // this->close(); // 移除这行代码
    }
}

void SettingsWidget::on_cancelButton_clicked()
{
    if (m_settingsChanged) {
        QMessageBox msgBox(this);
        msgBox.setWindowTitle("确认重置");
        msgBox.setText("您有未保存的更改，确定要重置到原始设置吗？");
        msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
        msgBox.setDefaultButton(QMessageBox::No);
        msgBox.setIcon(QMessageBox::Question);
        setupMessageBoxStyle(&msgBox);

        QMessageBox::StandardButton reply = static_cast<QMessageBox::StandardButton>(msgBox.exec());

        if (reply == QMessageBox::Yes) {
            loadSettings(); // 重新加载设置，恢复到原始状态

            // 显示重置成功消息
            QMessageBox successBox(this);
            successBox.setWindowTitle("重置完成");
            successBox.setText("设置已重置为保存的状态！");
            successBox.setIcon(QMessageBox::Information);
            setupMessageBoxStyle(&successBox);
            successBox.exec();
        }
    } else {
        // 如果没有更改，直接重新加载设置
        loadSettings();

        QMessageBox infoBox(this);
        infoBox.setWindowTitle("重置完成");
        infoBox.setText("设置已重置为默认状态！");
        infoBox.setIcon(QMessageBox::Information);
        setupMessageBoxStyle(&infoBox);
        infoBox.exec();
    }
}

// void SettingsWidget::on_closeButton_clicked()
// {
//     if (m_settingsChanged) {
//         QMessageBox msgBox(this);
//         msgBox.setWindowTitle("确认关闭");
//         msgBox.setText("您有未保存的更改，确定要关闭设置窗口吗？\n未保存的更改将会丢失。");
//         msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
//         msgBox.setDefaultButton(QMessageBox::No);
//         msgBox.setIcon(QMessageBox::Question);
//         setupMessageBoxStyle(&msgBox);

//         QMessageBox::StandardButton reply = static_cast<QMessageBox::StandardButton>(msgBox.exec());

//         if (reply == QMessageBox::Yes) {
//             this->close();
//         }
//     } else {
//         this->close();
//     }
// }

// 获取当前设置值的方法
QString SettingsWidget::getCurrentDataPath() const
{
    return ui->dataPathLineEdit->text();
}

QString SettingsWidget::getCurrentReportPath() const
{
    return ui->reportPathLineEdit->text();
}

QString SettingsWidget::getCurrentBackupPath() const
{
    return ui->backupPathLineEdit->text();
}

int SettingsWidget::getAutoSaveInterval() const
{
    return ui->autoSaveSpinBox->value();
}

bool SettingsWidget::isBackupEnabled() const
{
    return ui->backupEnabledCheckBox->isChecked();
}
