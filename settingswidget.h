#ifndef SETTINGSWIDGET_H
#define SETTINGSWIDGET_H

#include <QWidget>
#include <QSettings>
#include <QFileDialog>
#include <QMessageBox>
#include <QStandardPaths>
#include <QDir>
#include <QTimer>

namespace Ui {
class SettingsWidget;
}

/**
 * @brief 试井解释软件系统设置界面
 *
 * 提供系统级别的配置管理，包括文件路径和系统设置
 */
class SettingsWidget : public QWidget
{
    Q_OBJECT

public:
    explicit SettingsWidget(QWidget *parent = nullptr);
    ~SettingsWidget();

    // 获取当前设置值
    QString getCurrentDataPath() const;
    QString getCurrentReportPath() const;
    QString getCurrentBackupPath() const;
    int getAutoSaveInterval() const;
    bool isBackupEnabled() const;

signals:
    /**
     * @brief 系统设置改变信号
     */
    void systemSettingsChanged();

    /**
     * @brief 自动保存设置改变信号
     * @param interval 新的自动保存间隔（分钟）
     */
    void autoSaveIntervalChanged(int interval);

    /**
     * @brief 备份设置改变信号
     * @param enabled 是否启用备份
     */
    void backupSettingsChanged(bool enabled);

private slots:
    // 导航相关
    void on_navigationList_currentRowChanged(int currentRow);

    // 路径设置相关
    void on_browseDataPathButton_clicked();
    void on_browseReportPathButton_clicked();
    void on_browseBackupPathButton_clicked();

    // 系统设置相关
    void on_autoSaveSpinBox_valueChanged(int value);
    void on_backupEnabledCheckBox_toggled(bool checked);
    void on_maxBackupsSpinBox_valueChanged(int value);
    void on_cleanupOldLogsCheckBox_toggled(bool checked);

    // 应用设置
    void on_applyButton_clicked();
    void on_cancelButton_clicked();
    // void on_closeButton_clicked();

private:
    Ui::SettingsWidget *ui;
    QSettings *m_settings;

    // 初始化方法
    void initializeInterface();
    void setupConnections();
    void loadSettings();
    void applySettings();

    // 设置验证
    bool validatePaths();
    bool validateSettings();

    // 路径管理
    void createDirectoryIfNotExists(const QString &path);
    QString getDefaultPath(const QString &pathType);

    // 消息框样式设置
    void setupMessageBoxStyle(QMessageBox *msgBox);

    // 内部状态
    bool m_settingsChanged;
    QTimer *m_validationTimer;

    // 常量
    static const int DEFAULT_AUTO_SAVE_INTERVAL;
    static const int DEFAULT_MAX_BACKUPS;
};

#endif // SETTINGSWIDGET_H
