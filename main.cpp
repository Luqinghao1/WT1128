#include "mainwindow.h"
#include <QApplication >
#include <QStyleFactory>
#include <QMessageBox>
#include <QFileDialog>

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    // 设置全局样式，确保所有对话框和消息框的文本都显示为黑色
    QString styleSheet = R"(
        /* 全局黑色文字样式 */
        QLabel, QLineEdit, QComboBox, QPushButton, QToolButton,
        QTreeView, QHeaderView, QTableView, QTabBar, QRadioButton,
        QCheckBox, QGroupBox, QMenu, QMenuBar, QStatusBar,
        QListView, QListWidget, QTextEdit, QPlainTextEdit {
            color: black;
        }

        /* 消息框样式 */
        QMessageBox QLabel {
            color: black;
        }

        /* 文件对话框样式 */
        QFileDialog QLabel, QFileDialog QTreeView, QFileDialog QComboBox {
            color: black;
        }

        /* 数据编辑器样式 */
        DataEditorWidget, DataEditorWidget * {
            color: black;
        }

        QTableView {
            alternate-background-color: #f0f0f0;
            background-color: white;
            gridline-color: #d0d0d0;
        }

        QTableView::item {
            color: black;
        }

        QHeaderView::section {
            background-color: #e0e0e0;
            color: black;
            padding: 4px;
            border: 1px solid #c0c0c0;
        }

        QPushButton {
            background-color: #e0e0e0;
            border: 1px solid #c0c0c0;
            padding: 5px 15px;
            min-width: 80px;
        }

        QPushButton:hover {
            background-color: #d0d0d0;
        }

        QPushButton:pressed {
            background-color: #c0c0c0;
        }
    )";

    app.setStyleSheet(styleSheet);

    // 为了确保消息框和文件对话框的字体颜色，我们可以子类化QMessageBox和QFileDialog
    // 但这里我们使用更简单的方法，直接设置全局样式

    // 设置所有已存在的对话框类的调色板
    QPalette darkTextPalette;
    darkTextPalette.setColor(QPalette::WindowText, Qt::black);
    darkTextPalette.setColor(QPalette::Text, Qt::black);
    darkTextPalette.setColor(QPalette::ButtonText, Qt::black);

    QApplication::setPalette(darkTextPalette);

    MainWindow w;
    w.show();

    // 让qt程序进入事件循环，不返回
    return app.exec();
}
