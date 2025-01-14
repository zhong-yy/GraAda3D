import os
import typing
import pandas as pd
import numpy as np
import shutil
from PyQt5 import QtCore, uic
from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    QWidget,
    QFileDialog,
    QLabel,
    QPushButton,
    QMessageBox,
    QGraphicsScene,
    QMenu,
    QGraphicsPixmapItem,
    QAction,
    QGraphicsView,
)
from PyQt5.QtGui import QImage, QPixmap, QPainter
from PyQt5.QtCore import Qt, QRectF
from pathlib import Path
from functools import partial
import sys


class MyGUI(QMainWindow):
    def __init__(
        self,
        ui_path: str | Path,
    ) -> None:
        super(MyGUI, self).__init__()
        # self.win = uic.loadUi(ui_path)
        uic.loadUi(str(ui_path), self)
        self.setWindowTitle("GraInv3D")
        self.select_data_file_btn.clicked.connect(self.select_data_file)
        self.datafilepath_lineEdit.setReadOnly(True)
        self.dataconfig_lineEdit.setReadOnly(True)

        self.data_config_path = Path("config_data")
        self.dataconfig_lineEdit.setText(str(self.data_config_path))
        self.datafile = None
        self.component_check_boxes = {
            "V": self.checkBox_V,
            "gx": self.checkBox_gx,
            "gy": self.checkBox_gy,
            "gz": self.checkBox_gz,
            "Txx": self.checkBox_Txx,
            "Txy": self.checkBox_Txy,
            "Txz": self.checkBox_Txz,
            "Tyy": self.checkBox_Tyy,
            "Tyz": self.checkBox_Tyz,
            "Tzz": self.checkBox_Tzz,
        }

        self.abs_err_line_edits = {
            "V": self.line_abs_err_V,
            "gz": self.line_abs_err_gz,
            "gx": self.line_abs_err_gx,
            "gy": self.line_abs_err_gy,
            "Tzz": self.line_abs_err_Tzz,
            "Txz": self.line_abs_err_Txz,
            "Tyz": self.line_abs_err_Tyz,
            "Txx": self.line_abs_err_Txx,
            "Txy": self.line_abs_err_Txy,
            "Tyy": self.line_abs_err_Tyy,
        }

        self.rel_err_line_edits = {
            "V": self.line_rel_err_V,
            "gz": self.line_rel_err_gz,
            "gx": self.line_rel_err_gx,
            "gy": self.line_rel_err_gy,
            "Tzz": self.line_rel_err_Tzz,
            "Txz": self.line_rel_err_Txz,
            "Tyz": self.line_rel_err_Tyz,
            "Txx": self.line_rel_err_Txx,
            "Txy": self.line_rel_err_Txy,
            "Tyy": self.line_rel_err_Tyy,
        }
        self.comp_order_line_edits = {
            "V": self.line_order_V,
            "gz": self.line_order_gz,
            "gx": self.line_order_gx,
            "gy": self.line_order_gy,
            "Tzz": self.line_order_Tzz,
            "Txz": self.line_order_Txz,
            "Tyz": self.line_order_Tyz,
            "Txx": self.line_order_Txx,
            "Txy": self.line_order_Txy,
            "Tyy": self.line_order_Tyy,
        }

        self.component_marker_dict = {
            "V": "0",
            "gz": "1",
            "gx": "2",
            "gy": "3",
            "Tzz": "4",
            "Txz": "5",
            "Tyz": "6",
            "Txx": "7",
            "Txy": "8",
            "Tyy": "9",
        }
        for comp, check_box in self.component_check_boxes.items():
            print(check_box.text())
            # check_box.stateChanged.connect(
            #     lambda _: self.check_component(check_box.text())
            # )
            check_box.stateChanged.connect(
                partial(self.check_component, comp=check_box.text())
            )
            self.rel_err_line_edits[comp].setReadOnly(True)
            self.abs_err_line_edits[comp].setReadOnly(True)
            self.comp_order_line_edits[comp].setReadOnly(True)
            self.rel_err_line_edits[comp].setStyleSheet("background-color: lightgray;")
            self.abs_err_line_edits[comp].setStyleSheet("background-color: lightgray;")
            self.comp_order_line_edits[comp].setStyleSheet(
                "background-color: lightgray;"
            )
        # self.checkBox_gx.stateChanged.connect(self.check_component)
        # self.checkBox_gy.stateChanged.connect(self.check_component)
        # self.checkBox_gz.stateChanged.connect(self.check_component)
        # self.checkBox_Txx.stateChanged.connect(self.check_component)
        # self.checkBox_Txy.stateChanged.connect(self.check_component)
        # self.checkBox_Txz.stateChanged.connect(self.check_component)
        # self.checkBox_Tyy.stateChanged.connect(self.check_component)
        # self.checkBox_Tyz.stateChanged.connect(self.check_component)
        # self.checkBox_Tzz.stateChanged.connect(self.check_component)
        # self.checkBox_V.stateChanged.connect(self.check_component)

        self.V_abs_err = 1
        self.V_rel_err = 0

        self.gz_abs_err = 1
        self.gz_rel_err = 0

        self.gx_abs_err = 1
        self.gx_rel_err = 0

        self.gy_abs_err = 1
        self.gy_rel_err = 0

        self.Tzz_abs_err = 1
        self.Tzz_rel_err = 0

        self.Txz_abs_err = 1
        self.Txz_rel_err = 0

        self.Tyz_abs_err = 1
        self.Tyz_rel_err = 0

        self.Txx_abs_err = 1
        self.Txx_rel_err = 0

        self.Txy_abs_err = 1
        self.Txy_rel_err = 0

        self.Tyy_abs_err = 1
        self.Tyy_rel_err = 0

        self.component_markers = []
        self.n_components = 0

        # self.current_

        self.set_saved_data_config_btn.clicked.connect(
            self.select_saving_path_for_data_configuration
        )

        self.generate_btn.clicked.connect(self.write_config)

    def select_file(self):
        dialog = QFileDialog(self)
        # dialog.setFileMode(dialog.FileMode.ExistingFiles)
        file_path, _ = dialog.getOpenFileName(self, "Select File", "", "All Files (*)")

        return file_path

    def select_data_file(self):
        file_path = self.select_file()
        if file_path:
            print(file_path)
            self.datafile = Path(file_path)
            self.datafilepath_lineEdit.setText(str(self.datafile))
            print(self.datafile)

    def select_save_path(self):
        initial_directory = os.getcwd()
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Select Output File", initial_directory, "All Files (*)"
        )
        return file_path

    def select_saving_path_for_data_configuration(self):
        file_path = self.select_save_path()
        self.data_config_path = Path(file_path)
        if file_path:
            self.dataconfig_lineEdit.setText(str(self.data_config_path))

    def check_component(self, comp):
        print(comp)
        checkbox = self.component_check_boxes[comp]
        if checkbox.isChecked():
            self.abs_err_line_edits[comp].setReadOnly(False)
            self.rel_err_line_edits[comp].setReadOnly(False)
            self.comp_order_line_edits[comp].setReadOnly(False)
            self.rel_err_line_edits[comp].setStyleSheet("background-color: white;")
            self.abs_err_line_edits[comp].setStyleSheet("background-color: white;")
            self.comp_order_line_edits[comp].setStyleSheet("background-color: white;")
            self.abs_err_line_edits[comp].setText("1")
            self.rel_err_line_edits[comp].setText("0.0")
            self.component_markers.append(self.component_marker_dict[comp])
            self.n_components += 1
            self.comp_order_line_edits[comp].setText(str(self.n_components))
        else:
            current_number = int(self.comp_order_line_edits[comp].text())
            for c in ["V", "gz", "gx", "gy", "Tzz", "Txz", "Tyz", "Txx", "Txy", "Tyy"]:
                if self.component_check_boxes[c].isChecked():
                    n = int(self.comp_order_line_edits[c].text())
                    if n > current_number:
                        self.comp_order_line_edits[c].setText(str(n - 1))
            self.abs_err_line_edits[comp].setText("")
            self.rel_err_line_edits[comp].setText("")
            self.comp_order_line_edits[comp].setText("")
            self.abs_err_line_edits[comp].setReadOnly(True)
            self.rel_err_line_edits[comp].setReadOnly(True)
            self.comp_order_line_edits[comp].setReadOnly(True)
            self.abs_err_line_edits[comp].setStyleSheet("background-color: lightgray;")
            self.rel_err_line_edits[comp].setStyleSheet("background-color: lightgray;")
            self.comp_order_line_edits[comp].setStyleSheet(
                "background-color: lightgray;"
            )
            self.component_markers.remove(self.component_marker_dict[comp])
            self.n_components -= 1

        print(self.n_components)
        print(self.component_markers)
        assert self.n_components == len(self.component_markers)

    # def check_component(self):
    #     self.component_markers = []
    #     for comp in ["V", "gz", "gx", "gy", "Tzz", "Txz", "Tyz", "Txx", "Txy", "Tyy"]:
    #         checkbox = self.component_check_boxes[comp]
    #         if checkbox.isChecked():
    #             self.component_markers.append(self.component_marker_dict[comp])
    #             self.abs_err_line_edits[comp].setReadOnly(False)
    #             self.rel_err_line_edits[comp].setReadOnly(False)
    #             self.rel_err_line_edits[comp].setStyleSheet("background-color: white;")
    #             self.abs_err_line_edits[comp].setStyleSheet("background-color: white;")
    #         else:
    #             self.abs_err_line_edits[comp].setReadOnly(True)
    #             self.rel_err_line_edits[comp].setReadOnly(True)
    #             self.abs_err_line_edits[comp].setStyleSheet(
    #                 "background-color: lightgray;"
    #             )
    #             self.rel_err_line_edits[comp].setStyleSheet(
    #                 "background-color: lightgray;"
    #             )
    #     self.n_components = len(self.component_markers)
    #     print(self.n_components, self.component_markers)

    def write_data_configuration(self):
        ready_to_write = self.check_data_config()
        if not ready_to_write:
            return

        with open(self.data_config_path, "w") as f:
            f.write("#data file name\n")
            f.write(str(self.datafile) + "\n\n")
            f.write(
                "#observation height z. If z axis is pointing downwards, this value should be less than 0 for observation points above the ground\n"
            )

            f.write(self.observation_height_lineEdit.text() + "\n\n")

            self.component_markers = []
            order = []
            errs = []
            for comp in [
                "V",
                "gz",
                "gx",
                "gy",
                "Tzz",
                "Txz",
                "Tyz",
                "Txx",
                "Txy",
                "Tyy",
            ]:
                checkbox = self.component_check_boxes[comp]
                if checkbox.isChecked():
                    self.component_markers.append(self.component_marker_dict[comp])
                    order.append(int(self.comp_order_line_edits[comp].text()) - 1)
                    errs.append(
                        self.rel_err_line_edits[comp].text().strip()
                        + " "
                        + self.abs_err_line_edits[comp].text().strip()
                    )
            print(order)
            print(self.component_markers)
            # self.component_markers = [self.component_markers[i] for i in order]
            # errs = [errs[i] for i in order]
            tmp1 = ["" for i in order]
            tmp2 = ["" for i in order]
            for i, marker, err in zip(order, self.component_markers, errs):
                tmp1[i] = marker
                tmp2[i] = err
            self.component_markers = tmp1
            errs = tmp2

            self.n_components = len(self.component_markers)
            f.write("#How many components will be use?\n")
            f.write(f"{self.n_components}\n\n")
            f.write("#gravity component markers\n")
            f.write(
                "#0:V   1:gz   2:gx   3:gy   4:T_zz   5:T_zx   6:T_zy   7:T_xx   8:T_xy   9:T_yy\n"
            )
            f.write(
                "Please make sure the order of the components here is consistent with the order of data columns in the data file\n"
            )
            f.write(" ".join(self.component_markers) + "\n\n")
            f.write("# Assumed relative error, assumed absolute error \n")
            f.write(
                "# The total error is assumed as (relative error*data+absolute error), which is used to construct the data weighting matrix Wd. The weight is 1/error. \n#If abs err is 1, rel err is 0, then there is no weight.\n"
            )
            f.write("\n".join(errs) + "\n")

    def check_data_config(self):
        ready = True
        if not self.datafile:
            QMessageBox.warning(
                self,
                "Warning",
                "Please choose a data file",
                QMessageBox.Ok,
            )
            ready = False
        if not self.observation_height_lineEdit.text().strip():
            print("The observation height is empty")
            QMessageBox.warning(
                self,
                "Warning",
                "The observation height is empty. Please set the observation height",
                QMessageBox.Ok,
            )
            ready = False
        if self.n_components == 0:
            QMessageBox.warning(
                self,
                "Warning",
                "Please select at least one component",
                QMessageBox.Ok,
            )
            ready = False
        return ready

    def write_config(self):
        self.write_data_configuration()


if __name__ == "__main__":
    pwd = Path(__file__).parent
    app = QApplication(sys.argv)
    mygui = MyGUI(pwd / "UI" / "mainwindow.ui")
    mygui.show()
    app.exec_()
