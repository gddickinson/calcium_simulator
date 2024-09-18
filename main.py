#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 11:45:47 2024

@author: george
"""

import sys
from PyQt5.QtWidgets import QApplication
from calcium_model import CalciumModel
from gui import MainWindow

def main():
    app = QApplication(sys.argv)
    model = CalciumModel()
    window = MainWindow(model)
    window.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
