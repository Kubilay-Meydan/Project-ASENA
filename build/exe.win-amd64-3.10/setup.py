import sys
from cx_Freeze import setup, Executable

# Dependencies are automatically detected, but it might need fine tuning.
build_exe_options = {"packages": ["os", 'ete3', "Bio", "sys","matplotlib",'PyQt5.QtWidgets','PyQt5.QtCore','PyQt5.QtGui','IPython.display','random','subprocess','platform','src.main','ipywidgets']}

# GUI applications require a different base on Windows (the default is for a
# console application).
base = None
if sys.platform == "win32":
    base = "Win32GUI"

setup(  name = "GUI.py",
        version = "1.1",
        description = "ASENA v1.1",
        options = {"build_exe": build_exe_options},
        executables = [Executable("GUI.py", base=base)])