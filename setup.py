import sys
from cx_Freeze import setup, Executable

# Dependencies are automatically detected, but it might need fine tuning.
build_exe_options = {"packages": ["os", "Bio", "sys","matplotlib"]}

# GUI applications require a different base on Windows (the default is for a
# console application).
base = None
if sys.platform == "win32":
    base = "Win32GUI"

setup(  name = "main.py",
        version = "0.1",
        description = "Main KM",
        options = {"build_exe": build_exe_options},
        executables = [Executable("main.py", base=base)])