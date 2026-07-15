@echo off
setlocal
cd /d "%~dp0"
if exist ".venv\Scripts\python.exe" (
  ".venv\Scripts\python.exe" run_all.py --clean-output --with-comsol
) else (
  python run_all.py --clean-output --with-comsol
)
exit /b %ERRORLEVEL%
