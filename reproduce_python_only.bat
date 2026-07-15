@echo off
setlocal
cd /d "%~dp0"
if exist ".venv\Scripts\python.exe" (
  ".venv\Scripts\python.exe" run_all.py --clean-output
) else (
  python run_all.py --clean-output
)
exit /b %ERRORLEVEL%
