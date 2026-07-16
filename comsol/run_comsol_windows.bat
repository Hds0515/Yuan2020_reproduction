@echo off
setlocal
set "COMSOL_HOME=E:\COMSOL\COMSOL64\Multiphysics"
if not exist "%COMSOL_HOME%\bin\win64\comsolcompile.exe" exit /b 2
if not exist generated_v4 mkdir generated_v4
if not exist results_v4 mkdir results_v4
if not exist logs_v4 mkdir logs_v4
pushd "%~dp0"
"%COMSOL_HOME%\bin\win64\comsolcompile.exe" "Yuan2020Equivalent3D.java" > "logs_v4\stage1_compile.log" 2>&1
if exist Yuan2020Equivalent3D.class move /Y Yuan2020Equivalent3D.class generated_v4\Yuan2020Equivalent3D.class >nul
if not exist generated_v4\Yuan2020Equivalent3D.class (
  popd
  exit /b 3
)
"%COMSOL_HOME%\bin\win64\comsolbatch.exe" -inputfile generated_v4\Yuan2020Equivalent3D.class -batchlog "logs_v4\stage1_flow_0p1ms.log" -batchlogout
set "RC=%ERRORLEVEL%"
if not exist results_v4\flow_0p1ms_converged.mph set "RC=4"
popd
exit /b %RC%
