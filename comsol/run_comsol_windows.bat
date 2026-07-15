@echo off
setlocal
set "COMSOL_HOME=E:\COMSOL\COMSOL64\Multiphysics"
if not exist "%COMSOL_HOME%\bin\win64\comsolcompile.exe" exit /b 2
if not exist logs mkdir logs
pushd "%~dp0"
"%COMSOL_HOME%\bin\win64\comsolcompile.exe" "Yuan2020Equivalent3D.java" > "logs\compile.log" 2>&1
if not exist Yuan2020Equivalent3D.class (
  popd
  exit /b 3
)
"%COMSOL_HOME%\bin\win64\comsolbatch.exe" -inputfile Yuan2020Equivalent3D.class -batchlog "logs\solve.log" -batchlogout
set "RC=%ERRORLEVEL%"
popd
exit /b %RC%
