@echo off
REM Adjust COMSOL path if necessary.
set COMSOL="C:\Program Files\COMSOL\COMSOL64\Multiphysics\bin\win64\comsolbatch.exe"
%COMSOL% -inputfile Yuan2020Equivalent3D.java -outputfile comsol_build_log.txt
pause
