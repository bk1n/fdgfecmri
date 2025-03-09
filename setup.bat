@echo off
setlocal

:: Get the hostname
set HOSTNAME=%COMPUTERNAME%

echo Device name is %HOSTNAME%

set SCRIPT_DIR=%~dp0
cd /d %SCRIPT_DIR%

echo Working directory %cd%

:: Check the hostname and execute commands accordingly
if "%HOSTNAME%"=="BKPC" (
    echo Identified Home PC ^(BKPC^)
    mklink /D "data" "F:\Documents\OneDrive - University of Edinburgh\projects\MRes\data\" 
    mklink /D "outputs" "F:\Documents\OneDrive - University of Edinburgh\projects\MRes\outputs\"
) else if "%HOSTNAME%"=="MVM-IGC-D0060" (
    echo Identified Work PC ^(MVM-IGC-D0060^)
    mklink /D "data" "C:\Users\s1754085\OneDrive - University of Edinburgh\projects\MRes\data\" 
    mklink /D "outputs" "C:\Users\s1754085\OneDrive - University of Edinburgh\projects\MRes\outputs\"
) else if "%HOSTNAME%"=="BKLT" (
    echo Identified Laptop ^(BKLT^)
    mklink /D "data" "C:\Users\bgkin\OneDrive - University of Edinburgh\projects\MRes\data\" 
    mklink /D "outputs" "C:\Users\bgkin\OneDrive - University of Edinburgh\projects\MRes\outputs\"
) else (
    echo Device name not matched
    pause
)

timeout /t 10

endlocal