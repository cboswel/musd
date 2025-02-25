@echo off

REM Check if the "venv" directory exists
if not exist "venv" (
    exit /b 1
  )

call .\venv\Scripts\activate
if errorlevel 1 (
  exit /b 1
)

python -m musd
