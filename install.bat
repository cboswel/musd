@echo off

REM Check if the "venv" directory exists
if not exist "venv" (
  python -m venv venv
  if errorlevel 1 (
    exit /b 1
  )
)

call .\venv\Scripts\activate
if errorlevel 1 (
  exit /b 1
)

pip install freegsnke/.
pip install .
pip install freegs4e
