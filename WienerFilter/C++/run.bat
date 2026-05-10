@echo off
g++ -Wall -std=c++11 main.cpp wienerFilter.cpp -o main
if %errorlevel% neq 0 (
echo Compile failed!
pause
exit /b
)

main
pause
