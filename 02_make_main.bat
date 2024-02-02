@echo off
rem gcc 10.3.0 (tdm64) win10
g++ main.cpp -O3 -std=c++11 -Wall -Wextra -Wpedantic -pedantic -L./ -lorbit -o orbit.exe
pause