@echo off
rem gcc 10.3.0 (tdm64) win10
g++ orbit.cpp -O3 -std=c++11 -Wall -Wextra -Wpedantic -pedantic -DBUILD_LIB -shared -o orbit.dll
pause