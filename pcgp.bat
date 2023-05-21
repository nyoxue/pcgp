@echo off
set N=32
set K=8
set C=1

set program="pcgp.exe"
set workdir="pcgp-%N%-%K%-%C%"
set arguments="s %N% %K% %C%"
mkdir "%workdir%"
start /WAIT /B /D "%workdir%" "" "%program%" "%arguments%"
pause
