@echo off
if not exist arc mkdir arc
if not exist myvers.bat (set vnumber=0) else (call myvers.bat)
set /a vnumber=vnumber+1
echo set vnumber=%vnumber% > myvers.bat
copy %1.tpl arc\%1_v%vnumber%.tpl
echo archived  to %1.tpl to arc\%1_v%vnumber%.tpl
