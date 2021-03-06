@echo off

if not defined _MACHINE set _MACHINE=ultra3d
if not defined _EXTRUDER set _EXTRUDER=left
if not defined _MATERIAL set _MATERIAL=inland_pla

set _INFILE=%1
set _OUTFILE=%_INFILE:.stl=_s.gcode%

shift
set params=%1
:loop
shift
if [%1]==[] goto afterloop
set params=%params% %1
goto loop
:afterloop

shiv.exe -c configs\global -c configs\\%_MACHINE%\\%_MACHINE% -c configs\\%_MACHINE%\\%_EXTRUDER% -c configs\\%_MACHINE%\\%_MATERIAL% -o %_OUTFILE% %params% %_INFILE%
