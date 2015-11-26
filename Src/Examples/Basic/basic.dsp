# Microsoft Developer Studio Project File - Name="basic" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 5.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Application" 0x0101

CFG=basic - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "basic.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "basic.mak" CFG="basic - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "basic - Win32 Release" (based on "Win32 (x86) Application")
!MESSAGE "basic - Win32 Debug" (based on "Win32 (x86) Application")
!MESSAGE 

# Begin Project
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
MTL=midl.exe
RSC=rc.exe

!IF  "$(CFG)" == "basic - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /G5 /W3 /GX /Ot /Oa /Ow /Oi /Op /Oy /Ob2 /I "..\..\Inc" /I "..\Shared" /I "c:\mssdk\include" /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /FD /c
# SUBTRACT CPP /Og
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /o NUL /win32
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /o NUL /win32
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib d3dim.lib ddraw.lib dxguid.lib winmm.lib /nologo /subsystem:windows /machine:I386 /out:"../bin/basic_example.exe" /libpath:"c:\mssdk\lib"

!ELSEIF  "$(CFG)" == "basic - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /Zi /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /G5 /W3 /GX /Zi /Od /Op /Oy /I "..\..\Inc" /I "..\Shared" /I "c:\mssdk\include" /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /FD /c
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /o NUL /win32
# ADD MTL /nologo /D "_DEBUG" /mktyplib203 /o NUL /win32
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:windows /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib d3dim.lib ddraw.lib dxguid.lib winmm.lib /nologo /subsystem:windows /debug /machine:I386 /pdbtype:sept /libpath:"c:\mssdk\lib"
# SUBTRACT LINK32 /incremental:no

!ENDIF 

# Begin Target

# Name "basic - Win32 Release"
# Name "basic - Win32 Debug"
# Begin Source File

SOURCE=..\..\Cpp\actuator.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\actuator_fv.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\actuator_tv.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\bar.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\basic.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\bspline.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\bsplinesegment.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\circle.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\collision.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\connector.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\constraint.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\constraint_manager.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\containerlist.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\controller.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\cspline.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\csplinesegment.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\curve.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\cyl.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\doubleeuler.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\dyna.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\dyna_system.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\ellipsoid.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\euler.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\flatsurface.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\force_drawable.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\force_drawer.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\geo.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\largematrix.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\largevector.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\line.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\linehinge.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\linesegment.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\list.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\m_integrator.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\matrix.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\multibar.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\multirope.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\orientation.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\pid.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\plc.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\pointvector.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\pris.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\ptc.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\ptp.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\pts.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\rope.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\rungekutta2.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\rungekutta4.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\sensor.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\sensor_angle_v.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\sensor_avelo_v.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\sensor_dist_v.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\sensor_velo_v.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\spring.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\supvec.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\surface.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\torquespring.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\usr_force_drawable.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\vector4.cpp
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\vtv.cpp
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=..\..\Cpp\wheel.cpp
# PROP Exclude_From_Build 1
# End Source File
# End Target
# End Project
