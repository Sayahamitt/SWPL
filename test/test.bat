@echo off
for /f "usebackq tokens=*" %%i in (`cd`) do @set current_dir=%%i
@cd /d %~dp0

rem ビルド時にユーザー定義の環境変数を使いたい場合にはこのファイルと同じディレクトリにenvpath.batを作成し、その中で定義する。
if exist envpath.bat (
	call envpath.bat
	)
if not defined VisualStudioVersion (
	call "C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\vcvarsall.bat" amd64
)

if "%1"=="performance" (
	set exearg=256 256 0.25e-6 1e-6 10e-6 100e-6 obs.bin
) else if "%1"=="propagation" (
	set exearg=128 128 0.25e-6 1e-6 10e-6 100e-6 obs.bin
) else (
	echo SET TEST ARGUMENT
	@cd /d %current_dir%
	@exit /B
)

if exist %1.exe (
	del %1.exe
)
if exist %1.obj (
	del %1.obj
)


@echo on
cl /EHsc /arch:AVX2 /O2 /Oi /openmp /source-charset:utf-8 %1.cpp
@echo off

if exist %1.exe (
	echo %1.exe %exearg%
	%1.exe %exearg%
	del %1.exe
)
if exist %1.obj (
	del %1.obj
)

if not defined PAUSEFLAG (
	pause
)

@cd /d %current_dir%