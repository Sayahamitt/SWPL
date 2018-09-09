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
set clargs=/EHsc /arch:AVX2 /O2 /Oi /openmp /source-charset:utf-8

if "%1"=="performance" (
	set exearg=256 256 0.25e-6 1e-6 10e-6 100e-6 obs.bin
	set execmd=%1.exe
) else if "%1"=="propagation" (
	set exearg=128 128 0.25e-6 1e-6 10e-6 100e-6 obs.bin
	set execmd=%1.exe
) else if "%1"=="mpi" (
	set clargs=%clargs% /I "C:\Program Files (x86)\Microsoft SDKs\MPI\Include" /I "%MSMPI_INC%\x64" /link /LIBPATH:"%MSMPI_LIB64%" msmpi.lib
	set execmd=mpiexec -n 4 %1.exe
	set exearg=
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
cl %1.cpp %clargs%
@echo off

if exist %1.exe (
	echo %execmd% %exearg%
	%execmd% %exearg%
	del %1.exe
)
if exist %1.obj (
	del %1.obj
)

if not defined PAUSEFLAG (
	pause
)

@cd /d %current_dir%