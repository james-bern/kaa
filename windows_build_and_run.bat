@echo off
setlocal ENABLEDELAYEDEXPANSION

Taskkill /IM "executable.exe" /F  >nul 2>&1

for %%a in (%*) do set "argv[%%a]=1"

IF defined argv[--help] (

    echo build and   run in   debug mode: [36mwindows_build_and_run.bat[0m
    echo build and   run in release mode: [36mwindows_build_and_run.bat --release[0m
    echo build and debug in      VS Code: [36mwindows_build_and_run.bat --debug-vscode[0m
    echo build and debug in     remedyBG: [36mwindows_build_and_run.bat --debug-remedybg[0m

) ELSE (

    IF EXIST "*.dll" ( del *.dll )
    IF EXIST "*.exe" ( del *.exe )
    IF EXIST "*.exp" ( del *.exp )
    IF EXIST "*.ilk" ( del *.ilk )
    IF EXIST "*.lib" ( del *.lib )
    IF EXIST "*.obj" ( del main.obj )
    IF EXIST "*.pdb" ( del *.pdb )

    cls

    set O_SWITCH=
    IF defined argv[--release] (
        echo [36m[cow] compiling in release mode[0m
        set O_SWITCH=-O2
    ) ELSE (
        echo [36m[cow] compiling in debug mode[0m
        set O_SWITCH=-Od
    )

    set COMMON=main.cpp ^
    /openmp ^
    /I.\codebase\ext\ ^
    !O_SWITCH! ^
    /EHsc /MDd -fp:except -GR- -EHa- -FC -Zi ^
    -W4 -wd4201 -wd4127 ^
    /nologo ^
    /link /NODEFAULTLIB:MSVCRT /INCREMENTAL:NO ^
    OpenGL32.lib user32.lib gdi32.lib shell32.lib vcruntime.lib ^
    codebase\ext\windows_glfw3.lib 

    cl     /Femain.exe  !COMMON! 
    REM cl /DJIM_REDIRECT_PRINTF_TO_LOG_TXT /LD /Fesnake.dll !COMMON! 

    IF EXIST "main.exe" (
        IF defined argv[--debug-vscode] (
            echo [36m[cow] debugging in Visual Studio Code[0m
            _xplat_debug_vscode.bat
        ) ELSE IF defined argv[--debug-remedybg] (
            echo [36m[cow] debugging in remedyBG[0m
            call _windows_debug_remedybg.bat
        ) ELSE (
            echo [36m[cow] running executable[0m
            @echo on
            start main.exe
        )
    )

)


@echo off
endlocal

rem IF EXIST "*.exp" ( del *.exp )
rem IF EXIST "*.ilk" ( del *.ilk )
rem IF EXIST "*.lib" ( del *.lib )
rem IF EXIST "*.obj" ( del main.obj )

