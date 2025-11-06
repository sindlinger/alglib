# ALGLIB – Wave Runtime (Via Única)

Arquitetura mínima e única para acelerar FFT/Operações espectrais em GPU no MT5 usando um serviço Windows + DLL.

- Pipe único: `\\.\pipe\alglib-wave_pipe` (sem “2”).
- Serviço: `alglib_service.exe` (sem admin; DACL permissiva para clientes não elevados).
- DLL: `alglib.dll` (consumida por Terminal/Agents via `MQL5\Libraries`).
- Artefatos de build: `alglib\dist-wave`.
- Runtime consolidado: `C:\mql5\wave`.
- Script único: `alglib\wave.ps1`.

## Pré‑requisitos
- Windows 10/11 x64, Visual Studio 2022 (CMake).
- CUDA runtime (vendor embutido em `alglib\vendor\alglib-gpu\runtime\win64`).
- MetaEditor: `C:\mql5\MetaEditor64.exe`.

## .env (opcional)
Arquivo: `alglib\.env`
```
WAVE_USER_PROFILE=C:\Users\SEUUSUARIO
WAVE_TERMINAL_ID=PASTA-HEX-DO-TERMINAL
WAVE_METAEDITOR=C:\mql5\MetaEditor64.exe
WAVE_ROOT=C:\mql5\wave
WAVE_PIPE=\\.\\pipe\\alglib-wave_pipe
```
- O `wave.ps1` usa estes valores como default.

## Fluxo (sem admin)
1) Build → dist-wave
```
PS> cd alglib
PS> .\wave.ps1 -Action Build
```
2) Consolidar → C:\mql5\wave
```
PS> .\wave.ps1 -Action Consolidate
```
3) Deploy → Terminal/Agents (inicia serviço)
```
PS> .\wave.ps1 -Action Deploy -UserProfile "C:\Users\SEUUSUARIO" -TerminalId "PASTA-HEX"
```
4) Compilar Scripts/Indicador no Terminal
```
PS> .\wave.ps1 -Action Compile -UserProfile "C:\Users\SEUUSUARIO" -TerminalId "PASTA-HEX"
```
5) Iniciar serviço (apenas)
```
PS> .\wave.ps1 -Action Start
```

## Convenções
- Sem fallbacks (apenas `dist-wave` → `C:\mql5\wave`).
- Pipe único e nominal.
- Serviço imprime versão fixa (hardcoded) e aceita clientes não-elevados.

## Estrutura (alglib/)
- `CMakeLists.txt` → outputs em `dist-wave`
- `src/alglib_gpu/` → serviço/DLL (CUDA/OpenCL)
- `vendor/alglib-gpu/runtime/win64` → runtimes CUDA
- `wave.ps1` → script único (build/consolidate/deploy/compile/start)
- `.env` → parâmetros da instância atual
- `dist-wave/` → artefatos (após build)

