# ALGLIB – Wave Runtime (Via Única)

Arquitetura mínima e única para acelerar FFT/Operações espectrais em GPU no MT5 usando um serviço Windows. Não existe caminho direto de computação via DLL no MQL: o MQL fala com o serviço.

- Pipe único: `\\.\pipe\alglib-wave_pipe` (sem “2”).
- Serviço: `alglib_service.exe` (sem admin; DACL permissiva para clientes não elevados).
- Caminhos de integração do MQL (escolha 1):
  - Include gerado (`Include/pipe/<Preset>.mqh`) que fala direto com o serviço via Named Pipe; ou
  - DLL ponte `alglib_bridge.dll` (opcional), que apenas repassa as mensagens do MQL para o serviço (sem computar nada).
- `alglib.dll` é componente de runtime do serviço (GPU); não deve ser importada diretamente pelo MQL para computação.
- Artefatos de build: `alglib\dist-wave`.
- Runtime consolidado: `C:\mql5\MQL5\dist-wave` — o serviço e DLLs são sempre buscados em `dist-wave` dentro da pasta MQL5 (sem usar `C:\mql5\wave` e sem fallbacks para subpastas como `Release`).
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
2) Consolidar/Deploy → MQL5 (Libraries/Services) + Agents
```
PS> .\wave.ps1 -Action StageBins
```
3) Deploy → Terminal/Agents (sem fallbacks)
```
PS> .\wave.ps1 -Action Deploy -UserProfile "C:\Users\SEUUSUARIO" -TerminalId "PASTA-HEX"
```
4) Compilar Scripts/Indicador no Terminal
```
PS> .\wave.ps1 -Action Compile -UserProfile "C:\Users\SEUUSUARIO" -TerminalId "PASTA-HEX"
```
5) Iniciar serviço (apenas, sem fallback)
```
PS> .\wave.ps1 -Action Start
```

## Convenções
- Sem fallbacks: serviço deve ser iniciado apenas de `C:\mql5\MQL5\dist-wave\alglib_service.exe`.
- Pipe único e nominal: `\\.\\pipe\\alglib-wave_pipe`.
- Serviço imprime versão fixa (hardcoded) e aceita clientes não-elevados.
- Limpeza de binários duplicados: use `wave.ps1 -Action Clean` para mover executáveis legados para `C:\mql5\Lixeira\<timestamp>`.

## Estrutura (alglib/)
- `CMakeLists.txt` → outputs em `dist-wave`
- `src/alglib_gpu/` → serviço/DLL (CUDA/OpenCL)
- `vendor/alglib-gpu/runtime/win64` → runtimes CUDA
- `wave.ps1` → script único (clean/build/consolidate/deploy/compile/start/doctor)
- `docs/ALGLIB_BRIDGE.md` → documentação da DLL ponte (via #import, sem computação)
- `.env` → parâmetros da instância atual
- `dist-wave/` → artefatos (após build)
