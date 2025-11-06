param(
  [ValidateSet('Build','Consolidate','Deploy','Compile','Start','All')]
  [string]$Action = 'All',
  [string]$UserProfile,
  [string]$TerminalId,
  [string]$MetaEditor
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

# --- Constants / Paths (root = alglib/) ---
$AlglibRoot = Split-Path -Parent $PSCommandPath
$DistRoot   = Join-Path $AlglibRoot 'dist-wave'
$WaveRoot   = 'C:\\mql5\\wave'
$PipeName   = '\\.\pipe\alglib-wave_pipe'
$RuntimeSource = Join-Path $AlglibRoot 'vendor\alglib-gpu\runtime\win64'

$DllName     = 'alglib.dll'
$ServiceName = 'alglib_service.exe'

# --- .env (opcional) ---
$DotEnvPath = Join-Path $AlglibRoot '.env'
$DotEnv = @{}
if(Test-Path $DotEnvPath){
  Get-Content -Path $DotEnvPath | Where-Object { $_ -and ($_ -notmatch '^\s*#') } | ForEach-Object {
    if($_ -match '^(?<k>[^=\s]+)\s*=\s*(?<v>.*)$'){
      $DotEnv[$Matches.k] = $Matches.v
    }
  }
}

if(-not $UserProfile -and $DotEnv.ContainsKey('WAVE_USER_PROFILE')){ $UserProfile = $DotEnv['WAVE_USER_PROFILE'] }
if(-not $TerminalId  -and $DotEnv.ContainsKey('WAVE_TERMINAL_ID')) { $TerminalId  = $DotEnv['WAVE_TERMINAL_ID'] }
if(-not $MetaEditor  -and $DotEnv.ContainsKey('WAVE_METAEDITOR'))  { $MetaEditor  = $DotEnv['WAVE_METAEDITOR'] }
if($DotEnv.ContainsKey('WAVE_ROOT')){ $WaveRoot = $DotEnv['WAVE_ROOT'] }
if($DotEnv.ContainsKey('WAVE_PIPE')){ $PipeName = $DotEnv['WAVE_PIPE'] }

# --- Helpers ---
function Ensure-Folder([string]$Path) { if(-not (Test-Path $Path)) { New-Item -Path $Path -ItemType Directory -Force | Out-Null } }
function Copy-IfExists([string]$Source,[string]$Destination) { if(-not (Test-Path $Source)) { throw "Arquivo não encontrado: $Source" }; Ensure-Folder (Split-Path $Destination -Parent); Copy-Item -Path $Source -Destination $Destination -Force; Write-Host "[copy] $Source -> $Destination" }
function Stop-AlglibService { Get-Process -Name alglib_service -ErrorAction SilentlyContinue | Stop-Process -Force -ErrorAction SilentlyContinue }
function Start-AlglibService { param([string]$ExePath,[string]$Pipe=$PipeName) if(-not (Test-Path $ExePath)){throw "alglib_service não encontrado: $ExePath"}; Start-Process -FilePath $ExePath -ArgumentList "--pipe $Pipe" -ErrorAction SilentlyContinue | Out-Null }
function Get-TerminalRoot([string]$User,[string]$Tid){ if([string]::IsNullOrWhiteSpace($User) -or [string]::IsNullOrWhiteSpace($Tid)){ throw 'UserProfile/TerminalId obrigatórios' }; Join-Path $User ("AppData/Roaming/MetaQuotes/Terminal/{0}" -f $Tid) }

# --- Actions ---
function Do-Build {
  Write-Host '=== Build (dist-wave) ==='
  $buildDir = Join-Path $AlglibRoot 'build-win'
  Ensure-Folder $buildDir
  Push-Location $buildDir
  cmake -S .. -B . -G "Visual Studio 17 2022" -A x64
  if($LASTEXITCODE -ne 0){ throw 'cmake configure falhou' }
  cmake --build . --config Release
  if($LASTEXITCODE -ne 0){ throw 'cmake build falhou' }
  Pop-Location
  if(-not (Test-Path (Join-Path $DistRoot $ServiceName)) -or -not (Test-Path (Join-Path $DistRoot $DllName))) { throw "dist-wave não contém $ServiceName/$DllName. Verifique o build." }
}

function Do-Consolidate {
  Write-Host '=== Consolidate -> C:\mql5\wave ==='
  Stop-AlglibService
  Copy-IfExists (Join-Path $DistRoot $ServiceName) (Join-Path $WaveRoot $ServiceName)
  Copy-IfExists (Join-Path $DistRoot $DllName)     (Join-Path $WaveRoot $DllName)
  $rt = 'cudart64_13.dll','cufft64_12.dll','cufftw64_12.dll','cublas64_13.dll','cublasLt64_13.dll'
  foreach($n in $rt){ Copy-IfExists (Join-Path $RuntimeSource $n) (Join-Path $WaveRoot $n) }
}

function Do-Deploy {
  Write-Host '=== Deploy (Terminal/Agents) ==='
  if([string]::IsNullOrWhiteSpace($UserProfile) -or [string]::IsNullOrWhiteSpace($TerminalId)) { throw 'Informe -UserProfile e -TerminalId' }
  $termRoot = Join-Path (Get-TerminalRoot $UserProfile $TerminalId) 'MQL5'
  $lib     = Join-Path $termRoot 'Libraries'
  Ensure-Folder $lib
  Copy-IfExists (Join-Path $WaveRoot $DllName) (Join-Path $lib $DllName)
  $rt = 'cudart64_13.dll','cufft64_12.dll','cufftw64_12.dll','cublas64_13.dll','cublasLt64_13.dll'
  foreach($n in $rt){ Copy-IfExists (Join-Path $RuntimeSource $n) (Join-Path $lib $n) }
  # Agents
  $agents = Get-ChildItem -Path 'C:\\mql5\\Tester' -Directory -Filter 'Agent-*' -ErrorAction SilentlyContinue
  foreach($a in $agents){ $alib = Join-Path $a.FullName 'MQL5\Libraries'; Ensure-Folder $alib; Copy-IfExists (Join-Path $WaveRoot $DllName) (Join-Path $alib $DllName); foreach($n in $rt){ Copy-IfExists (Join-Path $RuntimeSource $n) (Join-Path $alib $n) } }
  # Serviço
  Stop-AlglibService
  Start-AlglibService -ExePath (Join-Path $WaveRoot $ServiceName) -Pipe $PipeName
}

function Do-Compile {
  Write-Host '=== Compile (Scripts/Indicator) ==='
  if([string]::IsNullOrWhiteSpace($UserProfile) -or [string]::IsNullOrWhiteSpace($TerminalId)) { throw 'Informe -UserProfile e -TerminalId' }
  $root = Join-Path (Get-TerminalRoot $UserProfile $TerminalId) 'MQL5'
  $outLog = 'C:\\mql5\\compile_alglib.log'
  function Compile-File([string]$Path){ if(-not (Test-Path $Path)){ Write-Warning "[compile] não encontrado: $Path"; return }; Write-Host "[compile] $Path"; & $MetaEditor "/compile:$Path" "/log:$outLog" | Out-Null }
  Compile-File (Join-Path $root 'Scripts\TestPipeCreateFile.mq5')
  Compile-File (Join-Path $root 'Scripts\TestAlglibPing.mq5')
  Compile-File (Join-Path $root 'Indicators\GPU_LegacyWave1.0.4.mq5')
}

function Do-Start { Write-Host '=== Start service (no-terminal) ==='; Stop-AlglibService; Start-AlglibService -ExePath (Join-Path $WaveRoot $ServiceName) -Pipe $PipeName }

switch($Action){
  'Build'       { Do-Build }
  'Consolidate' { Do-Consolidate }
  'Deploy'      { Do-Deploy }
  'Compile'     { Do-Compile }
  'Start'       { Do-Start }
  'All'         { Do-Build; Do-Consolidate; Do-Deploy; Do-Compile; Do-Start }
}

Write-Host '=== Done ==='
