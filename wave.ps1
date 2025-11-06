param(
  [ValidateSet('Build','Consolidate','Deploy','DeployMql','Compile','Start','Doctor','Undeploy','Generate','Connect','Clean','All')]
  [string]$Action = 'All',
  [string]$UserProfile,
  [string]$TerminalId,
  [string]$PortableRoot,
  [string]$MetaEditor = "C:\\mql5\\MetaEditor64.exe",
  [string]$LogPath,
  [string]$ReportPath,
  [string]$BackupRoot,
  [string]$BackupId,
  [string]$Preset = 'legacy_wave',
  [string]$IncludeName = 'GPU_LegacyWave1.0.4.mqh',
  [switch]$AutoFix,
  [switch]$TestService,
  [switch]$ZipBackup
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

# --- Constants / Paths (root = alglib/) ---
$AlglibRoot = Split-Path -Parent $PSCommandPath
$RepoRoot   = Split-Path -Parent $AlglibRoot
$DistRoot   = Join-Path $RepoRoot 'dist-wave'
$DistDir    = $DistRoot
$WaveRoot   = $DistRoot
$PipeName   = '\\.\pipe\alglib-wave_pipe'
$RuntimeSource = Join-Path $AlglibRoot 'vendor\alglib-gpu\runtime\win64'
if(-not $BackupRoot){ $BackupRoot = Join-Path $AlglibRoot 'backups' }
$TrashRoot  = Join-Path (Join-Path $env:SystemDrive 'mql5') 'Lixeira'

$DllName     = 'alglib.dll'
$ServiceName = 'alglib_service.exe'

# --- Logging ---
if(-not $LogPath){ $LogPath = Join-Path $AlglibRoot 'wave.log' }
function Write-Log([string]$Message,[string]$Level='INFO'){
  $ts = (Get-Date).ToString('yyyy-MM-dd HH:mm:ss.fff')
  $line = "[$ts][$Level] $Message"
  Add-Content -Path $LogPath -Value $line
  Write-Host $line
}

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

function Expand-Env([string]$s){ if([string]::IsNullOrWhiteSpace($s)){ return $s } else { return [Environment]::ExpandEnvironmentVariables($s) } }
if(-not $UserProfile -and $DotEnv.ContainsKey('WAVE_USER_PROFILE')){ $UserProfile = Expand-Env $DotEnv['WAVE_USER_PROFILE'] }
if(-not $TerminalId  -and $DotEnv.ContainsKey('WAVE_TERMINAL_ID')) { $TerminalId  = Expand-Env $DotEnv['WAVE_TERMINAL_ID'] }
if(-not $MetaEditor  -and $DotEnv.ContainsKey('WAVE_METAEDITOR'))  { $MetaEditor  = $DotEnv['WAVE_METAEDITOR'] }
if($DotEnv.ContainsKey('WAVE_ROOT')){ $WaveRoot = $DotEnv['WAVE_ROOT'] }
if($DotEnv.ContainsKey('WAVE_PIPE')){ $PipeName = $DotEnv['WAVE_PIPE'] }
if(-not $PortableRoot -and $DotEnv.ContainsKey('WAVE_PORTABLE_ROOT')){ $PortableRoot = $DotEnv['WAVE_PORTABLE_ROOT'] }
if($DotEnv.ContainsKey('WAVE_PRESET')){ $Preset = $DotEnv['WAVE_PRESET'] }
if($DotEnv.ContainsKey('WAVE_INCLUDE_NAME')){ $IncludeName = $DotEnv['WAVE_INCLUDE_NAME'] }

# --- Helpers ---
function Ensure-Folder([string]$Path) { if(-not (Test-Path $Path)) { New-Item -Path $Path -ItemType Directory -Force | Out-Null } }
function Copy-IfExists([string]$Source,[string]$Destination) {
  if(-not (Test-Path $Source)) { throw "Arquivo não encontrado: $Source" }
  Ensure-Folder (Split-Path $Destination -Parent)
  Copy-Item -Path $Source -Destination $Destination -Force
  Write-Log "copy $Source -> $Destination"
}
function Find-FirstFile([string[]]$Candidates){ foreach($c in $Candidates){ if(Test-Path $c){ return $c } } return $null }
function Get-PythonCmd {
  $cands = @('python3','py','python')
  foreach($c in $cands){ try { $cmd = (Get-Command $c -ErrorAction Stop).Source; return $cmd } catch {} }
  throw 'Python não encontrado (python3/py/python)'
}
function ConvertTo-WSLPath([string]$winPath){
  if([string]::IsNullOrWhiteSpace($winPath)){ return $winPath }
  $p = $winPath -replace '\\','/'
  if($p -match '^[A-Za-z]:/'){
    $drive = $p.Substring(0,1).ToLower()
    $rest = $p.Substring(2)
    return "/mnt/$drive$rest"
  }
  return $p
}
function Stop-AlglibService {
  $names = @('alglib_service','ALGLIB_Service','Gen4EngineService','Gen4Engine')
  foreach($n in $names){ Get-Process -Name $n -ErrorAction SilentlyContinue | Stop-Process -Force -ErrorAction SilentlyContinue }
}
function Start-AlglibService { param([string]$ExePath,[string]$Pipe=$PipeName) if(-not (Test-Path $ExePath)){throw "alglib_service não encontrado: $ExePath"}; Start-Process -FilePath $ExePath -ArgumentList "--pipe $Pipe" -ErrorAction SilentlyContinue | Out-Null }
function Get-TerminalRoot([string]$User,[string]$Tid){ if([string]::IsNullOrWhiteSpace($User) -or [string]::IsNullOrWhiteSpace($Tid)){ throw 'UserProfile/TerminalId obrigatórios' }; Join-Path $User ("AppData/Roaming/MetaQuotes/Terminal/{0}" -f $Tid) }
function Get-PortableRoot([string]$Root){ if([string]::IsNullOrWhiteSpace($Root)){ throw 'Informe -PortableRoot ou use TerminalId' }; return $Root }
# (auto-detecção de TerminalId removida para evitar parsing issues)

# --- Always pre-stop any running service before any action ---
Write-Host '=== Pre-Stop (kill legacy/new services) ==='
Stop-AlglibService
Start-Sleep -Milliseconds 150

# --- Clean legacy binaries (move to Lixeira) ---
function Get-TrashSessionPath {
  Ensure-Folder $TrashRoot
  $stamp = (Get-Date).ToString('yyyyMMdd-HHmmss')
  $dir = Join-Path $TrashRoot $stamp
  Ensure-Folder $dir
  return $dir
}
function Move-ToTrash([string]$File,[string]$SessionRoot){
  try{
    $drive = $File.Substring(0,1).ToUpper()
    $rel   = $File.Substring(2) -replace '^\\',''
    $dest  = Join-Path (Join-Path $SessionRoot $drive) $rel
    Ensure-Folder (Split-Path $dest -Parent)
    Move-Item -Path $File -Destination $dest -Force
    Write-Log ("trash {0} -> {1}" -f $File,$dest)
  } catch { Write-Log ("trash error {0}: {1}" -f $File, $_) 'ERROR' }
}
function Clean-OldServices {
  Write-Host '=== Clean (remove serviços duplicados) ==='
  $session = Get-TrashSessionPath
  $roots = @('C:\\mql5','C:\\mql5\\MQL5')
  if($UserProfile -and $TerminalId){ $roots += (Join-Path (Get-TerminalRoot $UserProfile $TerminalId) 'MQL5') }
  $keep = @((Join-Path $DistDir $ServiceName))
  $patterns = @('alglib_service.exe','ALGLIB_Service.exe','Gen4EngineService.exe','Gen4Engine.exe')
  foreach($r in $roots){
    if(-not (Test-Path $r)) { continue }
    Get-ChildItem -Path $r -Recurse -File -ErrorAction SilentlyContinue |
      Where-Object { $patterns -contains $_.Name } |
      ForEach-Object {
        $full = $_.FullName
        if($keep -contains $full){ return }
        Move-ToTrash -File $full -SessionRoot $session
      }
  }
}

# --- Backup helpers ---
function New-BackupSession {
  Ensure-Folder $BackupRoot
  $id = (Get-Date).ToString('yyyyMMdd_HHmmss')
  $dir = Join-Path $BackupRoot $id
  Ensure-Folder $dir
  return [pscustomobject]@{ Id=$id; Dir=$dir; Manifest=@() }
}
function Get-BackupPath([string]$Target,[string]$SessionDir){
  $full = [System.IO.Path]::GetFullPath($Target)
  $drive = [System.IO.Path]::GetPathRoot($full)
  $rel = $full.Substring($drive.Length)
  $driveTag = ($drive.TrimEnd('\')).Substring(0,1) + '_'
  $dest = Join-Path $SessionDir (Join-Path $driveTag $rel)
  return $dest
}
function Backup-File([string]$Target,[object]$Session){
  if(-not (Test-Path $Target)){ return }
  $dest = Get-BackupPath -Target $Target -SessionDir $Session.Dir
  Ensure-Folder (Split-Path $dest -Parent)
  Copy-Item -Path $Target -Destination $dest -Force
  $fi = Get-Item $Target
  $hash = (Get-FileHash -Algorithm SHA256 -Path $Target).Hash
  $Session.Manifest += [pscustomobject]@{ original=$Target; backup=$dest; size=$fi.Length; mtime=$fi.LastWriteTimeUtc.ToString('s'); sha256=$hash }
  Write-Log ("backup {0} -> {1}" -f $Target,$dest)
}
function Save-BackupManifest([object]$Session){
  $man = @{ id=$Session.Id; created=(Get-Date).ToString('s'); files=$Session.Manifest } | ConvertTo-Json -Depth 5
  $path = Join-Path $Session.Dir 'manifest.json'
  $man | Out-File -FilePath $path -Encoding UTF8
  Write-Log ("backup manifest: {0}" -f $path)
  if($ZipBackup){
    $zipPath = Join-Path $BackupRoot ("backup_{0}.zip" -f $Session.Id)
    if(Test-Path $zipPath){ Remove-Item -Path $zipPath -Force -ErrorAction SilentlyContinue }
    Compress-Archive -Path (Join-Path $Session.Dir '*') -DestinationPath $zipPath -Force
    Write-Log ("backup zip: {0}" -f $zipPath)
  }
}

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
  # Build concluído; artefatos devem estar em dist-wave (sem subpastas)
  $script:DistDir = $DistRoot
  $script:WaveRoot = $DistRoot
  Write-Log ("build OK em {0}" -f $script:DistDir)
}

function Do-Consolidate {
  Write-Host '=== Consolidate -> MQL5 (Libraries/Services) + Agents ==='
  Do-Deploy
}

function Do-Deploy {
  Write-Host '=== Deploy (Terminal/Agents) ==='
  Do-StageBins
  if($PortableRoot){ $termRoot = Join-Path (Get-PortableRoot $PortableRoot) 'MQL5' }
  else{
    if([string]::IsNullOrWhiteSpace($UserProfile)) { $UserProfile = [Environment]::GetFolderPath('UserProfile') }
    if([string]::IsNullOrWhiteSpace($TerminalId)) { throw 'Informe -TerminalId' }
    $termRoot = Join-Path (Get-TerminalRoot $UserProfile $TerminalId) 'MQL5'
  }
  $incPipe = Join-Path (Join-Path $termRoot 'Include') 'pipe'
  Ensure-Folder $incPipe
  # MQL-only: includes gerados + indicador
  $repoIncPipe = Join-Path (Join-Path $RepoRoot 'Include') 'pipe'
  foreach($n in @('Connect.mqh', $IncludeName)){
    $src = Join-Path $repoIncPipe $n
    if(Test-Path $src){ Copy-IfExists $src (Join-Path $incPipe $n) }
  }
  $repoInd = Join-Path $RepoRoot 'MQL5\Indicators\GPU_LegacyWave1.0.4.mq5'
  if(Test-Path $repoInd){ Copy-IfExists $repoInd (Join-Path (Join-Path $termRoot 'Indicators') 'GPU_LegacyWave1.0.4.mq5') }
  Write-Log "deploy OK (MQL-only)"
}

function Do-DeployMql {
  Write-Host '=== Deploy MQL (Indicator + Includes) ==='
  if($PortableRoot){ $termRoot = Join-Path (Get-PortableRoot $PortableRoot) 'MQL5' }
  else{
    if([string]::IsNullOrWhiteSpace($UserProfile)) { $UserProfile = [Environment]::GetFolderPath('UserProfile') }
    if([string]::IsNullOrWhiteSpace($TerminalId)) { throw 'Informe -TerminalId' }
    $termRoot = Join-Path (Get-TerminalRoot $UserProfile $TerminalId) 'MQL5'
  }
  $inc     = Join-Path $termRoot 'Include'
  $incPipe = Join-Path $inc 'pipe'
  $indDir  = Join-Path $termRoot 'Indicators'
  Ensure-Folder $incPipe
  Ensure-Folder $indDir
  # Copiar apenas os includes necessários
  $repoIncPipe = Join-Path (Join-Path $RepoRoot 'Include') 'pipe'
  foreach($n in @('Connect.mqh', $IncludeName)){
    $src = Join-Path $repoIncPipe $n
    if(Test-Path $src){ Copy-IfExists $src (Join-Path $incPipe $n) }
  }
  # Indicador
  $repoInd = Join-Path $RepoRoot 'MQL5\Indicators\GPU_LegacyWave1.0.4.mq5'
  if(Test-Path $repoInd){ Copy-IfExists $repoInd (Join-Path $indDir 'GPU_LegacyWave1.0.4.mq5') }
  Write-Log "deploy-mql OK"
}

function Do-Generate {
  Write-Host '=== Generate include from YAML preset ==='
  $ok = $false
  try{
    $py = Get-PythonCmd
    Push-Location $RepoRoot
    & $py 'tools/alglib_gen.py' 'config/alglib_ops.yaml' '--preset' $Preset '--non-interactive'
    if($LASTEXITCODE -ne 0){ throw 'py failed' }
    $ok = $true
    Pop-Location
  } catch {
    # Fallback: WSL python3
    $wsl = (Get-Command 'wsl.exe' -ErrorAction SilentlyContinue).Source
    if(-not $wsl){ throw 'Python não encontrado e WSL indisponível' }
    $repoWSL = ConvertTo-WSLPath $RepoRoot
    & $wsl -e bash -lc "cd '$repoWSL' && python3 tools/alglib_gen.py config/alglib_ops.yaml --preset '$Preset' --non-interactive"
    if($LASTEXITCODE -ne 0){ throw 'Falha ao gerar include via WSL/python3' }
    $ok = $true
  }
  if(-not $ok){ throw 'Falha ao gerar include via alglib_gen.py' }
  Write-Log ("generate OK preset={0}" -f $Preset)
}

function Do-Connect {
  Write-Host '=== Connect (service + include + compile) ==='
  # 1) Gera include atualizado
  Do-Generate
  # 2) Copia binários + includes + indicador
  Do-Deploy
  # 3) Start service (garante pipe correto)
  Stop-AlglibService
  if($PortableRoot){ $termRoot = Join-Path (Get-PortableRoot $PortableRoot) 'MQL5' } else { if([string]::IsNullOrWhiteSpace($UserProfile)) { $UserProfile = [Environment]::GetFolderPath('UserProfile') }; if([string]::IsNullOrWhiteSpace($TerminalId)) { throw 'Informe -TerminalId' }; $termRoot = Join-Path (Get-TerminalRoot $UserProfile $TerminalId) 'MQL5' }
  $svcPath = Join-Path (Join-Path $termRoot 'Services') $ServiceName
  Start-AlglibService -ExePath $svcPath -Pipe $PipeName
  Start-Sleep -Milliseconds 300
  # 4) Diagnóstico: tenta conectar ao pipe
  try{
    $client = New-Object System.IO.Pipes.NamedPipeClientStream('.', 'alglib-wave_pipe', [System.IO.Pipes.PipeDirection]::InOut, [System.IO.Pipes.PipeOptions]::Asynchronous)
    $client.Connect(1000)
    $ok = $client.IsConnected
    $client.Dispose()
    Write-Log ("pipe-connect ok={0}" -f $ok)
  } catch { Write-Log ("pipe-connect error: $_") 'ERROR' }
  # 5) Compila o indicador no Terminal
  Do-Compile
}

# Stage binários em MQL5\Libraries e MQL5\Services e espelhar alglib/ para C:\mql5\alglib

function Do-Compile {
  Write-Host '=== Compile (Scripts/Indicator) ==='
  if($PortableRoot){ $root = Join-Path (Get-PortableRoot $PortableRoot) 'MQL5' }
  else{
    if([string]::IsNullOrWhiteSpace($UserProfile)) { $UserProfile = [Environment]::GetFolderPath('UserProfile') }
    if([string]::IsNullOrWhiteSpace($TerminalId)) { throw 'Informe -TerminalId' }
    $root = Join-Path (Get-TerminalRoot $UserProfile $TerminalId) 'MQL5'
  }
  $outLog = 'C:\\mql5\\compile_alglib.log'
  function Compile-File([string]$Path){ if(-not (Test-Path $Path)){ Write-Warning "[compile] não encontrado: $Path"; return }; Write-Host "[compile] $Path"; & $MetaEditor "/compile:$Path" "/log:$outLog" | Out-Null }
  Compile-File (Join-Path $root 'Scripts\TestPipeCreateFile.mq5')
  Compile-File (Join-Path $root 'Scripts\TestAlglibPing.mq5')
  # Indicadores desabilitados por padrão; compile manualmente se necessário
  Write-Log "compile OK (scripts) em $root"
}

function Do-Start { Write-Host '=== Start service (no-terminal) ==='; Stop-AlglibService; Start-AlglibService -ExePath (Join-Path $DistDir $ServiceName) -Pipe $PipeName }

function Compare-Artifact {
  param([string]$Expected,[string]$Target)
  $obj = [ordered]@{ expected=$Expected; target=$Target; exists=$false; sameName=$false; sameSize=$false; sameHash=$false; hashExpected=$null; hashTarget=$null; mtimeExpected=$null; mtimeTarget=$null }
  if(-not (Test-Path $Expected)) { $obj.error = 'expected not found'; return [pscustomobject]$obj }
  if(Test-Path $Target){
    $obj.exists = $true
    $obj.sameName = ([System.IO.Path]::GetFileName($Expected) -eq [System.IO.Path]::GetFileName($Target))
    $e = Get-Item $Expected; $t = Get-Item $Target
    $obj.mtimeExpected = $e.LastWriteTimeUtc.ToString('s'); $obj.mtimeTarget = $t.LastWriteTimeUtc.ToString('s')
    $obj.sameSize = ($e.Length -eq $t.Length)
    $he = (Get-FileHash -Algorithm SHA256 -Path $Expected).Hash
    $obj.hashExpected = $he
    $ht = (Get-FileHash -Algorithm SHA256 -Path $Target).Hash
    $obj.hashTarget = $ht
    $obj.sameHash = ($he -eq $ht)
  }
  else { $obj.error = 'target not found' }
  return [pscustomobject]$obj
}

function Do-Doctor {
  Write-Host '=== Doctor (audit/verify) ==='
  $report = @()
  # dist-wave (fonte) vs instalados em Libraries/Services
  if($UserProfile -and $TerminalId){
    $mql = Join-Path (Get-TerminalRoot $UserProfile $TerminalId) 'MQL5'
    $lib = Join-Path $mql 'Libraries'
    $svc = Join-Path $mql 'Services'
    $report += (Compare-Artifact -Expected (Join-Path $DistDir $DllName) -Target (Join-Path $lib $DllName))
    $report += (Compare-Artifact -Expected (Join-Path $DistDir $ServiceName) -Target (Join-Path $svc $ServiceName))
  }
  # agents dlls
  $agents = Get-ChildItem -Path 'C:\\mql5\\Tester' -Directory -Filter 'Agent-*' -ErrorAction SilentlyContinue
  foreach($a in $agents){ $alib = Join-Path $a.FullName 'MQL5\Libraries'; $report += (Compare-Artifact -Expected (Join-Path $DistDir $DllName) -Target (Join-Path $alib $DllName)) }

  # cmake sanity
  $cmake = Get-Content -Path (Join-Path $AlglibRoot 'CMakeLists.txt') -Raw
  if($cmake -notmatch 'dist-wave'){ Write-Log 'CMakeLists sem dist-wave' 'WARN' }

  # service check (reinicia se -TestService)
  $svc = Get-Process -Name alglib_service -ErrorAction SilentlyContinue | Select-Object -First 1
  $svcOk = $false; $cmd = $null; $pidDisp = '-'
  if($TestService){ Start-AlglibService -ExePath (Join-Path $DistDir $ServiceName) -Pipe $PipeName; Start-Sleep -Milliseconds 500 }
  $svc = Get-Process -Name alglib_service -ErrorAction SilentlyContinue | Select-Object -First 1
  if($svc){ $pidDisp = $svc.Id; $proc = Get-CimInstance Win32_Process -Filter "ProcessId=$($svc.Id)"; $cmd = $proc.CommandLine; if($cmd -and ($cmd -match [Regex]::Escape($PipeName))){ $svcOk = $true } }
  Write-Log ("service pid={0} ok={1}" -f $pidDisp, $svcOk)

  if($AutoFix){
    foreach($r in $report){ if(-not $r.sameHash){ try{ Copy-IfExists $r.expected $r.target } catch { Write-Log "autofix error $($_)" 'ERROR' } } }
    if(-not $svcOk){ Stop-AlglibService; Start-AlglibService -ExePath (Join-Path $DistDir $ServiceName) -Pipe $PipeName }
  }

  if($TestService){ try{ $client = New-Object System.IO.Pipes.NamedPipeClientStream('.', 'alglib-wave_pipe', [System.IO.Pipes.PipeDirection]::InOut, [System.IO.Pipes.PipeOptions]::Asynchronous); $client.Connect(1000); $ok=$client.IsConnected; $client.Dispose(); Write-Log ("pipe-test connected={0}" -f $ok) } catch { Write-Log ("pipe-test error: $_") 'ERROR' } }

  if($ReportPath){ $report | ConvertTo-Json -Depth 5 | Out-File -FilePath $ReportPath -Encoding UTF8 }
  else { $report | Format-Table -AutoSize | Out-String | Write-Host }
}

function Do-Undeploy {
  Write-Host '=== Undeploy (restore from last backup) ==='
  Stop-AlglibService
  $dir = $null
  if($BackupId){ $dir = Join-Path $BackupRoot $BackupId }
  else {
    if(Test-Path $BackupRoot){ $dir = Get-ChildItem -Path $BackupRoot -Directory | Sort-Object Name -Descending | Select-Object -First 1 | ForEach-Object { $_.FullName } }
  }
  if(-not $dir){ throw 'Nenhum backup encontrado. Informe -BackupId ou crie um backup em Consolidate/Deploy.' }
  $manifestPath = Join-Path $dir 'manifest.json'
  if(-not (Test-Path $manifestPath)) { throw "Manifesto não encontrado: $manifestPath" }
  $data = Get-Content -Path $manifestPath -Raw | ConvertFrom-Json
  foreach($f in $data.files){
    Ensure-Folder (Split-Path $f.original -Parent)
    Copy-Item -Path $f.backup -Destination $f.original -Force
    Write-Log ("restore {0} -> {1}" -f $f.backup,$f.original)
  }
  Write-Log ("undeploy OK a partir de {0}" -f $dir)
}

switch($Action){
  'Clean'      { Clean-OldServices }
  'Build'       { Do-Build }
  'Consolidate' { Do-Consolidate }
  'Deploy'      { Do-Deploy }
  'DeployMql'   { Do-DeployMql }
  'Generate'    { Do-Generate }
  'Connect'     { Do-Connect }
  'Compile'     { Do-Compile }
  'Start'       { Do-Start }
  'Doctor'      { Do-Doctor }
  'Undeploy'    { Do-Undeploy }
  'All'         { Clean-OldServices; Do-Build; Do-Deploy; Do-Compile; Do-Start }
}

Write-Host '=== Done ==='
