param(
  [string]$UserProfile = "$env:USERPROFILE",
  [Parameter(Mandatory=$true)][string]$TerminalId,
  [string]$DestinationRoot = 'C:\\mql5\\Terminal',
  [string]$DestinationPath,
  [switch]$ForceKill,
  [switch]$DryRun
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

function Write-Log([string]$Message,[string]$Level='INFO'){
  $ts = (Get-Date).ToString('yyyy-MM-dd HH:mm:ss.fff')
  Write-Host "[$ts][$Level] $Message"
}

function Ensure-Folder([string]$Path){ if(-not (Test-Path $Path)){ New-Item -ItemType Directory -Path $Path -Force | Out-Null } }

$src = Join-Path (Join-Path $UserProfile 'AppData/Roaming/MetaQuotes/Terminal') $TerminalId
if([string]::IsNullOrWhiteSpace($DestinationPath)){
  $dst = Join-Path $DestinationRoot $TerminalId
} else {
  $dst = $DestinationPath
}

Write-Log ("Fonte: {0}" -f $src)
Write-Log ("Destino: {0}" -f $dst)

if(-not (Test-Path $src)){ throw "Pasta de origem não encontrada: $src" }
Ensure-Folder (Split-Path $dst -Parent)

# 1) Garantir que Terminal/MetaEditor/Agents estejam parados
$procs = @('terminal64','metaeditor64')
foreach($p in $procs){
  $list = Get-Process -Name $p -ErrorAction SilentlyContinue
  if($list){
    if($ForceKill){ $list | Stop-Process -Force -ErrorAction SilentlyContinue; Write-Log ("matado: {0}" -f $p) 'WARN' }
    else { throw "Processo em execução: $p (use -ForceKill ou feche manualmente)" }
  }
}

# 2) Copiar conteúdo com preservação de ACLs/horários
if($DryRun){ Write-Log '[DryRun] pular robocopy' }
else{
  Ensure-Folder $dst
  $rcmd = @('robocopy', '"{0}"' -f $src, '"{0}"' -f $dst, '/MIR','/COPYALL','/R:1','/W:1') -join ' '
  Write-Log ("Executando: {0}" -f $rcmd)
  $rob = Start-Process -FilePath robocopy.exe -ArgumentList @($src,$dst,'/MIR','/COPYALL','/R:1','/W:1') -NoNewWindow -PassThru -Wait
}

# 3) Renomear origem e criar junction no lugar
$bak = "{0}_bak_{1}" -f $src, (Get-Date).ToString('yyyyMMdd_HHmmss')
if($DryRun){ Write-Log ("[DryRun] renomear {0} -> {1}" -f $src,$bak) }
else{
  Rename-Item -Path $src -NewName (Split-Path $bak -Leaf)
  Write-Log ("renomeado: {0} -> {1}" -f $src,$bak)
  $mk = "cmd /c mklink /J \"{0}\" \"{1}\"" -f $src,$dst
  Write-Log ("Criando junction: {0}" -f $mk)
  $p = Start-Process -FilePath cmd.exe -ArgumentList @('/c','mklink','/J',"$src","$dst") -PassThru -Wait -NoNewWindow
}

Write-Log 'Concluído. O caminho original agora aponta para C:\mql5.'
