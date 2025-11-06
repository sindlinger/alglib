param(
  [string]$ExePath = 'C:\\mql5\\MQL5\\dist-wave\\Release\\alglib_service.exe',
  [string]$Pipe    = '\\\\.\\pipe\\alglib-wave_pipe',
  [switch]$ShowLog
)

$ErrorActionPreference = 'Stop'

function Stop-AllAlglib {
  $names = @('alglib_service','ALGLIB_Service','Gen4EngineService','Gen4Engine')
  foreach($n in $names){ Get-Process -Name $n -ErrorAction SilentlyContinue | Stop-Process -Force -ErrorAction SilentlyContinue }
}

if(-not (Test-Path $ExePath)){ throw "alglib_service não encontrado: $ExePath" }

Write-Host "[start] Killing any previous alglib_service..."
Stop-AllAlglib
Start-Sleep -Milliseconds 200

Write-Host "[start] Starting: $ExePath --pipe $Pipe"
$p = Start-Process -FilePath $ExePath -ArgumentList "--pipe $Pipe" -PassThru -WindowStyle Normal
Start-Sleep -Milliseconds 300

# Report process path via WMI (robusto)
$wmi = Get-CimInstance Win32_Process -Filter "ProcessId=$($p.Id)" -ErrorAction SilentlyContinue
if($wmi){ Write-Host ("[start] Running PID={0} PATH='{1}'" -f $p.Id, $wmi.ExecutablePath) }
else{ Write-Host ("[start] Running PID={0}" -f $p.Id) }

if($ShowLog){
  $log1 = Join-Path (Split-Path $ExePath -Parent) 'alglib_service.log'
  if(Test-Path $log1){ Get-Content -Path $log1 -Tail 20 -Wait }
}

