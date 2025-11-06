param(
  [string]$Pipe = '\\.\pipe\alglib-wave_pipe',
  [int]$TimeoutMs = 3000
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

function New-LEBinaryWriter([System.IO.Stream]$s){
  return New-Object System.IO.BinaryWriter($s,[System.Text.Encoding]::ASCII)
}

function Write-UInt32LE([System.IO.BinaryWriter]$bw,[uint32]$v){ $bw.Write([byte[]][BitConverter]::GetBytes([uint32]$v)) }
function Write-Int32LE ([System.IO.BinaryWriter]$bw,[int32] $v){ $bw.Write([byte[]][BitConverter]::GetBytes([int32] $v)) }
function Write-UInt64LE([System.IO.BinaryWriter]$bw,[uint64]$v){ $bw.Write([byte[]][BitConverter]::GetBytes([uint64]$v)) }

function Read-Exact([System.IO.Stream]$s,[int]$len){
  $buf = New-Object byte[] $len
  $off = 0
  while($off -lt $len){
    $n = $s.Read($buf, $off, $len-$off)
    if($n -le 0){ throw "EOF while reading ($off/$len)" }
    $off += $n
  }
  return $buf
}

function Read-UInt32LE([byte[]]$b,[ref]$o){ $val=[BitConverter]::ToUInt32($b,$o.Value); $o.Value += 4; return $val }
function Read-Int32LE ([byte[]]$b,[ref]$o){ $val=[BitConverter]::ToInt32($b,$o.Value);  $o.Value += 4; return $val }
function Read-UInt64LE([byte[]]$b,[ref]$o){ $val=[BitConverter]::ToUInt64($b,$o.Value); $o.Value += 8; return $val }

Write-Host ("[ping] Connecting to {0}" -f $Pipe)
$client = New-Object System.IO.Pipes.NamedPipeClientStream('.', ($Pipe -replace '^\\\\\.\\pipe\\',''), [System.IO.Pipes.PipeDirection]::InOut, [System.IO.Pipes.PipeOptions]::None)
$client.Connect($TimeoutMs)
if(-not $client.IsConnected){ throw 'connect failed' }
Write-Host '[ping] Connected.'

# Build PROCESS_SUBMIT (PING) header
$bw = New-LEBinaryWriter $client
$MAGIC   = 0x4C574650
$VER     = 1
$CMD_SUB = 20
$tag     = [uint64](Get-Random) -bor ([uint64]([DateTimeOffset]::UtcNow.ToUnixTimeMilliseconds()))
Write-UInt32LE $bw $MAGIC
Write-UInt32LE $bw $VER
Write-UInt32LE $bw $CMD_SUB
Write-UInt64LE $bw $tag
Write-UInt32LE $bw 0   # operation = PING
Write-UInt32LE $bw 0   # flags
Write-Int32LE  $bw 0   # primary_len
Write-Int32LE  $bw 0   # secondary_len
Write-UInt32LE $bw 0   # param_len
$bw.Flush()
Write-Host '[ping] Submit sent.'

# Read ProcessResponse, handle PENDING with PROCESS_FETCH until OK
while($true){
  $hdr = Read-Exact $client 48
  $o = [ref]0
  $magic = Read-UInt32LE $hdr $o
  $ver   = Read-UInt32LE $hdr $o
  $cmd   = Read-UInt32LE $hdr $o
  $rtag  = Read-UInt64LE $hdr $o
  $stat  = Read-Int32LE  $hdr $o
  $pcount= Read-UInt32LE $hdr $o
  $scount= Read-UInt32LE $hdr $o
  $ccount= Read-UInt32LE $hdr $o
  $icount= Read-UInt32LE $hdr $o
  $tcount= Read-UInt32LE $hdr $o
  $psize = Read-UInt32LE $hdr $o
  Write-Host ("[ping] hdr magic=0x{0:X} ver={1} cmd={2} tag={3} status={4} payload={5}" -f $magic,$ver,$cmd,$rtag,$stat,$psize)
  if($magic -ne $MAGIC -or $ver -ne $VER -or $rtag -ne $tag){ throw 'invalid response header' }
  if($stat -eq 1){
    # PENDING: send PROCESS_FETCH
    $CMD_FETCH = 21
    Write-Host '[ping] status=PENDING -> sending FETCH'
    Write-UInt32LE $bw $MAGIC
    Write-UInt32LE $bw $VER
    Write-UInt32LE $bw $CMD_FETCH
    Write-UInt64LE $bw $tag
    $bw.Flush()
    continue
  }
  if($stat -ne 0){ throw ("service returned status {0}" -f $stat) }
  $payload = if($psize -gt 0){ Read-Exact $client $psize } else { @() }
  break
}
$client.Dispose()
if($psize -gt 0){
  try{ $json = [System.Text.Encoding]::UTF8.GetString($payload) } catch { $json = [System.Text.Encoding]::ASCII.GetString($payload) }
  Write-Host "[ping] info: $json"
} else {
  Write-Host '[ping] info: <empty>'
}
