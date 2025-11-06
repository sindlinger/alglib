param(
  [string]$Pipe = '\\.\pipe\alglib-wave_pipe',
  [int]$TimeoutMs = 5000
)

# Constants
$MAGIC   = [UInt32]0x4C574650
$VERSION = [UInt32]1
$CMD_SUBMIT = [UInt32]20
$CMD_FETCH  = [UInt32]21

Add-Type -AssemblyName System.Core
Add-Type -AssemblyName System.IO.Compression.FileSystem

function Write-UInt32([System.IO.BinaryWriter]$bw, [UInt32]$v){ $bw.Write($v) }
function Write-Int32 ([System.IO.BinaryWriter]$bw, [Int32] $v){ $bw.Write($v) }
function Write-UInt64([System.IO.BinaryWriter]$bw, [UInt64]$v){ $bw.Write($v) }

function Read-Exact([System.IO.Stream]$s, [byte[]]$buf, [int]$len){
  $off = 0
  while($off -lt $len){
    $r = $s.Read($buf, $off, $len - $off)
    if($r -le 0){ throw "stream closed while reading ($off/$len)" }
    $off += $r
  }
}

try{
  $client = New-Object System.IO.Pipes.NamedPipeClientStream('.', ($Pipe -replace '^.*\\',''), [System.IO.Pipes.PipeDirection]::InOut, [System.IO.Pipes.PipeOptions]::None)
  Write-Host "[ping-min] Connecting to $Pipe"
  $client.Connect($TimeoutMs)
  Write-Host "[ping-min] Connected."

  $tag = [UInt64]([DateTimeOffset]::UtcNow.ToUnixTimeMilliseconds())

  $ms = New-Object System.IO.MemoryStream
  $bw = New-Object System.IO.BinaryWriter($ms)
  Write-UInt32 $bw $MAGIC
  Write-UInt32 $bw $VERSION
  Write-UInt32 $bw $CMD_SUBMIT
  Write-UInt64 $bw $tag
  Write-UInt32 $bw 0     # operation = PING
  Write-UInt32 $bw 0     # flags
  Write-Int32  $bw 0     # primary_len
  Write-Int32  $bw 0     # secondary_len
  Write-UInt32 $bw 0     # param_len
  $bw.Flush()
  $bytes = $ms.ToArray()
  $client.Write($bytes, 0, $bytes.Length)
  $client.Flush()
  Write-Host "[ping-min] Submit sent (tag=$tag)."

  $deadline = [DateTimeOffset]::UtcNow.ToUnixTimeMilliseconds() + $TimeoutMs
  $payload = $null
  $psize   = 0
  while($true){
    # Read 48-byte ProcessResponse header
    $hdr = New-Object byte[] 48
    Read-Exact -s $client -buf $hdr -len 48
    $o = 0
    $r_magic = [BitConverter]::ToUInt32($hdr, $o); $o+=4
    $r_ver   = [BitConverter]::ToUInt32($hdr, $o); $o+=4
    $r_cmd   = [BitConverter]::ToUInt32($hdr, $o); $o+=4
    $r_tag   = [BitConverter]::ToUInt64($hdr, $o); $o+=8
    $status  = [BitConverter]::ToInt32($hdr, $o); $o+=4
    $pcount  = [BitConverter]::ToUInt32($hdr, $o); $o+=4
    $scount  = [BitConverter]::ToUInt32($hdr, $o); $o+=4
    $ccount  = [BitConverter]::ToUInt32($hdr, $o); $o+=4
    $icount  = [BitConverter]::ToUInt32($hdr, $o); $o+=4
    $tcount  = [BitConverter]::ToUInt32($hdr, $o); $o+=4
    $psize   = [BitConverter]::ToUInt32($hdr, $o); $o+=4

    Write-Host ([string]::Format('[ping-min] hdr magic=0x{0:X8} ver={1} cmd={2} tag={3} status={4} payload={5}', $r_magic, $r_ver, $r_cmd, $r_tag, $status, $psize))

    if($status -eq 0){
      if($psize -gt 0){
        $payload = New-Object byte[] $psize
        Read-Exact -s $client -buf $payload -len $psize
      }
      break
    }
    elseif($status -eq 1){
      if([DateTimeOffset]::UtcNow.ToUnixTimeMilliseconds() -gt $deadline){ throw 'timeout awaiting OK (PENDING loop)' }
      # Send FETCH for same tag
      $ms2 = New-Object System.IO.MemoryStream
      $bw2 = New-Object System.IO.BinaryWriter($ms2)
      Write-UInt32 $bw2 $MAGIC
      Write-UInt32 $bw2 $VERSION
      Write-UInt32 $bw2 $CMD_FETCH
      Write-UInt64 $bw2 $r_tag
      $bw2.Flush()
      $bytes2 = $ms2.ToArray()
      $client.Write($bytes2, 0, $bytes2.Length)
      $client.Flush()
      Start-Sleep -Milliseconds 10
      continue
    }
    else{
      throw "service returned status $status"
    }
  }

  $json = if($psize -gt 0){ [System.Text.Encoding]::UTF8.GetString($payload) } else { '{}' }
  Write-Host "[ping-min] JSON: $json"
  $client.Dispose()
  exit 0
}
catch{
  Write-Host "[ping-min][ERROR] $_" -ForegroundColor Red
  exit 1
}
