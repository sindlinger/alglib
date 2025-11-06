ALGLIB Bridge (DLL de Ponte)
================================

Resumo
- Objetivo: A DLL `alglib_bridge.dll` NÃO executa computação. Ela apenas faz a ponte de transporte entre MQL5 e o serviço `alglib_service.exe` via Named Pipe (\\.\\pipe\\alglib-wave_pipe), suportando rotas síncronas, assíncronas e batch.
- Motivação: manter um único caminho (Via A) para GPU (serviço), sem import direto de DLLs de computação no MQL. A ponte facilita integração (via `#import`) quando preferir evitar MQL puro de pipe.

Binários e caminhos
- Serviço: `C:\mql5\MQL5\dist-wave\alglib_service.exe`
- Ponte:   `C:\mql5\MQL5\dist-wave\alglib_bridge.dll`
- Deploy padrão (instância): `...\MQL5\Libraries\alglib_bridge.dll`
- Pipe: `\\.\\pipe\\alglib-wave_pipe` (sem “2”).

Protocolo (resumo)
- Mensagens: PROCESS_SUBMIT (20) / PROCESS_FETCH (21) e BATCH_SUBMIT (30) / BATCH_FETCH (31).
- Operações: FFT_REAL/COMPLEX/SINE/COSINE/TWO_REAL, SPECTRAL_* (phase_unwrap, instant, filters, analyze, detect_transitions), PING=0.
- Respostas: cabeçalho com status (OK=0, PENDING=1, …) + payload com vetores/results conforme operação.

Exports (C)
```
int AlglibBridge_Ping(unsigned char* out_json, int cap, int* out_len, int timeout_ms);

int AlglibBridge_Process(
  unsigned int op,
  const double* primary,  int primary_len,
  const double* secondary,int secondary_len,
  const unsigned char* params, int params_len,
  double* out_primary, int out_primary_cap, int* out_primary_len,
  double* out_secondary,int out_secondary_cap,int* out_secondary_len,
  unsigned char* out_payload,int out_payload_cap,int* out_payload_len,
  int timeout_ms);

int AlglibBridge_ProcessSubmit(/* mesmos parâmetros sem os buffers de saída */ unsigned long long* out_handle, int timeout_ms);
int AlglibBridge_ProcessFetch(unsigned long long handle,
  double* out_primary,int cap_primary,int* out_primary_len,
  double* out_secondary,int cap_secondary,int* out_secondary_len,
  unsigned char* out_payload,int cap_payload,int* out_payload_len,
  int timeout_ms);

int AlglibBridge_BatchSubmit(/* descritores e blobs */ unsigned long long* out_handle, int timeout_ms);
int AlglibBridge_BatchFetch(unsigned long long handle, unsigned char* out_payload, int cap, int* out_len, int* out_result_count, int timeout_ms);
```

Características
- Zero computação: toda operação é encaminhada ao serviço; a DLL apenas empacota/desempacota mensagens.
- Síncrono/Assíncrono/Batch: aderente ao protocolo do serviço.
- Thread-safe: cada chamada abre/usa handle próprio (stateless do ponto de vista de computação).

Uso em MQL5 (exemplo mínimo)
```
#import "alglib_bridge.dll"
  int  AlglibBridge_Ping(uchar &out_json[], int cap, int &out_len, int timeout_ms);
#import

void OnStart(){
  uchar buf[]; ArrayResize(buf, 1024);
  int out_len=0;
  int st = AlglibBridge_Ping(buf, ArraySize(buf), out_len, 5000);
  if(st==0){ string json = CharArrayToString(buf,0,out_len); Print("[Bridge] ", json); }
  else{ Print("[Bridge][ERR] status=", st); }
}
```

Boas práticas
- Sempre inicie o serviço com o script estrito (sem fallback): `alglib\start_service_strict.ps1`.
- Não importe DLLs de computação no indicador; use somente a ponte (quando preferir DLL) ou o include gerado de pipe.
- Tenha um único include gerado por preset (ex.: `GPU_LegacyWave1.0.4.mqh`) quando optar por MQL puro.

Validação rápida
- Serviço: `start_service_strict.ps1` → mostra o PATH exato do binário em execução.
- Ping sem CLI: `alglib\test_ping_min.ps1` → imprime JSON do PING.

Histórico
- 2025-11-06: build verificado em `dist-wave\alglib_bridge.dll` e copiado para a instância (Libraries).
