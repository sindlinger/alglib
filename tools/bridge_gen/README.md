Bridge Generator (MetaTrader <-> ALGLIB service)

- bridge_wizard.py: TUI para escolher operações/ordem e backend (pipe ou DLL);
- opgen.py: gera exports_ops_generated.cpp (DLL) e AlglibOps.mqh (MQL5 alto nível);
- ops_catalog.yaml: catálogo inicial de operações.

Presets por indicador
- Use presets para gerar um .mqh específico por indicador e manter uma única DLL com a união de todas as operações necessárias:
  - Salvar e usar preset (nome do indicador):
    python tools/bridge_gen/opgen.py --root C:\mql5\MQL5\alglib --preset GPU_LegacyWave1.0.4 --all
    # Saídas:
    #  - src/alglib_gpu/exports_ops_generated.cpp (união de todos os presets encontrados)
    #  - MQL5/Include/pipe/AlglibOps_GPU_LegacyWave1.0.4.mqh
    #  - tools/bridge_gen/presets/GPU_LegacyWave1.0.4.yaml

  - O arquivo exports_ops_generated.cpp é gerado com a união das operações de todos os presets em tools/bridge_gen/presets/*.yaml, garantindo que uma única alglib_bridge.dll atenda a todos os indicadores.
  - Para desativar a união e gerar exports apenas da seleção atual, use --no-union-presets.

  - Para customizar o nome/locais do .mqh:
    python tools/bridge_gen/opgen.py --root ... --preset MeuIndicador \
      --mqh-name AlglibOps_MeuIndicador.mqh --include-dir C:\\mql5\\MQL5\\Include\\pipe


Uso rápido (PowerShell / VS Dev Prompt):
  python tools/bridge_gen/bridge_wizard.py
  # depois:
  cmake -S . -B build-win -A x64
  cmake --build build-win --config Release
