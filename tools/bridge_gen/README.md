Bridge Generator (MetaTrader <-> ALGLIB service)

- bridge_wizard.py: TUI para escolher operações/ordem e backend (pipe ou DLL);
- opgen.py: gera exports_ops_generated.cpp (DLL) e AlglibOps.mqh (MQL5 alto nível);
- ops_catalog.yaml: catálogo inicial de operações.

Uso rápido (PowerShell / VS Dev Prompt):
  python tools/bridge_gen/bridge_wizard.py
  # depois:
  cmake -S . -B build-win -A x64
  cmake --build build-win --config Release
