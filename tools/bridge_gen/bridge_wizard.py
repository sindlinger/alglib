#!/usr/bin/env python3
r"""
Bridge Wizard (terminal) — seleciona operações de domínio da frequência
Gera:
 - bridge_mql5/bridge/exports_generated.cpp (exports tipados no DLL)
 - bridge_mql5/mqh/bridge_api.mqh (wrapper único p/ MQL5)
Opcional: compila a DLL e copia para pastas do MetaTrader.
"""
import os, sys, json, subprocess, re

def load_spec(path):
    try:
        import yaml
        with open(path, 'r', encoding='utf-8') as f:
            return yaml.safe_load(f)
    except ImportError:
        # tenta JSON
        with open(path, 'r', encoding='utf-8') as f:
            return json.load(f)

def save_spec(obj, path):
    try:
        import yaml
        with open(path, 'w', encoding='utf-8') as f:
            yaml.safe_dump(obj, f, sort_keys=False, allow_unicode=True)
    except ImportError:
        with open(path, 'w', encoding='utf-8') as f:
            json.dump(obj, f, ensure_ascii=False, indent=2)

def clear():
    os.system('cls' if os.name=='nt' else 'clear')

def prompt(msg):
    try:
        return input(msg)
    except EOFError:
        return ''

def choose_operations(catalog):
    ops = catalog.get('functions', [])
    selected = []
    while True:
        clear()
        print('Selecione operações (digite números separados por vírgula para alternar seleção).')
        print('c = confirmar | o = ordenar | a = selecionar todos | n = limpar | s = scan C++ | q = sair')
        for i, op in enumerate(ops, 1):
            mark = '[x]' if op in selected else '[ ]'
            print(f"{i:2d}. {mark} {op['name']}({', '.join(p['name']+':'+p['type'] for p in op.get('params',[]))})")
        cmd = prompt('> ').strip()
        if cmd.lower() == 'q':
            sys.exit(0)
        if cmd.lower() == 'c':
            return selected
        if cmd.lower() == 'a':
            selected = ops[:]
            continue
        if cmd.lower() == 'n':
            selected = []
            continue
        if cmd.lower() == 'o':
            if not selected:
                prompt('Nada selecionado. ENTER para continuar.'); continue
            clear()
            print('Ordenação: informe a nova ordem pelos índices atuais (ex.: 2,1,3).')
            for i, op in enumerate(selected, 1):
                print(f"{i:2d}. {op['name']}")
            new = prompt('> ').strip()
            try:
                idx = [int(x) for x in new.split(',') if x.strip()]
                selected = [selected[i-1] for i in idx]
            except Exception:
                prompt('Entrada inválida. ENTER para continuar.')
            continue
        if cmd.lower() == 's':
            base = prompt('Diretório base para scan (headers .h/.hpp/.cpp): ').strip()
            if os.path.isdir(base):
                scanned = scan_cpp_headers(base)
                # merge sem duplicar
                names_existing = {x['name'] for x in ops}
                for op in scanned:
                    if op['name'] not in names_existing:
                        ops.append(op)
                prompt(f"Encontradas {len(scanned)} funções. ENTER para continuar.")
            else:
                prompt('Diretório inválido. ENTER para continuar.')
            continue
        # toggle selection
        try:
            idx = [int(x) for x in cmd.split(',') if x.strip()]
            for i in idx:
                if i<1 or i>len(ops):
                    continue
                op = ops[i-1]
                if op in selected:
                    selected.remove(op)
                else:
                    selected.append(op)
        except Exception:
            pass

def map_ctype_to_kind(t: str):
    t = t.strip()
    t = re.sub(r'const\s+', '', t)
    t = t.replace('*', '').replace('&','').strip()
    t = t.lower()
    if t in ('int', 'long', 'int32_t', 'uint32_t'):
        return 'int'
    if t in ('double', 'float', 'float64_t'):
        return 'double'
    if 'wchar_t' in t or 'lpcwstr' in t or 'lpwstr' in t:
        return 'string'  # wide strings
    if 'char' in t or 'lpcstr' in t or 'lpstr' in t:
        return 'string'
    return 'string'

def scan_cpp_headers(base_dir: str):
    pattern = re.compile(r'extern\s+"C"\s+__declspec\(dllexport\)\s+[\w:]+\s+__stdcall\s+(\w+)\s*\(([^)]*)\)', re.IGNORECASE)
    ops = []
    for root, _, files in os.walk(base_dir):
        for fn in files:
            if not fn.lower().endswith(('.h', '.hpp', '.cpp')):
                continue
            fp = os.path.join(root, fn)
            try:
                text = open(fp, 'r', encoding='utf-8', errors='ignore').read()
            except Exception:
                continue
            for m in pattern.finditer(text):
                name = m.group(1)
                params = m.group(2).strip()
                plist = []
                if params and params.lower() != 'void':
                    for p in params.split(','):
                        ptoks = p.strip().split()
                        if not ptoks:
                            continue
                        # last token as name
                        pname = ptoks[-1].replace('*','').replace('&','')
                        ptype = ' '.join(ptoks[:-1]) or 'const char*'
                        kind = map_ctype_to_kind(ptype)
                        plist.append({'name': pname, 'type': kind})
                ops.append({'name': name, 'params': plist})
    # dedup by name, keep first
    seen = set()
    unique = []
    for op in ops:
        if op['name'] in seen:
            continue
        seen.add(op['name'])
        unique.append(op)
    return unique

def main():
    root = os.path.dirname(os.path.abspath(__file__))
    catalog_path = os.path.join(root, 'ops_catalog.yaml')
    spec = load_spec(catalog_path)
    clear()
    print('Backend (pipe/dll)? [pipe]')
    backend = (prompt('> ') or 'pipe').strip().lower()
    if backend == 'pipe':
        print(f"Nome do pipe [{spec.get('pipe_name','\\\\.\\pipe\\alglib-wave_pipe')}] ")
        pn = (prompt('> ') or spec.get('pipe_name','\\\\.\\pipe\\alglib-wave_pipe')).strip()
        spec['pipe_name'] = pn
    else:
        print('Caminho da DLL backend (ex.: C:\\path\\service.dll)')
        dll = prompt('> ').strip()
        spec['backend_dll'] = dll
        spec['pipe_name'] = ''

    selected = choose_operations(spec)
    spec['functions'] = selected

    # Geração via opgen.py (interativo ou --all)
    project_root = os.path.abspath(os.path.join(root, '..', '..'))  # .../alglib
    print("\nGerando arquivos com opgen.py...")
    # opgen próprio já implementa seleção interativa; se o usuário quiser tudo, pressione apenas ENTER
    subprocess.check_call([sys.executable, os.path.join(root, 'opgen.py'), '--root', project_root])
    print('\nArquivos gerados:')
    print(' - src/alglib_gpu/exports_ops_generated.cpp')
    print(' - MQL5/Include/pipe/AlglibOps.mqh')

    # compilar DLL?
    ans = prompt('\nCompilar DLL agora? [S/n] ').strip().lower()
    if ans in ('', 's', 'y', 'yes'):
        br = project_root
        build = os.path.join(br, 'build-win')
        os.makedirs(build, exist_ok=True)
        try:
            subprocess.check_call(['cmake', '-S', br, '-B', build, '-A', 'x64'])
        except Exception:
            # fallback: ambiente Unix/WSL
            subprocess.check_call(['cmake', '-S', br, '-B', build])
        subprocess.check_call(['cmake', '--build', build, '--config', 'Release'])
        print('DLL compilada em dist-wave/alglib_bridge.dll')

    print('\nConcluído. Copie (se não automatizado):')
    print(' - dist-wave/alglib_bridge.dll → MQL5/ Libraries')
    print(' - MQL5/Include/pipe/AlglibOps.mqh → MQL5/ Include/pipe')

if __name__ == '__main__':
    main()
