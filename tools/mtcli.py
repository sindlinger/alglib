#!/usr/bin/env python3
r"""
mtcli.py — MetaTrader helper CLI for ALGLIB bridge

Commands:
  detect                 Detect Terminal folders (MQL5\Libraries) for current user
  list                   List dist-wave artifacts with times and hashes
  build                  Configure+build via CMake; optional preset regen
  deploy                 Kill MT/Editor/Service, confirm freshness, copy DLLs to Terminals, backup old files (and optional Recycle Bin), print table
  undeploy               Restore previous version from backup (default) or attempt Recycle Bin restore
  kill                   Kill running processes (terminal/metaeditor/alglib_service)

Examples (PowerShell, one per line):
  python tools/mtcli.py detect
  python tools/mtcli.py list
  python tools/mtcli.py build --preset GPU_LegacyWave1.0.4 --use-preset-selection
  python tools/mtcli.py deploy --all --files alglib_bridge.dll --yes

Notes:
  - On Windows, deploy requires admin. Use an elevated shell.
  - Default root is the repository directory that contains this script (.. up to alglib).
"""
import argparse, os, sys, subprocess, hashlib, shutil, datetime, platform


def is_windows():
    return platform.system().lower().startswith('win')

def is_wsl():
    try:
        # typical WSL indicator
        return 'microsoft' in platform.uname().release.lower() or 'WSL_DISTRO_NAME' in os.environ
    except Exception:
        return False


def repo_root_from_here():
    here = os.path.abspath(os.path.dirname(__file__))
    # tools -> alglib
    return os.path.abspath(os.path.join(here, '..'))


def dist_wave_dir(root):
    # Repo layout: ..\dist-wave
    parent = os.path.abspath(os.path.join(root, '..'))
    # In this repo, dist-wave is under C:\mql5\MQL5\dist-wave
    dw = os.path.join(parent, 'dist-wave')
    # If not exist, also check local under root
    return dw if os.path.isdir(dw) else os.path.join(root, 'dist-wave')


def terminal_dirs():
    dirs = []
    if is_windows():
        appdata = os.environ.get('APPDATA')
        if appdata:
            base = os.path.join(appdata, 'MetaQuotes', 'Terminal')
            if os.path.isdir(base):
                for name in os.listdir(base):
                    p = os.path.join(base, name)
                    if os.path.isdir(p):
                        libs = os.path.join(p, 'MQL5', 'Libraries')
                        if os.path.isdir(libs):
                            dirs.append({'id': name, 'root': p, 'libs': libs})
    elif is_wsl():
        # Fallback scan on Windows drive from WSL
        base_users = '/mnt/c/Users'
        if os.path.isdir(base_users):
            for user in os.listdir(base_users):
                tbase = os.path.join(base_users, user, 'AppData', 'Roaming', 'MetaQuotes', 'Terminal')
                if not os.path.isdir(tbase):
                    continue
                for name in os.listdir(tbase):
                    p = os.path.join(tbase, name)
                    if os.path.isdir(p):
                        libs = os.path.join(p, 'MQL5', 'Libraries')
                        if os.path.isdir(libs):
                            dirs.append({'id': name, 'root': p, 'libs': libs})
    return dirs


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def fmt_dt(ts):
    return datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')


def check_admin_windows():
    if not is_windows():
        return True
    try:
        import ctypes
        return bool(ctypes.windll.shell32.IsUserAnAdmin())
    except Exception:
        return False


def run(cmd, cwd=None, check=True):
    print('> ' + ' '.join(cmd))
    return subprocess.run(cmd, cwd=cwd, check=check)


def cmd_detect(_args):
    tdirs = terminal_dirs()
    if not tdirs:
        print('No Terminal directories found under %APPDATA%/MetaQuotes/Terminal')
        return 0
    print('Detected Terminal directories:')
    for d in tdirs:
        print(f"- id={d['id']} libs={d['libs']}")
    return 0


def list_file(path):
    if not os.path.isfile(path):
        print(f'MISSING: {path}')
        return
    st = os.stat(path)
    age_min = int((datetime.datetime.now().timestamp() - st.st_mtime) / 60)
    print(f"{path:80s} | {st.st_size:10d} B | {fmt_dt(st.st_mtime)} | age ~{age_min:3d}m | {sha256_of(path)[:16]}")


def cmd_list(args):
    root = args.root or repo_root_from_here()
    dw = dist_wave_dir(root)
    files = [
        os.path.join(dw, 'alglib_bridge.dll'),
        os.path.join(dw, 'alglib.dll'),
        os.path.join(dw, 'alglib_service.exe'),
    ]
    print('Artifacts:')
    for f in files:
        list_file(f)
    return 0


def cmd_build(args):
    root = args.root or repo_root_from_here()
    build = args.build_dir or os.path.join(root, 'build-win')
    arch = args.arch or 'x64'
    config = args.config or 'Release'

    # Optional preset regeneration
    if args.preset:
        opgen = os.path.join(root, 'tools', 'bridge_gen', 'opgen.py')
        if not os.path.isfile(opgen):
            print('ERROR: opgen.py not found at', opgen)
            return 2
        regen_cmd = [sys.executable, opgen, '--root', root, '--preset', args.preset]
        if args.use_preset_selection:
            regen_cmd.append('--use-preset-selection')
        if args.no_union_presets:
            regen_cmd.append('--no-union-presets')
        run(regen_cmd)

    # Configure
    cfg_cmd = ['cmake', '-S', root, '-B', build]
    if is_windows():
        cfg_cmd += ['-A', arch]
    run(cfg_cmd)

    # Build
    build_cmd = ['cmake', '--build', build, '--config', config]
    if args.targets:
        build_cmd += ['--target'] + args.targets
    if args.parallel:
        build_cmd += ['-j', str(args.parallel)]
    # Pass-through extras
    if args.extra:
        build_cmd += args.extra
    run(build_cmd)
    return 0


def kill_processes(names):
    if not is_windows():
        return
    for n in names:
        # taskkill ignores if not found
        try:
            run(['taskkill', '/F', '/IM', n], check=False)
        except Exception:
            pass


def cmd_kill(_args):
    if not is_windows():
        print('kill only implemented for Windows.')
        return 1
    kill_processes(['terminal64.exe', 'terminal.exe', 'metaeditor64.exe', 'metaeditor.exe', 'alglib_service.exe'])
    return 0


def print_table(rows, headers):
    widths = [len(h) for h in headers]
    for row in rows:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(str(cell)))
    fmt = '  '.join('{:' + str(w) + '}' for w in widths)
    print(fmt.format(*headers))
    print(fmt.format(*['-' * w for w in widths]))
    for row in rows:
        print(fmt.format(*row))


def cmd_deploy(args):
    # On Windows we check admin; on WSL we proceed with best-effort
    if is_windows() and not check_admin_windows():
        print('ERROR: deploy requires Administrator (elevated) shell.')
        return 2

    root = args.root or repo_root_from_here()
    dw = dist_wave_dir(root)
    files = args.files or ['alglib_bridge.dll']
    srcs = []
    for fn in files:
        p = fn if os.path.isabs(fn) else os.path.join(dw, fn)
        if not os.path.isfile(p):
            print('ERROR: missing file', p)
            return 3
        srcs.append(p)

    # freshness check
    if not args.yes:
        print('About to deploy these files:')
        for s in srcs:
            st = os.stat(s)
            print(f"- {os.path.basename(s)}  {fmt_dt(st.st_mtime)}  size={st.st_size}")
        ans = input('Proceed to deploy? [y/N] ').strip().lower()
        if ans not in ('y', 'yes'):
            print('Aborted.')
            return 0

    # Kill processes
    if not args.no_kill:
        print('Stopping Terminal/MetaEditor/Service...')
        if is_windows():
            kill_processes(['terminal64.exe', 'terminal.exe', 'metaeditor64.exe', 'metaeditor.exe', 'alglib_service.exe'])
        elif is_wsl():
            # try via powershell.exe
            try:
                run(['powershell.exe','-NoProfile','-Command','taskkill /F /IM terminal64.exe; taskkill /F /IM terminal.exe; taskkill /F /IM metaeditor64.exe; taskkill /F /IM metaeditor.exe; taskkill /F /IM alglib_service.exe'], check=False)
            except Exception:
                pass

    # Destinations
    tdirs = terminal_dirs()
    if not tdirs:
        print('No Terminal dirs found.')
        return 4
    targets = []
    if args.all:
        targets = tdirs
    elif args.terminal_id:
        for d in tdirs:
            if d['id'] == args.terminal_id:
                targets = [d]; break
        if not targets:
            print('Terminal id not found:', args.terminal_id)
            return 5
    else:
        # default: all
        targets = tdirs

    # Prepare backup folder
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    backup_root = os.path.join(root, 'tools', 'mtcli_backups', stamp)

    # Optional trash support
    def send_to_trash(path):
        try:
            import send2trash  # optional dependency
            send2trash.send2trash(path)
            return True
        except Exception:
            return False

    # Copy
    rows = []
    for d in targets:
        for s in srcs:
            dst = os.path.join(d['libs'], os.path.basename(s))
            os.makedirs(d['libs'], exist_ok=True)
            # backup existing
            if os.path.isfile(dst):
                bdst_dir = os.path.join(backup_root, d['id'])
                os.makedirs(bdst_dir, exist_ok=True)
                bdst = os.path.join(bdst_dir, os.path.basename(dst))
                shutil.copy2(dst, bdst)
                if getattr(args, 'trash', False):
                    send_to_trash(dst)
            shutil.copy2(s, dst)
            st = os.stat(dst)
            rows.append([
                d['id'], os.path.basename(dst), fmt_dt(st.st_ctime), fmt_dt(st.st_mtime), str(st.st_size), sha256_of(dst)[:16], d['libs']
            ])

    print('\nDeployed files:')
    print_table(rows, headers=['TerminalID','File','Created','Modified','Bytes','SHA256(16)','Destination'])
    return 0


def latest_backup_for(root, terminal_id, filename):
    base = os.path.join(root, 'tools', 'mtcli_backups')
    if not os.path.isdir(base):
        return None
    stamps = sorted([n for n in os.listdir(base) if os.path.isdir(os.path.join(base, n))], reverse=True)
    for st in stamps:
        p = os.path.join(base, st, terminal_id, filename)
        if os.path.isfile(p):
            return p
    return None


def restore_from_recycle_bin(dst_path):
    # Best effort: use PowerShell COM to find and restore a recycled item matching original name and destination folder
    if not is_windows():
        return False
    name = os.path.basename(dst_path)
    folder = os.path.dirname(dst_path)
    script = f"$s=New-Object -ComObject Shell.Application; $bin=$s.NameSpace(0xA); $items=$bin.Items(); $restored=$false; foreach($i in $items){{ $on=$bin.GetDetailsOf($i,0); $ol=$bin.GetDetailsOf($i,1); if($on -eq '{name}' -and $ol -like '*MQL5\\Libraries*'){{ try{{$i.InvokeVerb('RESTORE')}}catch{{}} $restored=$true; break }} }}; if($restored){{exit 0}} else {{exit 1}}"
    try:
        r = subprocess.run(['powershell','-NoProfile','-Command',script], capture_output=True)
        return r.returncode == 0
    except Exception:
        return False


def cmd_undeploy(args):
    # On Windows we check admin; on WSL we proceed best-effort
    if is_windows() and not check_admin_windows():
        print('ERROR: undeploy may require Administrator (elevated) shell for process kill/copy.')
    root = args.root or repo_root_from_here()
    files = args.files or ['alglib_bridge.dll']
    # Kill processes if requested
    if not args.no_kill:
        kill_processes(['terminal64.exe', 'terminal.exe', 'metaeditor64.exe', 'metaeditor.exe', 'alglib_service.exe'])

    tdirs = terminal_dirs()
    if not tdirs:
        print('No Terminal dirs found.')
        return 2
    targets = []
    if args.all:
        targets = tdirs
    elif args.terminal_id:
        for d in tdirs:
            if d['id'] == args.terminal_id:
                targets = [d]; break
        if not targets:
            print('Terminal id not found:', args.terminal_id)
            return 3
    else:
        targets = tdirs

    rows = []
    for d in targets:
        for fn in files:
            dst = os.path.join(d['libs'], fn)
            restored = False
            reason = ''
            if args.from_recycle:
                restored = restore_from_recycle_bin(dst)
                reason = 'recycle'
            if not restored:
                # try backup
                bk = latest_backup_for(root, d['id'], fn)
                if bk:
                    os.makedirs(d['libs'], exist_ok=True)
                    shutil.copy2(bk, dst)
                    restored = True
                    reason = 'backup'
            st = os.stat(dst) if os.path.isfile(dst) else None
            rows.append([d['id'], fn, 'OK' if restored else 'MISS', reason, fmt_dt(st.st_mtime) if st else '-', str(st.st_size) if st else '-', sha256_of(dst)[:16] if st else '-', d['libs']])

    print('\nUndeploy result:')
    print_table(rows, headers=['TerminalID','File','Restored','Source','Modified','Bytes','SHA256(16)','Destination'])
    return 0


def main():
    ap = argparse.ArgumentParser(prog='mtcli', description='MetaTrader CLI for ALGLIB bridge')
    sub = ap.add_subparsers(dest='cmd', required=True)

    p = sub.add_parser('detect', help='Detect Terminal directories')
    p.set_defaults(func=cmd_detect)

    p = sub.add_parser('list', help='List dist-wave artifacts')
    p.add_argument('--root', help='Repo root (default: auto)')
    p.set_defaults(func=cmd_list)

    p = sub.add_parser('build', help='Configure+build via CMake; optional preset regen')
    p.add_argument('--root', help='Repo root (default: auto)')
    p.add_argument('--build-dir', help='Build directory (default: <root>/build-win)')
    p.add_argument('--arch', default='x64', help='Windows arch for CMake -A (default: x64)')
    p.add_argument('--config', default='Release', help='Build config (default: Release)')
    p.add_argument('--targets', nargs='+', help='Build targets (default: ALL)')
    p.add_argument('--parallel', type=int, default=os.cpu_count() or 8, help='-j (default: cpu count)')
    p.add_argument('--preset', help='Run opgen.py to (re)generate preset include and exports')
    p.add_argument('--use-preset-selection', action='store_true', help='Use exactly preset ops without prompt')
    p.add_argument('--no-union-presets', action='store_true', help='Do not union presets for exports')
    p.add_argument('extra', nargs=argparse.REMAINDER, help='Extra args after -- passed to cmake --build')
    p.set_defaults(func=cmd_build)

    p = sub.add_parser('kill', help='Kill Terminal/MetaEditor/Service processes')
    p.set_defaults(func=cmd_kill)

    p = sub.add_parser('deploy', help='Deploy DLL(s) to Terminal MQL5\\Libraries (admin only)')
    p.add_argument('--root', help='Repo root (default: auto)')
    p.add_argument('--files', nargs='+', help='Files to copy (default: dist-wave/alglib_bridge.dll)')
    p.add_argument('--terminal-id', help='Specific Terminal ID (default: all)')
    p.add_argument('--all', action='store_true', help='Deploy to all Terminals')
    p.add_argument('--yes', action='store_true', help='Do not prompt for confirmation')
    p.add_argument('--no-kill', action='store_true', help='Do not kill processes before copying')
    p.add_argument('--trash', action='store_true', help='Send replaced files to Recycle Bin (requires send2trash, optional)')
    p.set_defaults(func=cmd_deploy)

    p = sub.add_parser('undeploy', help='Restore previous files from backup or Recycle Bin')
    p.add_argument('--root', help='Repo root (default: auto)')
    p.add_argument('--files', nargs='+', help='Files to restore (default: alglib_bridge.dll)')
    p.add_argument('--terminal-id', help='Specific Terminal ID (default: all)')
    p.add_argument('--all', action='store_true', help='Apply to all Terminals')
    p.add_argument('--from-recycle', action='store_true', help='Attempt restore from Recycle Bin (best effort)')
    p.add_argument('--no-kill', action='store_true', help='Do not kill processes before restore')
    p.set_defaults(func=cmd_undeploy)

    args = ap.parse_args()
    return args.func(args)


if __name__ == '__main__':
    sys.exit(main())
