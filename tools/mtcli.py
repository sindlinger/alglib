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
import json


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


def to_winpath(p: str) -> str:
    """Convert /mnt/c/... to C:\... when running under WSL for PowerShell tools."""
    if not is_wsl() or not isinstance(p, str):
        return p
    if p.startswith('/mnt/') and len(p) > 6 and p[6] == '/':
        drive = p[5].upper()
        rest = p[7:].replace('/', '\\')
        return f"{drive}:\\{rest}"
    return p


# ---------------- Project persistence -----------------
def _projects_path(root: str) -> str:
    return os.path.join(root, 'tools', 'mtcli_projects.json')

def load_projects(root: str):
    path = _projects_path(root)
    if not os.path.isfile(path):
        return {"last_project": None, "projects": {}}
    try:
        with open(path, 'r', encoding='utf-8') as f:
            return json.load(f)
    except Exception:
        return {"last_project": None, "projects": {}}

def save_projects(root: str, data):
    path = _projects_path(root)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + '.tmp'
    with open(tmp, 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    os.replace(tmp, path)


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

    # Always show source file metadata (name, bytes, modified, sha) and keep for verification
    src_meta = []
    src_hash_by_name = {}
    for s in srcs:
        st = os.stat(s)
        h = sha256_of(s)
        src_meta.append([os.path.basename(s), str(st.st_size), fmt_dt(st.st_mtime), h[:16], s])
        src_hash_by_name[os.path.basename(s)] = h
    print('\nSource files:')
    print_table(src_meta, headers=['File','Bytes','Modified','SHA256(16)','Path'])

    # freshness check/prompt unless forced
    if not args.yes:
        ans = input('\nProceed to deploy? [y/N] ').strip().lower()
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
    if getattr(args, 'libs', None):
        tdirs = [{'id': f'manual{idx+1}', 'root': os.path.dirname(os.path.dirname(p)), 'libs': p} for idx, p in enumerate(args.libs)]
    else:
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
            d_hash = sha256_of(dst)
            s_hash = src_hash_by_name.get(os.path.basename(s), '')
            verified = 'OK' if d_hash == s_hash else 'FAIL'
            rows.append([
                d['id'], os.path.basename(dst), fmt_dt(st.st_ctime), fmt_dt(st.st_mtime), str(st.st_size), d_hash[:16], verified, d['libs']
            ])

    print('\nDeployed files (verification):')
    print_table(rows, headers=['TerminalID','File','Created','Modified','Bytes','SHA256(16)','Verified','Destination'])
    # Summary
    fails = sum(1 for r in rows if r[6] != 'OK')
    if fails == 0:
        print(f"\nSummary: {len(rows)} file copies verified OK across {len(targets)} Terminal(s).")
    else:
        print(f"\nSummary: {fails} verification failure(s) out of {len(rows)} copies.")
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
    p.add_argument('--libs', action='append', help='Explicit MQL5\\Libraries directory (repeatable). Overrides detection when present')
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

    # ---- project (save/list/use/show) ----
    proj = sub.add_parser('project', help='Manage saved projects (paths, presets, defaults)')
    proj_sub = proj.add_subparsers(dest='pcmd', required=True)
    def _cmd_proj_save(args):
        root = args.root or repo_root_from_here()
        data = load_projects(root)
        pid = args.id
        if not pid:
            print('ERROR: --id is required'); return 2
        now = datetime.datetime.now().isoformat(timespec='seconds')
        proj = data.get('projects', {}).get(pid, {})
        proj.update({
            'project': args.name or proj.get('project') or pid,
            'libs': args.libs or proj.get('libs',''),
            'metaeditor': args.metaeditor or proj.get('metaeditor',''),
            'updated_at': now,
        })
        if 'created_at' not in proj:
            proj['created_at'] = now
        data.setdefault('projects', {})[pid] = proj
        if args.set_default:
            data['last_project'] = pid
        save_projects(root, data)
        print('Saved project:', pid)
        return 0
    ps = proj_sub.add_parser('save', help='Save or update a project')
    ps.add_argument('--root', help='Repo root (default: auto)')
    ps.add_argument('--id', required=True, help='Project id (e.g., legacy-wave)')
    ps.add_argument('--name', help='Preset/name for generator (e.g., GPU_LegacyWave1.0.4)')
    ps.add_argument('--libs', help='Path to MQL5\\Libraries for this project')
    ps.add_argument('--metaeditor', help='Path to metaeditor64.exe for this project')
    ps.add_argument('--set-default', action='store_true', help='Set as last_project')
    ps.set_defaults(func=_cmd_proj_save)

    def _cmd_proj_list(args):
        root = args.root or repo_root_from_here()
        data = load_projects(root)
        rows = []
        last = data.get('last_project')
        for pid, proj in data.get('projects', {}).items():
            mark = '*' if pid == last else ' '
            rows.append([mark+pid, proj.get('project',''), proj.get('libs',''), proj.get('metaeditor',''), proj.get('updated_at','')])
        if rows:
            print_table(rows, headers=['ProjectID','Name','Libraries','MetaEditor','Updated'])
        else:
            print('No saved projects.')
        return 0
    pl = proj_sub.add_parser('list', help='List projects')
    pl.add_argument('--root', help='Repo root (default: auto)')
    pl.set_defaults(func=_cmd_proj_list)

    def _cmd_proj_use(args):
        root = args.root or repo_root_from_here()
        data = load_projects(root)
        if args.id not in data.get('projects', {}):
            print('ERROR: unknown project id:', args.id); return 2
        data['last_project'] = args.id
        save_projects(root, data)
        print('Now using project:', args.id)
        return 0
    pu = proj_sub.add_parser('use', help='Set default (last) project')
    pu.add_argument('--root', help='Repo root (default: auto)')
    pu.add_argument('--id', required=True)
    pu.set_defaults(func=_cmd_proj_use)

    def _cmd_proj_show(args):
        root = args.root or repo_root_from_here()
        data = load_projects(root)
        pid = args.id or data.get('last_project')
        proj = data.get('projects', {}).get(pid)
        if not proj:
            print('No project selected or not found.'); return 1
        rows = [[pid, proj.get('project',''), proj.get('libs',''), proj.get('metaeditor',''), proj.get('updated_at','')]]
        print_table(rows, headers=['ProjectID','Name','Libraries','MetaEditor','Updated'])
        return 0
    psh = proj_sub.add_parser('show', help='Show current or specific project')
    psh.add_argument('--root', help='Repo root (default: auto)')
    psh.add_argument('--id', help='Project id (default: last)')
    psh.set_defaults(func=_cmd_proj_show)

    # ---- start (interactive wizard) ----
    p = sub.add_parser('start', help='Interactive setup: pick Terminal, prepare EA, generate preset, build, deploy, deploy agents')
    p.add_argument('--project-id', help='Saved project id to use (default: last saved)')
    p.add_argument('--project', help='Project/preset name (default: GPU_LegacyWave1.0.4)')
    p.add_argument('--ea-name', default='CommandListener', help='EA name without extension (default: CommandListener)')
    p.add_argument('--libs', help='Target MQL5\\Libraries path (skips selection if provided)')
    p.add_argument('--metaeditor', help='Full path to metaeditor64.exe (skips auto-detect/prompt)')
    p.add_argument('--yes', action='store_true', help='Assume yes to prompts when safe')
    def cmd_start(args):
        # 1) Project name
        root_cfg = repo_root_from_here()
        pdata = load_projects(root_cfg)
        chosen_id = args.project_id or pdata.get('last_project')
        project = None
        default_name = 'GPU_LegacyWave1.0.4'
        if chosen_id and chosen_id in pdata.get('projects', {}):
            projinfo = pdata['projects'][chosen_id]
            project = projinfo.get('project') or default_name
            args.libs = args.libs or projinfo.get('libs')
            args.metaeditor = args.metaeditor or projinfo.get('metaeditor')
            print(f"Using saved project '{chosen_id}': name={project}")
        project = project or args.project or (input(f'Project name [{default_name}]: ').strip() or default_name)

        # 2) Pick Terminal (MQL5\Libraries)
        if args.libs:
            libs = args.libs
            if libs.startswith('~'):
                libs = os.path.expanduser(libs)
            if not os.path.isdir(libs):
                print('ERROR: invalid Libraries path:', libs)
                return 2
            terminal = {'id':'manual','root': os.path.dirname(os.path.dirname(libs)), 'libs': libs}
        else:
            detected = terminal_dirs()
            print('\nTerminals detected:')
            for i, d in enumerate(detected, 1):
                print(f" {i}. id={d['id']} libs={d['libs']}")
            print(f" {len(detected)+1}. Other... (enter a custom MQL5\\Libraries path)")
            choice = input(f'Select [1-{len(detected)+1}]: ').strip()
            try:
                idx = int(choice)
            except Exception:
                idx = len(detected)+1
            if idx==len(detected)+1:
                libs = input('Enter full path to MQL5\\Libraries: ').strip()
                if libs.startswith('~'):
                    libs = os.path.expanduser(libs)
                if not os.path.isdir(libs):
                    print('ERROR: invalid Libraries path. Aborting.')
                    return 2
                terminal = {'id':'manual','root': os.path.dirname(os.path.dirname(libs)), 'libs': libs}
            else:
                terminal = detected[idx-1]
        mql5_root = os.path.dirname(terminal['libs'])
        experts_dir = os.path.join(mql5_root, 'Experts', project)
        os.makedirs(experts_dir, exist_ok=True)

        # 3) Prepare and compile EA
        ea_base = args.ea_name
        ea_mq5 = os.path.join(experts_dir, ea_base + '.mq5')
        if not os.path.isfile(ea_mq5):
            # minimal EA template
            tpl = (
                '#property strict\n'
                f'// Auto-generated by mtcli start for project {project}\n'
                'int OnInit(){ Print("CommandListener init"); return(INIT_SUCCEEDED);}\n'
                'void OnTick(){}\n'
                'void OnChartEvent(const int id,const long &l,const double &d,const string &s){ PrintFormat("EVT %d %s", id, s); }\n'
            )
            with open(ea_mq5, 'w', encoding='utf-8') as f:
                f.write(tpl)
            print('Created EA:', ea_mq5)

        # locate MetaEditor
        meta = args.metaeditor or os.environ.get('METAEDITOR','')
        if not meta:
            # try powershell
            try:
                r = subprocess.run(['powershell','-NoProfile','-Command', "(Get-Command metaeditor64.exe -ErrorAction SilentlyContinue).Source"], capture_output=True, text=True)
                cand = (r.stdout or '').strip()
                if not cand:
                    r = subprocess.run(['powershell','-NoProfile','-Command', "Get-ChildItem 'C:\\Program Files*' -Recurse -Filter metaeditor64.exe -ErrorAction SilentlyContinue | Select-Object -First 1 -ExpandProperty FullName"], capture_output=True, text=True)
                    cand = (r.stdout or '').strip()
                meta = cand
            except Exception:
                meta = ''
        def _meta_exists(p):
            if not p:
                return False
            if os.path.isfile(p):
                return True
            if is_wsl():
                try:
                    r = subprocess.run(['powershell','-NoProfile','-Command', f"if(Test-Path '{p}'){ 'OK' }"], capture_output=True, text=True)
                    return 'OK' in (r.stdout or '')
                except Exception:
                    return False
            return False
        if args.metaeditor:
            pass
        elif not _meta_exists(meta):
            print('WARN: MetaEditor not found automatically. Provide full path (ex.: C:\\Program Files\\MetaTrader 5\\metaeditor64.exe)')
            meta = input('MetaEditor path: ').strip()
        if not (args.metaeditor or _meta_exists(meta)):
            print('ERROR: MetaEditor path invalid. Aborting.')
            return 3

        # compile EA
        log_dir = os.path.join(repo_root_from_here(), 'tools', 'mtcli_build')
        os.makedirs(log_dir, exist_ok=True)
        log_path = os.path.join(log_dir, f'compile_{datetime.datetime.now().strftime("%Y%m%d_%H%M%S")}.log')
        # Use powershell for consistency
        ps = 'powershell.exe' if is_wsl() else 'powershell'
        meta_p = to_winpath(meta)
        log_p = to_winpath(log_path)
        ea_p = to_winpath(ea_mq5)
        cmd = [ps,'-NoProfile','-Command', f'& "{meta_p}" /log:"{log_p}" /compile:"{ea_p}"']
        run(cmd, check=False)
        ea_ex5 = os.path.join(experts_dir, ea_base + '.ex5')
        if not os.path.isfile(ea_ex5):
            print('ERROR: compile failed. See log:', log_path)
            return 4
        print('EA compiled OK:', ea_ex5)

        # create template and config
        tpl_dir = os.path.join(mql5_root, 'Profiles', 'Templates')
        os.makedirs(tpl_dir, exist_ok=True)
        tpl_path = os.path.join(tpl_dir, project + '.tpl')
        if not os.path.isfile(tpl_path):
            with open(tpl_path, 'w', encoding='utf-8') as f:
                f.write('; template placeholder for '+project+'\n')
        cfg_dir = os.path.join(terminal['root'], 'config')
        os.makedirs(cfg_dir, exist_ok=True)
        cfg_path = os.path.join(cfg_dir, f'terminal_{project}.ini')
        if not os.path.isfile(cfg_path):
            with open(cfg_path, 'w', encoding='utf-8') as f:
                f.write('[Common]\nEnableAutoTrading=1\nAllowDllImport=1\n')
        print('Template and config prepared.')

        # 4) Generate preset via opgen
        root = repo_root_from_here()
        opgen = os.path.join(root, 'tools', 'bridge_gen', 'opgen.py')
        if os.path.isfile(opgen):
            run([sys.executable, opgen, '--root', root, '--preset', project, '--use-preset-selection'])
        else:
            print('WARN: opgen.py not found; skipping preset regeneration.')

        # 5) Build (skip in WSL; use existing DLL)
        if is_wsl():
            print('WSL environment: skipping CMake build (use existing dist-wave/alglib_bridge.dll).')
        else:
            bargs = argparse.Namespace(root=root, build_dir=None, arch='x64', config='Release', targets=['alglib_bridge'], parallel=os.cpu_count() or 8, preset=project, use_preset_selection=True, no_union_presets=False, extra=[])
            cmd_build(bargs)

        # 6) Deploy bridge to selected Terminal only
        dw = dist_wave_dir(root)
        bridge = os.path.join(dw, 'alglib_bridge.dll')
        dargs = argparse.Namespace(root=root, files=[bridge], terminal_id=None, all=False, libs=[terminal['libs']], yes=True, no_kill=False, trash=False)
        cmd_deploy(dargs)

        # 7) Deploy agents (copy service)
        service = os.path.join(dw, 'alglib_service.exe')
        if os.path.isfile(service):
            dest_service_dir = os.path.join(terminal['root'], 'alglib_service')
            os.makedirs(dest_service_dir, exist_ok=True)
            shutil.copy2(service, os.path.join(dest_service_dir, 'alglib_service.exe'))
            print('Service copied to:', dest_service_dir)
            if args.yes or (input('Start service now? [y/N] ').strip().lower() in ('y','yes')):
                try:
                    run(['powershell','-NoProfile','-Command', f'Start-Process -WindowStyle Hidden "{os.path.join(dest_service_dir, "alglib_service.exe")}"'], check=False)
                except Exception:
                    pass
        else:
            print('WARN: service not found; skip agents deploy.')

        # persist as last_project
        pdata = load_projects(root_cfg)
        # if we used a saved id, update it, else save under derived id
        pid = chosen_id or (project.lower().replace(' ', '-'))
        pdata.setdefault('projects', {}).setdefault(pid, {})
        pdata['projects'][pid].update({'project': project, 'libs': terminal['libs'], 'metaeditor': meta, 'updated_at': datetime.datetime.now().isoformat(timespec='seconds')})
        pdata['last_project'] = pid
        save_projects(root_cfg, pdata)
        print('\nStart sequence completed for project:', project, 'as id:', pid)
        return 0
    p.set_defaults(func=cmd_start)

    args = ap.parse_args()
    return args.func(args)


if __name__ == '__main__':
    sys.exit(main())
