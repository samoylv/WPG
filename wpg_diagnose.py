"""Diagnostics of a WPG installation for bug reports.

Run and paste the output into a GitHub issue:

    python -m wpg_diagnose            # environment, WPG, SRW build, OpenMP, GPU, dependencies (about 1 s)
    python -m wpg_diagnose --test     # plus a numerical self-test of SRW on CPU, OpenMP and GPU (about 10 s)
    python -m wpg_diagnose --json     # the same as JSON

This module is deliberately outside the ``wpg`` package and imports neither ``wpg`` nor ``srwpy`` at
module level: it must work when ``import wpg`` fails. Every section catches its own errors and reports them.
The home directory in paths is replaced by ``~`` unless ``--full-paths`` is given; the host name and
the full environment are never printed, only variables that change SRW behaviour.
"""
import argparse
import datetime
import hashlib
import importlib.metadata as metadata
import importlib.util
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
import traceback

MIN_SRWPY = (4, 2)
BUILD_INFO_NAME = 'srw_build_info.json'
DEFAULT_SRW_REPO = 'https://github.com/ochubar/SRW.git'
ISSUES_URL = 'https://github.com/samoylv/WPG/issues'
ENV_PREFIXES = ('OMP_', 'KMP_', 'GOMP_', 'CUDA_', 'MKL_', 'OPENBLAS_', 'NVIDIA_')
ENV_NAMES = ('PYTHONPATH', 'PYTHONHOME', 'MPLBACKEND', 'LD_LIBRARY_PATH', 'LD_PRELOAD', 'DYLD_LIBRARY_PATH',
             'CONDA_DEFAULT_ENV', 'VIRTUAL_ENV')
DEPENDENCIES = ('numpy', 'scipy', 'h5py', 'matplotlib', 'xraylib', 'phenom-xfel', 'ocelot-collab',
                'threadpoolctl', 'pytest')
# markers of the OpenMP and CUDA runtimes in the srwlpy binary (Linux: GNU and LLVM OpenMP, Windows: MSVC)
OPENMP_MARKERS = (b'GOMP_parallel', b'libgomp', b'__kmpc_fork_call', b'libomp', b'VCOMP140.DLL', b'vcomp140.dll')
CUDA_MARKERS = (b'cufftPlan', b'libcufft', b'cudaLaunchKernel', b'cudart')


# ------------------------------------------------------------------ helpers

def _version_tuple(text):
    return tuple(int(v) for v in re.findall(r'\d+', text or '')[:3])


def _dist(name):
    try:
        return metadata.distribution(name)
    except metadata.PackageNotFoundError:
        return None


def _direct_url(dist):
    """PEP 610 record written by pip: where the package was installed from (None for a package index)."""
    try:
        text = dist.read_text('direct_url.json')
    except Exception:
        return None
    return json.loads(text) if text else None


def _run(cmd, cwd=None, timeout=10):
    try:
        res = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)
    except (OSError, subprocess.SubprocessError):
        return None
    return res.stdout.strip() if res.returncode == 0 else None


def _sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def _error(exc):
    return f'{type(exc).__name__}: {exc}'


# ------------------------------------------------------------------ sections

def python_info():
    prefix = sys.prefix
    if os.path.isdir(os.path.join(prefix, 'conda-meta')):
        kind = f"conda, environment '{os.environ.get('CONDA_DEFAULT_ENV') or os.path.basename(prefix)}'"
    elif prefix != sys.base_prefix:
        kind = 'venv'
        cfg = os.path.join(prefix, 'pyvenv.cfg')
        if os.path.exists(cfg) and 'uv' in open(cfg, encoding='utf-8', errors='replace').read():
            kind = 'venv (uv)'
    else:
        kind = 'system or base interpreter'
    in_notebook = 'ipykernel' in sys.modules
    return {
        'version': sys.version.split()[0],
        'implementation': platform.python_implementation(),
        'executable': sys.executable,
        'environment': kind,
        'prefix': prefix,
        'platform': f'{platform.system()} {platform.release()} ({platform.version()})',
        'machine': platform.machine(),
        'cpu_count': os.cpu_count(),
        'jupyter_kernel': in_notebook,
    }


def _wpg_distribution():
    """Installed WPG metadata; WPG.egg-info in the source tree (on sys.path when run from it) is a fallback."""
    spec = importlib.util.find_spec('wpg')
    source = os.path.abspath(os.path.dirname(os.path.dirname(spec.origin))) if spec else None
    dists = [d for d in metadata.distributions() if (d.metadata['Name'] or '').lower() == 'wpg']
    installed = [d for d in dists if os.path.abspath(str(d.locate_file(''))) != source]
    return (installed or dists or [None])[0]


def wpg_info():
    info = {}
    dist = _wpg_distribution()
    info['version'] = dist.version if dist else None
    url = _direct_url(dist) if dist else None
    if url:
        info['installed_from'] = url.get('url')
        info['editable'] = bool(url.get('dir_info', {}).get('editable'))
    spec = importlib.util.find_spec('wpg')
    if spec is None:
        info['location'] = None
        return info
    pkg_dir = os.path.dirname(spec.origin)
    info['location'] = pkg_dir
    src = os.path.dirname(pkg_dir)
    if os.path.isdir(os.path.join(src, '.git')) or _run(['git', '-C', src, 'rev-parse', '--git-dir']):
        info['git_commit'] = _run(['git', '-C', src, 'rev-parse', 'HEAD'])
        info['git_branch'] = _run(['git', '-C', src, 'rev-parse', '--abbrev-ref', 'HEAD'])
        status = _run(['git', '-C', src, 'status', '--porcelain', '--untracked-files=no'])
        info['git_modified_files'] = len(status.splitlines()) if status else 0
    try:
        import wpg  # noqa: F401  (imports srwpy)
        info['import'] = 'ok'
    except Exception as exc:
        info['import'] = _error(exc)
        info['import_traceback'] = traceback.format_exc()
    return info


def _installed_wpg_dir():
    """Directory of the wpg package that pip installed (editable: in the source tree it points to)."""
    from urllib.parse import urlparse
    from urllib.request import url2pathname
    dist = _wpg_distribution()
    if dist is None:
        return None
    url = _direct_url(dist)
    if url and url.get('dir_info', {}).get('editable') and url.get('url', '').startswith('file:'):
        return os.path.join(url2pathname(urlparse(url['url']).path), 'wpg')
    return str(dist.locate_file('wpg'))


def wpg_copies(wpg):
    """Other WPG and SRW copies that can shadow the installed ones."""
    found = []
    here = wpg.get('location')
    installed = _installed_wpg_dir()
    if here and installed and os.path.normcase(os.path.abspath(installed)) != os.path.normcase(os.path.abspath(here)):
        found.append(f'wpg is imported from {here}, but pip installed WPG in {installed} '
                     '(PYTHONPATH or the current directory?)')
    source = os.path.abspath(os.path.dirname(here)) if here else None
    # WPG.egg-info in the source tree of an editable install is the same installation, not a copy
    dists = [d for d in metadata.distributions() if (d.metadata['Name'] or '').lower() == 'wpg'
             and os.path.abspath(str(d.locate_file(''))) != source]
    if len(dists) > 1:
        found.append(f'{len(dists)} installed distributions named WPG: '
                     + ', '.join(str(d.locate_file('')) for d in dists))
    for entry in sys.path:
        candidate = os.path.join(entry or os.getcwd(), 'wpg', '__init__.py')
        if here and os.path.exists(candidate) and os.path.dirname(os.path.abspath(candidate)) != os.path.abspath(here):
            found.append(f'another wpg package on sys.path (not imported): {os.path.dirname(candidate)}')
    if here:
        old_srw = os.path.join(here, 'srw')
        if os.path.isdir(old_srw):
            names = sorted(n for n in os.listdir(old_srw) if not n.startswith('__'))
            found.append(f'leftover of the bundled SRW of WPG before 2026.11 in {old_srw}: {", ".join(names) or "empty"}')
    for name in ('srwlpy', 'srwlib'):
        try:
            spec = importlib.util.find_spec(name)
        except (ImportError, ValueError):
            spec = None
        if spec is not None:
            found.append(f'top-level module {name} (an old SRW installation, not srwpy): {spec.origin}')
    return found


def _wheel_tags(dist):
    try:
        text = dist.read_text('WHEEL') or ''
    except Exception:
        return []
    return [line.split(':', 1)[1].strip() for line in text.splitlines() if line.startswith('Tag:')]


def _build_info_candidates(binary):
    """Files written by build/build_srw.py: next to srwlpy and in build/ of the WPG checkout (last build only)."""
    candidates = [os.path.join(os.path.dirname(binary), BUILD_INFO_NAME)]
    spec = importlib.util.find_spec('wpg')
    if spec is not None:
        candidates.append(os.path.join(os.path.dirname(os.path.dirname(spec.origin)), 'build', BUILD_INFO_NAME))
    return candidates


def _build_info(binary, digest):
    """Description of a local SRW build written by build/build_srw.py, if it belongs to this binary."""
    stale = []
    for path in _build_info_candidates(binary):
        if not os.path.exists(path):
            continue
        try:
            with open(path, encoding='utf-8') as f:
                info = json.load(f)
        except (OSError, ValueError) as exc:
            stale.append(f'{path}: {_error(exc)}')
            continue
        if info.get('sha256') == digest:
            info['file'] = path
            return info, stale
        stale.append(f'{path} describes another binary ({info.get("srwlpy")}, mode {info.get("mode")}, '
                     f'{info.get("date")})')
    return None, stale


def srw_info():
    info = {}
    dist = _dist('srwpy')
    info['version'] = dist.version if dist else None
    if dist is None:
        info['import'] = 'srwpy is not installed'
        return info
    try:
        spec = importlib.util.find_spec('srwpy.srwlpy')
    except Exception as exc:
        spec = None
        info['import'] = _error(exc)
    if spec is None or not spec.origin:
        info.setdefault('import', 'srwpy.srwlpy not found')
        return info
    binary = spec.origin
    info['srwlpy'] = binary
    info['srwlpy_size_mb'] = round(os.path.getsize(binary) / 2 ** 20, 1)
    info['srwlpy_date'] = datetime.datetime.fromtimestamp(os.path.getmtime(binary)).isoformat(timespec='seconds')
    with open(binary, 'rb') as f:
        data = f.read()
    info['openmp'] = any(m in data for m in OPENMP_MARKERS)
    info['cuda'] = any(m in data for m in CUDA_MARKERS)
    digest = _sha256(binary)
    info['sha256'] = digest

    url = _direct_url(dist)
    tags = _wheel_tags(dist)
    build, stale = _build_info(binary, digest)
    if stale:
        info['build_info_not_matching'] = stale
    if build:
        info['origin'] = 'local build with build/build_srw.py'
        info['build'] = {k: build.get(k) for k in ('mode', 'repo', 'ref', 'commit', 'src', 'date', 'file')}
        ref = build.get('commit') or build.get('ref') or 'master'
        repo = build.get('repo')
        if repo and repo != DEFAULT_SRW_REPO:
            info['reproduce'] = (f'python build/build_srw.py --mode {build.get("mode")} --repo {repo} --ref {ref}'
                                 ' && python -m pip install -e .')
        elif repo:
            info['reproduce'] = f'make install SRW_MODE={build.get("mode")} SRW_REF={ref}'
        else:
            info['reproduce'] = (f'python build/build_srw.py --mode {build.get("mode")} --src <SRW checkout at '
                                 f'{build.get("commit")}>')
    elif url is None:
        info['origin'] = 'package index (PyPI wheel)'
        info['reproduce'] = f'python -m pip install srwpy=={dist.version}'
    else:
        info['origin'] = f'local install from {url.get("url")}'
        info['reproduce'] = 'unknown: no build description next to the binary (built without build/build_srw.py?)'
    info['wheel_tags'] = tags
    try:
        from srwpy import srwlpy  # noqa: F401
        info['import'] = 'ok'
    except Exception as exc:
        info['import'] = _error(exc)
        info['import_traceback'] = traceback.format_exc()
    return info


def openmp_runtime():
    """OpenMP libraries loaded into this process after importing srwlpy, with their thread numbers."""
    out = {}
    try:
        from threadpoolctl import threadpool_info
        # BLAS pools of numpy and scipy do not concern SRW, only OpenMP runtimes are listed
        out['threadpools'] = [f"{p.get('internal_api')} {p.get('version')}, {p.get('num_threads')} threads: "
                              f"{p.get('filepath')}" for p in threadpool_info() if p.get('user_api') == 'openmp']
    except ImportError:
        out['threadpools'] = 'threadpoolctl is not installed'
    except Exception as exc:
        out['threadpools'] = _error(exc)
    if os.path.exists('/proc/self/maps'):
        libs = set()
        with open('/proc/self/maps') as f:
            for line in f:
                path = line.split()[-1]
                if re.search(r'lib(g|i)?omp|libcufft|libcudart', os.path.basename(path)):
                    libs.add(path)
        out['loaded_libraries'] = sorted(libs)
    return out


def gpu_info():
    smi = shutil.which('nvidia-smi')
    if not smi:
        return {'nvidia_smi': 'not found'}
    out = {}
    query = _run([smi, '--query-gpu=name,driver_version,memory.total', '--format=csv,noheader'])
    out['gpus'] = query.splitlines() if query else 'nvidia-smi failed'
    head = _run([smi]) or ''
    # 'CUDA Version: 12.4' in older drivers, 'CUDA UMD Version: 13.4' in newer ones
    m = re.search(r'CUDA (?:UMD )?Version:\s*([\d.]+)', head)
    out['driver_cuda_version'] = m.group(1) if m else None
    return out


def dependencies():
    out = {}
    for name in DEPENDENCIES:
        d = _dist(name)
        out[name] = d.version if d else None
    try:
        import h5py
        out['hdf5'] = h5py.version.hdf5_version
    except Exception as exc:
        out['hdf5'] = _error(exc)
    try:
        import matplotlib
        out['matplotlib_backend'] = matplotlib.get_backend()
    except Exception as exc:
        out['matplotlib_backend'] = _error(exc)
    return out


def environment():
    return {k: v for k, v in sorted(os.environ.items()) if k.startswith(ENV_PREFIXES) or k in ENV_NAMES}


# ------------------------------------------------------------------ self-test

def _gauss(n, half, sig, d2waist, nz=1):
    """Gaussian beam of 9 keV: one slice in frequency domain or nz slices of a 2 fs pulse in time domain."""
    from wpg import Wavefront
    from wpg.generators import build_gauss_wavefront, build_gauss_wavefront_xy
    if nz == 1:
        return Wavefront(build_gauss_wavefront_xy(n, n, 9., -half, half, -half, half, sig, sig, d2waist))
    return Wavefront(build_gauss_wavefront(n, n, nz, 9., -half, half, -half, half, 2e-15, sig, sig, d2waist))


def _propagate(wf, length, semi, gpu=None):
    from wpg import Beamline
    from wpg.optical_elements import Drift, Use_PP
    bl = Beamline()
    bl.append(Drift(length), Use_PP(semi_analytical_treatment=semi))
    t0 = time.perf_counter()
    bl.propagate(wf, gpu=gpu)
    return time.perf_counter() - t0


def _rms_x(wf):
    import numpy as np
    m = wf.params.Mesh
    x = np.linspace(m.xMin, m.xMax, m.nx)
    p = np.asarray(wf.get_intensity(), dtype=float).sum(axis=(0, 2))  # intensity is (ny, nx, slices)
    p = p / p.sum()
    c = (p * x).sum()
    return float(np.sqrt((p * (x - c) ** 2).sum()))


def _bench_slices(threads):
    """Time of a drift of 32 slices (frequency domain), run in a child process with OMP_NUM_THREADS set."""
    import numpy as np
    from srwpy.srwlib import srwl
    wf = _gauss(256, 7e-4, 5e-6, 50., nz=32)
    srwl.SetRepresElecField(wf._srwl_wf, 'f')
    dt = _propagate(wf, 20., 1)
    total = float(np.asarray(wf.get_intensity(), dtype=float).sum())
    return {'threads': threads, 'seconds': round(dt, 3), 'total_intensity': total}


def _child_bench(threads):
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    res = subprocess.run([sys.executable, '-m', 'wpg_diagnose', '--_bench', str(threads)], env=env,
                         capture_output=True, text=True, timeout=300,
                         cwd=os.path.dirname(os.path.abspath(__file__)))
    if res.returncode != 0:
        raise RuntimeError(res.stderr.strip().splitlines()[-1] if res.stderr.strip() else 'child failed')
    return json.loads(res.stdout.strip().splitlines()[-1])


def self_test(srw):
    """Drift of a Gaussian beam against the analytic size, OpenMP speed-up and GPU against CPU."""
    import numpy as np
    out = []
    # 9 keV, waist 5 um RMS 50 m upstream (Rayleigh length 2.28 m, so the SRW radius of curvature is exact
    # to 0.2 %), window +-0.7 mm = 6.4 RMS sizes, drift of 20 m with the semi-analytical treatment
    lam = 1.239841984e-6 / 9e3
    sig0, z0, length = 5e-6, 50., 20.
    z_r = 4 * np.pi * sig0 ** 2 / lam
    expected = sig0 * np.sqrt(1 + ((z0 + length) / z_r) ** 2)
    try:
        wf = _gauss(256, 7e-4, sig0, z0)
        dt = _propagate(wf, length, 1)
        rms = _rms_x(wf)
        dev = rms / expected - 1
        out.append({'test': 'CPU: 20 m drift of a 9 keV Gaussian beam, RMS size against the analytic one',
                    'ok': abs(dev) < 0.02, 'rms_um': round(rms * 1e6, 2), 'analytic_um': round(expected * 1e6, 2),
                    'deviation_percent': round(100 * dev, 2), 'seconds': round(dt, 3)})
    except Exception as exc:
        out.append({'test': 'CPU drift', 'ok': False, 'error': _error(exc), 'traceback': traceback.format_exc()})

    if srw.get('openmp'):
        try:
            one = _child_bench(1)
            many = _child_bench(os.cpu_count())
            same = abs(many['total_intensity'] / one['total_intensity'] - 1) < 1e-6
            out.append({'test': 'OpenMP: 32 slices, 1 thread against all threads (same result, speed-up)',
                        'ok': same, 'seconds_1_thread': one['seconds'],
                        f'seconds_{many["threads"]}_threads': many['seconds'],
                        'speed_up': round(one['seconds'] / max(many['seconds'], 1e-9), 2)})
        except Exception as exc:
            out.append({'test': 'OpenMP speed-up', 'ok': False, 'error': _error(exc)})

    if srw.get('cuda'):
        # a single slice agrees with CPU on SRW 4.2.1, several slices without the semi-analytical treatment do not
        # (known SRW error: the drift length is ignored on GPU), see tests/test_wpg.py::test_gpu_option
        for label, nz, semi, known in (('1 slice', 1, 0, False),
                                       ('3 slices, no semi-analytical treatment', 3, 0, True)):
            try:
                from srwpy.srwlib import srwl
                res, times = [], []
                for gpu in (None, 1):
                    wf = _gauss(200, 7e-4, sig0, z0, nz=nz)
                    if nz > 1:
                        srwl.SetRepresElecField(wf._srwl_wf, 'f')
                    times.append(_propagate(wf, length, semi, gpu=gpu))
                    res.append(np.asarray(wf.get_intensity(), dtype=float))
                diff = float(np.abs(res[1] - res[0]).max() / res[0].max())
                out.append({'test': f'GPU against CPU, {label}', 'ok': diff < 1e-4, 'known_srw_error': known,
                            'max_difference_of_peak': f'{diff:.1e}', 'seconds_cpu': round(times[0], 3),
                            'seconds_gpu': round(times[1], 3)})
            except Exception as exc:
                out.append({'test': f'GPU against CPU, {label}', 'ok': False, 'error': _error(exc)})
    return out


# ------------------------------------------------------------------ report

def _warnings(rep):
    w = []
    srw, wpg = rep['srw'], rep['wpg']
    if srw.get('version') is None:
        w.append('srwpy is not installed: python -m pip install srwpy')
    elif _version_tuple(srw['version']) < MIN_SRWPY:
        w.append(f'srwpy {srw["version"]} is older than {".".join(map(str, MIN_SRWPY))} required by WPG')
    if srw.get('import') not in (None, 'ok'):
        w.append(f'srwpy.srwlpy does not import: {srw["import"]}')
    if wpg.get('import') not in (None, 'ok'):
        w.append(f'wpg does not import: {wpg["import"]}')
    if wpg.get('location') is None:
        w.append('wpg is not found on sys.path')
    if wpg.get('git_modified_files'):
        w.append(f'the WPG checkout has {wpg["git_modified_files"]} modified files')
    w.extend(rep['copies'])
    # build/ of a WPG checkout describes only its last build (several environments can share one checkout),
    # a description next to srwlpy that does not match it is left from a previous installation
    next_to_binary = os.path.dirname(srw['srwlpy']) if srw.get('srwlpy') else None
    for item in srw.get('build_info_not_matching', []):
        if next_to_binary and item.startswith(next_to_binary):
            w.append(f'SRW build description left from a previous installation: {item}')
    if srw.get('openmp') and srw.get('cuda'):
        w.append('srwlpy contains both OpenMP and CUDA markers, build/build_srw.py builds one of them')
    if srw.get('openmp') and os.environ.get('OMP_NUM_THREADS') == '1':
        w.append('SRW is built with OpenMP, but OMP_NUM_THREADS=1')
    for t in rep.get('self_test', []):
        if not t.get('ok') and not t.get('known_srw_error'):
            w.append(f'self-test failed: {t["test"]}')
    return w


def collect(test=False):
    rep = {'date': datetime.datetime.now().isoformat(timespec='seconds'),
           'command': ' '.join(['python -m wpg_diagnose'] + sys.argv[1:])}
    for key, func in (('python', python_info), ('wpg', wpg_info), ('srw', srw_info)):
        try:
            rep[key] = func()
        except Exception as exc:
            rep[key] = {'error': _error(exc), 'traceback': traceback.format_exc()}
    for key, func in (('copies', lambda: wpg_copies(rep['wpg'])), ('openmp_runtime', openmp_runtime),
                      ('gpu', gpu_info), ('dependencies', dependencies), ('environment', environment)):
        try:
            rep[key] = func()
        except Exception as exc:
            rep[key] = [f'error: {_error(exc)}'] if key == 'copies' else {'error': _error(exc)}
    if test:
        if rep['wpg'].get('import') == 'ok':
            rep['self_test'] = self_test(rep['srw'])
        else:
            rep['self_test'] = [{'test': 'self-test', 'ok': False, 'error': 'skipped: wpg does not import'}]
    rep['warnings'] = _warnings(rep)
    return rep


def _anonymize(obj, home):
    if isinstance(obj, str):
        return obj.replace(home, '~') if home else obj
    if isinstance(obj, dict):
        return {k: _anonymize(v, home) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_anonymize(v, home) for v in obj]
    return obj


def _anonymizer():
    home = os.path.expanduser('~')
    if not home or home in ('/', '~'):
        return lambda obj: obj
    variants = {home, home.replace('\\', '/'), home.replace('/', '\\')}

    def run(obj):
        for h in sorted(variants, key=len, reverse=True):
            obj = _anonymize(obj, h)
        return obj
    return run


def _md_value(v):
    if isinstance(v, bool):
        return 'yes' if v else 'no'
    if v is None:
        return '—'
    if isinstance(v, (list, dict)):
        return json.dumps(v)
    return ' '.join(str(v).splitlines())  # one line per list item


def to_markdown(rep):
    lines = ['## WPG diagnostics', '', f'Generated {rep["date"]} by `{rep["command"]}`.', '']
    lines += ['### Problems found', '']
    lines += [f'- {_md_value(w)}' for w in rep['warnings']] or ['- none']
    sections = (('Python', 'python'), ('WPG', 'wpg'), ('SRW (srwpy)', 'srw'), ('OpenMP and CUDA libraries', 'openmp_runtime'),
                ('GPU', 'gpu'), ('Dependencies', 'dependencies'), ('Environment variables', 'environment'))
    tracebacks = []
    for title, key in sections:
        lines += ['', f'### {title}', '']
        sec = rep.get(key) or {}
        if not sec:
            lines.append('- none')
        for k, v in sec.items():
            if k.endswith('traceback'):
                tracebacks.append((f'{title}: {k}', v))
                continue
            if k == 'build' and isinstance(v, dict):
                lines.append('- build:')
                lines += [f'  - {bk}: {_md_value(bv)}' for bk, bv in v.items()]
            elif k in ('threadpools', 'loaded_libraries', 'gpus', 'build_info_not_matching') and isinstance(v, list):
                lines.append(f'- {k}:' + ('' if v else ' none'))
                lines += [f'  - {_md_value(i)}' for i in v]
            else:
                lines.append(f'- {k}: {_md_value(v)}')
    if 'self_test' in rep:
        lines += ['', '### Self-test', '']
        for t in rep['self_test']:
            details = ', '.join(f'{k} {v}' for k, v in t.items()
                                if k not in ('test', 'ok', 'traceback', 'known_srw_error'))
            if t.get('known_srw_error'):  # error of SRW 4.2.1 on GPU, PASS here means SRW has fixed it
                status = 'PASS (known SRW error is fixed)' if t['ok'] else 'KNOWN SRW ERROR'
            else:
                status = 'PASS' if t['ok'] else 'FAIL'
            lines.append(f'- {status}: {t["test"]}' + (f' ({details})' if details else ''))
            if 'traceback' in t:
                tracebacks.append((t['test'], t['traceback']))
    for title, tb in tracebacks:
        lines += ['', '<details><summary>Traceback: ' + title + '</summary>', '', '```', tb.rstrip(), '```',
                  '', '</details>']
    lines += ['', f'Paste this text into an issue at {ISSUES_URL}.', '']
    return '\n'.join(lines)


def report(test=False, full_paths=False, as_json=False):
    """Diagnostics as Markdown text (or JSON) ready to paste into a bug report."""
    rep = collect(test=test)
    if not full_paths:
        rep = _anonymizer()(rep)
    return json.dumps(rep, indent=2) if as_json else to_markdown(rep)


def main(argv=None):
    parser = argparse.ArgumentParser(prog='python -m wpg_diagnose',
                                     description='Diagnostics of a WPG installation to paste into a bug report.')
    parser.add_argument('--test', action='store_true',
                        help='also run a numerical self-test of SRW on CPU, OpenMP and GPU (about 10 s)')
    parser.add_argument('--json', action='store_true', help='print JSON instead of Markdown')
    parser.add_argument('--full-paths', action='store_true', help='do not replace the home directory by ~')
    parser.add_argument('--_bench', type=int, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args._bench:
        print(json.dumps(_bench_slices(args._bench)))
        return 0
    text = report(test=args.test, full_paths=args.full_paths, as_json=args.json)
    try:
        print(text)
    except UnicodeEncodeError:  # Windows console without UTF-8
        print(text.encode('ascii', 'replace').decode('ascii'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
