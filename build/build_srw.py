#!/usr/bin/env python3
"""
Build SRW (the ``srwpy`` package) from sources and install it into a Python environment.

The official ``srwpy`` wheels from PyPI are serial builds. This script builds SRW
from https://github.com/ochubar/SRW (or another fork / local checkout) in one of the modes:

* ``cpu``  - serial build (FFTW3), the same as the official wheels but from the latest sources
* ``omp``  - OpenMP build (FFTW2), parallel over wavefront slices
* ``cuda`` - CUDA build (FFTW3 + cuFFT)

OpenMP and CUDA cannot be combined: SRW (cpp/cmake/CMakeLists.txt) uses FFTW2 for
the OpenMP build and FFTW3 with cuFFT for the CUDA build.

Only Linux is supported. On Windows and Mac OS X use the official wheels (``pip install srwpy``).

Examples::

    python build/build_srw.py --mode omp
    python build/build_srw.py --mode cuda --ref master --cuda-archs all
    python build/build_srw.py --mode omp --src ~/src/SRW
"""

import argparse
import datetime
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_REPO = 'https://github.com/ochubar/SRW.git'
DEFAULT_WORKDIR = os.path.join(HERE, 'SRW')
BUILD_INFO = os.path.join(HERE, 'srw_build_info.json')
# CMAKE_CUDA_ARCHITECTURES=native in SRW CMakeLists needs CMake >= 3.24
BUILD_REQUIREMENTS = ['packaging', 'setuptools', 'wheel', 'cmake>=3.24']


def run(cmd, **kwargs):
    print('+', ' '.join(cmd), flush=True)
    subprocess.check_call(cmd, **kwargs)


def output(cmd, **kwargs):
    return subprocess.check_output(cmd, text=True, **kwargs).strip()


def checkout(repo, ref, workdir):
    """Shallow fetch of ``ref`` (branch, tag or commit) into ``workdir``."""
    if not os.path.isdir(os.path.join(workdir, '.git')):
        os.makedirs(workdir, exist_ok=True)
        run(['git', 'init', '-q'], cwd=workdir)
        run(['git', 'remote', 'add', 'origin', repo], cwd=workdir)
    else:
        run(['git', 'remote', 'set-url', 'origin', repo], cwd=workdir)
    run(['git', 'fetch', '-q', '--depth', '1', 'origin', ref], cwd=workdir)
    run(['git', 'checkout', '-q', '--force', '--detach', 'FETCH_HEAD'], cwd=workdir)


def build_env(mode, cuda_archs, jobs, python):
    """Environment variables understood by SRW env/python/setup.py."""
    env = os.environ.copy()
    env.pop('MODE', None)
    cmake_args = env.pop('CMAKE_ARGS', '').split()
    if mode == 'omp':
        env['MODE'] = 'omp'  # setup.py: MODE=omp -> -DUSE_OPENMP=ON
    elif mode == 'cuda':
        cmake_args.append('-DUSE_CUDA=ON')
        if cuda_archs == 'all':
            cmake_args.append('-DCUDA_BUILDALL_ARCHS=ON')
        cuda_path = env.get('CUDA_PATH', '/usr/local/cuda')
        nvcc = os.path.join(cuda_path, 'bin', 'nvcc')
        if not os.path.exists(nvcc):
            sys.exit(f'nvcc not found at {nvcc}, set CUDA_PATH to the CUDA toolkit')
        env['CUDA_PATH'] = cuda_path
        env['CUDACXX'] = nvcc
        env['PATH'] = os.path.join(cuda_path, 'bin') + os.pathsep + env['PATH']
    if cmake_args:
        env['CMAKE_ARGS'] = ' '.join(cmake_args)  # setup.py splits it by single spaces
    env['CMAKE_BUILD_PARALLEL_LEVEL'] = str(jobs)

    # setup.py calls "cmake3" first, then "cmake": put the CMake from the target
    # environment first in PATH under both names (system cmake3 may be too old)
    bin_dir = os.path.dirname(python)
    cmake = os.path.join(bin_dir, 'cmake')
    shim_dir = tempfile.mkdtemp(prefix='srw_cmake_')
    os.symlink(cmake, os.path.join(shim_dir, 'cmake3'))
    env['PATH'] = os.pathsep.join([shim_dir, bin_dir, env['PATH']])
    return env, shim_dir


def check_build(python, mode):
    """Import the installed srwlpy and check that it is linked with OpenMP / CUDA as requested."""
    info = json.loads(output([python, '-c',
                              'import json, srwpy.srwlpy as s, importlib.metadata as m;'
                              'print(json.dumps({"srwlpy": s.__file__, "version": m.version("srwpy")}))'],
                             cwd=tempfile.gettempdir()))
    libs = output(['ldd', info['srwlpy']])
    with open(info['srwlpy'], 'rb') as f:
        binary = f.read()
    info['openmp'] = 'libgomp' in libs or b'GOMP_' in binary
    info['cuda'] = 'libcufft' in libs or b'cufftPlan' in binary
    expected = {'cpu': (False, False), 'omp': (True, False), 'cuda': (False, True)}[mode]
    if (info['openmp'], info['cuda']) != expected:
        sys.exit(f'srwlpy is built with openmp={info["openmp"]}, cuda={info["cuda"]}, '
                 f'expected openmp={expected[0]}, cuda={expected[1]} for mode {mode}')
    return info


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--mode', choices=['cpu', 'omp', 'cuda'], default='cpu', help='SRW build mode (default: cpu)')
    parser.add_argument('--src', help='existing SRW checkout to build (no cloning)')
    parser.add_argument('--repo', default=DEFAULT_REPO, help=f'SRW git repository (default: {DEFAULT_REPO})')
    parser.add_argument('--ref', default='master', help='branch, tag or commit of SRW (default: master)')
    parser.add_argument('--workdir', default=DEFAULT_WORKDIR, help='where to clone SRW (default: build/SRW)')
    parser.add_argument('--cuda-archs', choices=['native', 'all'], default='native',
                        help='CUDA architectures: GPU of this machine or 80 86 89 90 120 (default: native)')
    parser.add_argument('--python', default=sys.executable, help='Python of the target environment')
    parser.add_argument('--jobs', type=int, default=os.cpu_count(), help='parallel compile jobs')
    args = parser.parse_args()

    if platform.system() != 'Linux':
        sys.exit('Building SRW from sources is supported on Linux only, '
                 'on this system install the official wheels: python -m pip install srwpy')
    python = os.path.abspath(args.python)

    if args.src:
        src = os.path.abspath(os.path.expanduser(args.src))
    else:
        src = os.path.abspath(args.workdir)
        checkout(args.repo, args.ref, src)
    pkg_dir = os.path.join(src, 'env', 'python')
    if not os.path.exists(os.path.join(pkg_dir, 'setup.py')):
        sys.exit(f'{pkg_dir}/setup.py not found, is {src} an SRW checkout?')
    commit = output(['git', 'rev-parse', 'HEAD'], cwd=src)
    print(f'SRW {commit} from {src}, mode {args.mode}', flush=True)

    run([python, '-m', 'pip', 'install', '-q', '--upgrade'] + BUILD_REQUIREMENTS)

    # CMake cache of a previous build in another mode breaks the build
    old_build = os.path.join(pkg_dir, 'build')
    if os.path.isdir(old_build):
        print(f'removing previous build directory {old_build}', flush=True)
        shutil.rmtree(old_build)

    env, shim_dir = build_env(args.mode, args.cuda_archs, args.jobs, python)
    try:
        run([python, '-m', 'pip', 'install', '--no-build-isolation', '--no-deps', '--force-reinstall',
             '--verbose', pkg_dir], env=env)
    finally:
        shutil.rmtree(shim_dir, ignore_errors=True)

    info = check_build(python, args.mode)
    info.update(mode=args.mode, commit=commit,
                repo=None if args.src else args.repo, ref=None if args.src else args.ref, src=src,
                python=python, date=datetime.datetime.now().isoformat(timespec='seconds'))
    with open(BUILD_INFO, 'w') as f:
        json.dump(info, f, indent=2)
    print(json.dumps(info, indent=2))
    print(f'build info saved to {BUILD_INFO}')


if __name__ == '__main__':
    main()
