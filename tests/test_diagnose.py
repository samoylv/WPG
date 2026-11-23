"""python -m wpg_diagnose: report for bug reports."""
import json
import os
import subprocess
import sys

import pytest

import wpg_diagnose

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_markdown_sections():
    text = wpg_diagnose.report()
    for title in ('### Problems found', '### Python', '### WPG', '### SRW (srwpy)', '### OpenMP and CUDA libraries',
                  '### GPU', '### Dependencies', '### Environment variables'):
        assert title in text
    assert 'import: ok' in text


def test_json_report():
    rep = json.loads(wpg_diagnose.report(as_json=True))
    assert rep['srw']['import'] == 'ok'
    assert rep['srw']['version']
    assert rep['srw']['reproduce']
    assert rep['wpg']['version']
    assert isinstance(rep['warnings'], list)


def test_home_directory_hidden():
    home = os.path.expanduser('~')
    text = wpg_diagnose.report()
    assert home not in text and home.replace('\\', '/') not in text
    assert home in wpg_diagnose.report(full_paths=True) or home not in sys.executable


def test_environment_filtered(monkeypatch):
    monkeypatch.setenv('OMP_NUM_THREADS', '3')
    monkeypatch.setenv('MY_SECRET_TOKEN', 'abc')
    env = wpg_diagnose.environment()
    assert env['OMP_NUM_THREADS'] == '3'
    assert 'MY_SECRET_TOKEN' not in env


def test_self_test_cpu():
    rep = wpg_diagnose.collect(test=True)
    cpu = rep['self_test'][0]
    assert cpu['test'].startswith('CPU'), cpu
    assert cpu['ok'], cpu


def test_build_info_of_another_binary_is_reported(tmp_path, monkeypatch):
    """A description left from another build (reinstalled srwpy) must not be taken for this binary."""
    import srwpy.srwlpy as s
    info = tmp_path / wpg_diagnose.BUILD_INFO_NAME
    info.write_text(json.dumps({'sha256': '0' * 64, 'mode': 'omp', 'srwlpy': 'elsewhere', 'date': 'then'}))
    monkeypatch.setattr(wpg_diagnose, '_build_info_candidates', lambda binary: [str(info)])
    build, stale = wpg_diagnose._build_info(s.__file__, wpg_diagnose._sha256(s.__file__))
    assert build is None
    assert any('describes another binary' in item for item in stale)


def test_build_info_of_this_binary_gives_reproduce_command(tmp_path, monkeypatch):
    import srwpy.srwlpy as s
    info = tmp_path / wpg_diagnose.BUILD_INFO_NAME
    info.write_text(json.dumps({'sha256': wpg_diagnose._sha256(s.__file__), 'mode': 'omp', 'commit': 'abc123',
                                'repo': wpg_diagnose.DEFAULT_SRW_REPO, 'ref': 'master'}))
    monkeypatch.setattr(wpg_diagnose, '_build_info_candidates', lambda binary: [str(info)])
    srw = wpg_diagnose.srw_info()
    assert srw['origin'].startswith('local build')
    assert srw['reproduce'] == 'make install SRW_MODE=omp SRW_REF=abc123'


@pytest.mark.parametrize('blocked', ['srwpy'])
def test_works_when_srwpy_does_not_import(blocked):
    """The usual problem to diagnose: srwpy is missing or broken, import wpg fails."""
    code = (f'import sys; sys.modules[{blocked!r}] = None\n'
            'import wpg_diagnose\n'
            'print(wpg_diagnose.report())\n'
            'try:\n'
            '    import wpg\n'
            'except ImportError as e:\n'
            '    print("WPG_IMPORT_ERROR", e)\n')
    res = subprocess.run([sys.executable, '-X', 'utf8', '-c', code], capture_output=True, text=True, cwd=ROOT,
                         timeout=120, encoding='utf-8')
    assert res.returncode == 0, res.stderr
    assert '### SRW (srwpy)' in res.stdout
    assert 'wpg does not import' in res.stdout
    assert 'python -m wpg_diagnose' in res.stdout.split('WPG_IMPORT_ERROR', 1)[1]


def test_command_line():
    res = subprocess.run([sys.executable, '-X', 'utf8', '-m', 'wpg_diagnose', '--json'], capture_output=True,
                         text=True, cwd=ROOT, timeout=120, encoding='utf-8')
    assert res.returncode == 0, res.stderr
    assert json.loads(res.stdout)['command'] == 'python -m wpg_diagnose --json'


def _rep(**srw):
    return {'python': wpg_diagnose.python_info(), 'srw': srw, 'openmp_runtime': {}}


def test_activated_environment_is_not_this_python(monkeypatch, tmp_path):
    """Maxwell: VIRTUAL_ENV set, but module load put the shared Python earlier in PATH."""
    monkeypatch.setenv('VIRTUAL_ENV', str(tmp_path / 'venv'))
    w = wpg_diagnose.interpreter_problems(_rep())
    assert any('VIRTUAL_ENV=' in item and 'is activated, but this Python is' in item for item in w), w
    monkeypatch.setenv('VIRTUAL_ENV', sys.prefix)
    assert not any('VIRTUAL_ENV=' in item for item in wpg_diagnose.interpreter_problems(_rep()))


def test_srwpy_in_user_site(monkeypatch, tmp_path):
    """pip fell back to ~/.local because the shared Python is not writable."""
    monkeypatch.delenv('VIRTUAL_ENV', raising=False)
    user_site = tmp_path / 'local' / 'site-packages'
    rep = _rep(srwlpy=str(user_site / 'srwpy' / 'srwlpy.so'))
    rep['python']['user_site_packages'] = str(user_site)
    w = wpg_diagnose.interpreter_problems(rep)
    assert any(item.startswith('srwpy is installed in the user site-packages') for item in w), w
    rep['python']['user_site_packages'] = 'disabled'
    assert not any('user site-packages' in item for item in wpg_diagnose.interpreter_problems(rep))


def test_openmp_runtime_from_another_prefix(monkeypatch, tmp_path):
    monkeypatch.delenv('VIRTUAL_ENV', raising=False)
    rep = _rep()
    rep['openmp_runtime'] = {'loaded_libraries': [str(tmp_path / 'mamba' / 'lib' / 'libgomp.so.1.0.0'),
                                                  os.path.join(sys.prefix, 'lib', 'libgomp.so.1')]}
    w = [item for item in wpg_diagnose.interpreter_problems(rep) if 'OpenMP runtime' in item]
    assert len(w) == 1 and str(tmp_path) in w[0], w


def test_wpg_imported_from_another_directory(monkeypatch, tmp_path):
    monkeypatch.setattr(wpg_diagnose, '_installed_wpg_dir', lambda: str(tmp_path / 'wpg'))
    found = wpg_diagnose.wpg_copies(wpg_diagnose.wpg_info())
    assert any('pip installed WPG in' in f for f in found), found
