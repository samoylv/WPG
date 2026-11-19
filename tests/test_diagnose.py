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


def test_wpg_imported_from_another_directory(monkeypatch, tmp_path):
    monkeypatch.setattr(wpg_diagnose, '_installed_wpg_dir', lambda: str(tmp_path / 'wpg'))
    found = wpg_diagnose.wpg_copies(wpg_diagnose.wpg_info())
    assert any('pip installed WPG in' in f for f in found), found
