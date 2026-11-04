"""Execute sample notebooks and check that no cell raised an error (slow: run with ``pytest -m notebooks``)."""
import os
import subprocess
import sys

import nbformat
import pytest

pytestmark = pytest.mark.notebooks

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

NOTEBOOKS = [
    ('samples', 'beamlines', 'S1_SPB_CRL_simplified', 'S1_SPB_CRL_simplified.ipynb'),
    ('samples', 'Tutorials', 'Tutorial_intro.ipynb'),
    ('samples', 'Tutorials', 'Tutorial_case_1.ipynb'),
    ('samples', 'Tutorials', 'Tutorial_case_2.ipynb'),
    ('samples', 'Tutorials', 'Tutorial_case_2_new.ipynb'),
    ('samples', 'Tutorials', 'Tutorial_case_3.ipynb'),
    ('samples', 'Tutorials', 'Tutorial_case_3_new.ipynb'),
]


def _notebook_run(path, out_dir):
    """Execute a notebook via nbconvert in its own folder (notebooks use relative data paths).

    :returns (executed nb object, execution errors)
    """
    out_path = os.path.join(out_dir, os.path.basename(path))
    args = [sys.executable, '-m', 'nbconvert', '--execute', '--allow-errors',
            '--ExecutePreprocessor.timeout=-1', '--to', 'notebook',
            '--output', out_path, path]
    env = dict(os.environ, MPLBACKEND='Agg')
    subprocess.check_call(args, cwd=os.path.dirname(path), env=env)

    nb = nbformat.read(out_path, nbformat.current_nbformat)
    errors = [output for cell in nb.cells if "outputs" in cell
              for output in cell["outputs"]
              if output.output_type == "error"]
    return nb, errors


@pytest.mark.parametrize('parts', NOTEBOOKS, ids=lambda parts: parts[-1])
def test_notebook(parts, tmp_path):
    _, errors = _notebook_run(os.path.join(ROOT, *parts), str(tmp_path))
    assert errors == [], '\n'.join(f"{e.ename}: {e.evalue}" for e in errors)
