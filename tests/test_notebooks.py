import subprocess
import nbformat
import os

# import sys
# sys.path.insert(0, '..')


def _notebook_run(path):
    """Execute a notebook via nbconvert and collect output.
       :returns (parsed nb object, execution errors)
    """

    notebook_dir = os.path.dirname(path)
    test_ipynb = os.path.split(path)[-1] + '.test.ipynb'

    args = ["jupyter", "nbconvert", "--execute", "--allow-errors",
            "--ExecutePreprocessor.timeout=-1",
            "--to", "notebook", '--output', test_ipynb, path]
    subprocess.check_call(args)

    args = ["jupyter", "nbconvert", "--to", "html",
            os.path.join(notebook_dir, test_ipynb)]
    subprocess.check_call(args)

    nb = nbformat.read(path, nbformat.current_nbformat)
    errors = [output for cell in nb.cells if "outputs" in cell
              for output in cell["outputs"]
              if output.output_type == "error"]
    return nb, errors


def _test_notebook(path):
    nb, errors = _notebook_run(path)
    assert errors == []


def test_beamline_s1_simple():
    _test_notebook(
            os.path.join('samples', 'beamlines', 'S1_SPB_CRL_simplified', 'S1_SPB_CRL_simplified.ipynb'))

def test_tutoral_intro():
    _test_notebook(
        os.path.join('samples', 'Tutorials', 'Tutorial_intro.ipynb'))


def test_tutoral_1():
    _test_notebook(
        os.path.join('samples', 'Tutorials', 'Tutorial_case_1.ipynb'))


def test_tutoral_2():
    _test_notebook(
        os.path.join('samples', 'Tutorials', 'Tutorial_case_2.ipynb'))


def test_tutoral_2_new():
    _test_notebook(
        os.path.join('samples', 'Tutorials', 'Tutorial_case_2_new.ipynb'))


def test_tutoral_3():
    _test_notebook(
        os.path.join('samples', 'Tutorials', 'Tutorial_case_3.ipynb'))


def test_tutoral_3_new():
    _test_notebook(
        os.path.join('samples', 'Tutorials', 'Tutorial_case_3_new.ipynb'))