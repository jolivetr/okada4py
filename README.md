# okada4py
Okada implementation in Python

## Citation
This is a python implementation of the solution proposed by Okada in 1992. Please cite:

Okada, Y. (1992), Internal deformation due to shear and tensile faults in a half-space, Bulletin of the Seismological Society of America, 82(2), 1018–1040.

and this implementation: [![DOI](https://zenodo.org/badge/212802180.svg)](https://doi.org/10.5281/zenodo.14170826) 

[Dr Leah Langer](https://en-exact-sciences.tau.ac.il/profile/llanger) kindly implemented her receiver elevation correction to provide a first order correction of the effect of topography on the evaluation of Green's functions. If you use this, please cite (and read) this paper: [doi: 10.1016/j.tecto.2020.228566](https://doi.org/10.1016/j.tecto.2020.228566)

## Install okada4py:

### Install with pip

Install directly from the repository root:

```bash
python -m pip install .
```

This builds the extension with `meson-python` and installs `okada4py` into the active Python environment.

### Build requirements

`pip` will install the Python build requirements declared in `pyproject.toml`. You only need a working C++ toolchain available on your system.

### Publish on PyPI

The repository now includes a GitHub Actions workflow that builds the source distribution and wheel on Linux and macOS for multiple Python versions.

For PyPI trusted publishing, register this GitHub Actions publisher on PyPI:

1. Owner: `jolivetr`
2. Repository: `okada4py`
3. Workflow: `.github/workflows/python-package.yml`
4. Environment: `pypi`

To publish a release on PyPI:

1. Create a GitHub release from a version tag.
2. Configure PyPI trusted publishing for this repository using the values above.
3. Publish the GitHub release to trigger the upload job.

Before tagging a release locally, you can verify the packaging with:

```bash
python -m build
python -m pip install --force-reinstall dist/*.whl
```

### Install on a local directory

```
meson setup builddir --prefix /My/complete/path/to/the/install/dir
meson compile -C builddir
meson install -C builddir
```

Then update your PYTHONPATH variable to have it visible for python.
In your .bashrc, it would look like :
export PYTHONPATH=/My/complete/path/to/the/install/dir:$PYTHONPATH

### Install to the python root directory

```
meson setup builddir --prefix=$(python3 -c "import site; print(site.getsitepackages()[0])")
meson compile -C builddir
meson install -C builddir
```


