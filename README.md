# okada4py
Okada implementation in Python

## Citation
This is a python implementation of the solution proposed by Okada in 1992. Please cite:

Okada, Y. (1992), Internal deformation due to shear and tensile faults in a half-space, Bulletin of the Seismological Society of America, 82(2), 1018–1040.

## Install okada4py:

This installation requires meson. Install meson and meson-python first (check it out with pip).

### Install on a local directory

```
meson setup builddir --prefix /My/complete/path/to/the/install/dir
meson compile -C builddir
meson install -C builddir
```

Then update your PYTHONPATH variable to have it visible for python.
In your .bashrc, it would look like :
export PYTHONPATH=/My/complete/path/to/the/install/dir:$PYTHONPATH

and I am too lazy to find out how to install this in the root of python (it does not work yet for me)


