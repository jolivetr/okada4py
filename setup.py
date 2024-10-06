from setuptools import setup
from setuptools.build_meta import build_sdist, build_wheel

setup(
    name='okada4py',
    version='1.0.1',
    setup_requires=['meson-python'],
    build_sdist=build_sdist,
    build_wheel=build_wheel,
)