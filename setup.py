from os.path import join
import numpy as np
import sys

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.system_info import get_info
    
    # Get the numpy direction
    nppath=np.__path__[0]
    npinclude = join(nppath, 'core/include/numpy')

    # Create a config object
    config = Configuration('okada4py', parent_package, top_path)
    
    # Where are the librairies
    library_dir = []
    
    # Where are the includes
    include_dirs = ['src/', npinclude]
    
    # Which are the sources
    sources = ['src/dc3d.cpp', 'src/disloc3d.cpp']

    majorVersion = sys.version_info[0]
    if majorVersion == 2:
        sources.append('src/okada92.cpp')
    elif majorVersion == 3:
        sources.append('src/okada92_py3.cpp')
    else:
        print ('Unknown version of Python: Version {}'.format(majorVersion))
        sys.exit(1)

    # Additional flags
    CFLAGS = []

    # Create an extension
    config.add_extension('_okada92', 
            sources=sources, 
            include_dirs=include_dirs,
            extra_compile_args=CFLAGS)

    # All done
    return config

# Main
if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration, version='12.0.2')


