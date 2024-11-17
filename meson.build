project('okada4py', 'cpp', version : '1.0.1', default_options : ['cpp_std=c++11'])

py3 = import('python').find_installation('python')

numpy_include = run_command(py3, '-c', 'import numpy; print(numpy.get_include())').stdout().strip()

src = ['src/dc3d.cpp', 'src/disloc3d.cpp', 'src/okada92_py3.cpp']

okada92 = py3.extension_module('_okada92',
  sources : src,
  include_directories : include_directories('src', numpy_include),
  install : true,
  install_dir : 'okada4py'
)

install_data('src/__init__.py', install_dir : 'okada4py')
