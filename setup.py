"""
QM_Project
Refactoring of the QM project from MolSSI Software Summer School 2019
"""
import setuptools
import os, sys
from os import path
import qm_project.semi_empirical_model

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md')) as f:
    long_description = f.read()

#################################################################
# Build our C++ module
# NOTE: Pybind11/Eigen were installed into CONDA_PREFIX
#       so we need to add that to the include paths
conda_prefix = os.environ['CONDA_PREFIX']
eigen_path = os.path.join(conda_prefix, 'include', 'eigen3')

cpp_module = setuptools.Extension('qm_project.semi_empirical_model',
                        include_dirs = [eigen_path],
                        extra_compile_args = ['-std=c++11'],
                        sources = ['qm_project/fock_matrix.cpp',
                                   'qm_project/export_qm2.cpp'])
#################################################################

if __name__ == "__main__":
    setuptools.setup(
        name='semi_empirical_model',
#        version=qm_project.__version__,
        author='Gaurav Vishwakarma',
        author_email='gvishwak@buffalo.edu',
        project_urls={
            'Source': 'https://github.com/MolSSI-Education/qm_2019_sss_2',
        },
        description=
        'Semi-Empirical Model for Noble Gases',
        long_description=long_description,
        # scripts=['lib/'],
        license='BSD-3C',
        packages=setuptools.find_packages(),

        install_requires=[
            'numpy'
        ],
        extras_require={
            'docs': [
                'sphinx',
                'sphinxcontrib-napoleon',
                'sphinx_rtd_theme',
                'numpydoc',
            ],
            'tests': [
                'pytest',
                'pytest-cov',
                'pytest-pep8',
                'tox',
            ],
        },
        tests_require=[
            'pytest',
            'pytest-cov',
            'pytest-pep8',
            'tox',
        ],
        classifiers=[
            'Natural Language :: English',
            'Intended Audience :: Science/Research',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
        ],
        
        # Include the compiled extension
        ext_modules = [cpp_module],

        zip_safe=False,
    )
    
