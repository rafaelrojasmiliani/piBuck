from setuptools import setup, find_packages

setup(
    name='pibuck',
    version='0.1.0',
    packages=find_packages(),
    install_requires=['sympy', 'numpy', 'scipy'],
    python_requires='>=3.8',
)
