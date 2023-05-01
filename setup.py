from setuptools import setup 
import sys 
import os 
import shutil
import distutils.cmd


VERSION = "0.0.1"

class PypiCommand(distutils.cmd.Command):
    
    description = "Build and upload for PyPI."
    user_options = []
    
    def initialize_options(self):
        pass
    
    
    def finalize_options(self):
        pass
    
    
    def run(self):
        try:
            shutil.rmtree("dist/")
        except FileNotFoundError:
            pass
        
        wheel_file = "sde-{}-py3-none-any.whl".format(VERSION)
        tar_file = "sde-{}.tar.gz".format(VERSION)
        
        os.system("{} setup.py sdist bdist_wheel".format(sys.executable))
        os.system("twine upload dist/{} dist/{}".format(wheel_file, tar_file))


setup(
    name="sde",
    version=VERSION,
    description="A simple tool for numerical stochastic differential equations",
    author="Yuqiu Yang",
    author_email="yuqiuy@smu.edu",
    long_description_content_type="text/markdown",
    long_description=open("README.md").read(),
    packages=["sde"],
    python_requires=">=3.9",
    install_requires=[
        'matplotlib==3.7.1',
        'numpy==1.24.2',
        'tqdm==4.65.0',
        'pandas==2.0.0',
        'seaborn==0.12.2'
    ],
    test_requires=[
        'pytest==7.1.2',
        'coverage==6.3.2',
        'pytest-cov==3.0.0',
        'pytest-mock==3.10.0'
    ],
    classifiers=[
        "Programming Language :: Python :: 3 :: Only",
        "Natural Language :: English",
        "Intended Audience :: Science/Research",
    ],
    cmdclass={
        "pypi": PypiCommand
    }
)