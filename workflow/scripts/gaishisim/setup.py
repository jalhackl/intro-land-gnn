import os.path
from setuptools import setup, find_packages

setup(
    name="gaishisim",
    version="1.0.0",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "demes",
        "msprime",
        "numpy",
        "pandas",
        "scikit-allel",
        "scikit-learn",
        "scipy",
    ],
    entry_points={"console_scripts": ["gaishisim=gaishisim.__main__:main"]},
)
