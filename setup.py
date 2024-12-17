from setuptools import setup

setup(
    name="cctk",
    python_requires='>=3.6',
    version="1.0.3",
    author="Alan Collins",
    author_email="crisprtoolkit@gmail.com",
    url="https://github.com/Alan-Collins/CRISPR_comparison_toolkit",
    license="GPL-3.0",
    install_requires=[
        "pip",
        "dendropy",
        "matplotlib",
        "setuptools",
    ],
    scripts=[
        "CRISPR_comparison_toolkit/cctk",
    ],
    packages = ["cctkpkg"],
    package_dir={"": "CRISPR_comparison_toolkit"}
)

