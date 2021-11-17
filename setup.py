from setuptools import setup

setup(
    name="cctk",
    author="Alan Collins",
    license="GPL-3.0",
    install_requires=[
        "pip",
        "dendropy",
        "matplotlib",
        "setuptools",
    ],
    scripts=[
        'CRISPR_comparison_toolkit/CRISPRtree.py',
        "CRISPR_comparison_toolkit/CRISPRcasfinder2arrays.py",
        "CRISPR_comparison_toolkit/CRISPRdiff.py",
        "CRISPR_comparison_toolkit/CRISPRspacers2network.py",
        "CRISPR_comparison_toolkit/CRISPRtree.py",
        "CRISPR_comparison_toolkit/CRtree_constrained.py",
        "CRISPR_comparison_toolkit/minced2arrays.py",
        "CRISPR_comparison_toolkit/narbl2arrays.py",
        "CRISPR_comparison_toolkit/process_CRISPR_blast.py",
        "CRISPR_comparison_toolkit/Process_CRISPR_cas_finder_out.py",
        "CRISPR_comparison_toolkit/reps2arrays.py",
    ],
    packages = ["CRISPR_comparison_toolkit/cctk"]
)

