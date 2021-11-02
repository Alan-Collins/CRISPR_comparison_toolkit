from setuptools import setup

setup(
    name="cctk",
    install_requires=[
        "pip",
        "dendropy",
        "matplotlib",
        "setuptools",
    ],
    scripts=['CRISPR_mp_by_topology.py',
        'CRISPR_nj.py',
        'CRISPR_mp.py',
        'CRtree_constrained.py',
        'run_CRISPR_mp_all.py']
)

