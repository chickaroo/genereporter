from setuptools import setup, find_packages

setup(
    name='genereporter',
    version='0.1',
    author='Samantha Bening',
    author_email='samantha.bening@helmholtz-munich.de',
    description='focus on only one jean guys',
    url="https://github.com/chickaroo/genereporter",
    license='MIT',
    packages=['genereporter'],
    #package_data={'speos': ['LICENSE.md', "speos/adjacencies.json", "speos/mapping.json", "speos/utils/config_default.yaml"],
                  #"extensions": ["extensions/mapping.json", "extensions/adjacencies.json", "extensions/datasets.json"]},
    #package_dir={"": "coregenes"},
    install_requires=['anndata', 'anndata2ri', 'biopython', 'ca-certificates', 'celltypist', 'certifi', 'ipykernel', 'jupyter', 'loompy', 'matplotlib', 'multicore-tsne', 'numpy', 'openssl', 'pandas', 'plotly', 'pytables', 'pyvis', 'scanpy', 'scipy', 'scvi-tools', 'seaborn', 'xmltodict']



)