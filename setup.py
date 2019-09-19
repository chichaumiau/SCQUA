from setuptools import setup, find_packages

def readme():
    with open('README.md') as f:
        return f.read()


print(find_packages())

setup(
        name='SCQUA',
        version='0.0.2',
        description='Single-Cell sequencing QUAlity assessment',
        long_description=readme(),
        packages=find_packages(),
        install_requires=['numpy', 'pandas', 'scanpy'],
        author='Chichau Miau',
        author_email='zmiao@ebi.ac.uk',
        license='MIT'
    )
