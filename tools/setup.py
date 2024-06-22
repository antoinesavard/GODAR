# setup.py

from setuptools import setup, find_packages

setup(
    name='tools',
    version='0.1',
    description='A package that is used in combination with godar to create inputs and plot/analyse the outputs.'
    author='Antoine Savard'
    author_email='antoine.savard@mail.mcgill.ca'
    packages=find_packages(),
    install_requires=[
        # Add any dependencies here
    ],
    entry_points={
        'console_scripts': [
            # Define any command line scripts here
        ],
    },
)