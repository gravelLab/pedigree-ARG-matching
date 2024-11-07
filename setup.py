from setuptools import setup, find_packages


def parse_requirements(filename):
    with open(filename, 'r') as f:
        return f.read().splitlines()


# Setup configuration
setup(
    name='pedigree-arg-matching',
    version='0.1.0',
    author="Andrii Serdiuk, Simon Gravel",
    author_email="andrii_serdiuk@outlook.com, simon.gravel@mcgill.ca",
    packages=find_packages(),  # Automatically find all packages and subpackages
    license='MIT',
    python_requires=">=3.9",
    install_requires=parse_requirements('requirements.txt'),
)
