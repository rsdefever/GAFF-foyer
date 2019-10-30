from setuptools import setup

setup(
    name="gafffoyer",
    version="0.0",
    description="GAFF plugin for foyer",
    url="https://github.com/rsdefever/GAFF-foyer",
    author="Ryan S. DeFever",
    author_email="rdefever@nd.edu",
    license="MIT",
    install_requires="foyer",
    entry_points={
        'foyer.forcefields': [
            "GAFF = gafffoyer.gafffoyer:GAFF"]
    },
    packages=["gafffoyer"],
    zip_safe=False,
)
