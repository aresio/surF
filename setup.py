import setuptools


setuptools.setup(
    name="aresio", # Replace with your own username
    version="0.0.1",
    author="Marco S. Nobile",
    author_email="m.s.nobile@tue.nl",
    description="surF - Fourier surrogate modeling",
    long_description=open('README.md', encoding='utf-8').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/surfer",
    packages=setuptools.find_packages(),
     license='LICENSE.txt',
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)