import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="bacphy",
    version="0.1",
    author="Chuhao-Li",
    author_email="a879942613@qq.com",
    description="bacterial phylogenetic analysis tools",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="git@github.com:Chuhao-Li/bacphy.git",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)

