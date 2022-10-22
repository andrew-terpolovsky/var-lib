import setuptools

setuptools.setup(
    version="1.1.0",
    name="variants_lib",
    packages=setuptools.find_packages(),
    long_description="",
    long_description_content_type="text",
    license="",
    install_requires=["boto3", "pyarrow==9.0.0"],
    python_requires=">=3.8",
)
