import os
from pathlib import Path
from setuptools import setup, find_packages
from setuptools.command.install import install
import subprocess


def compile():
	src_dir = os.path.join(os.path.dirname(__file__), "pybiqcrunch", "biqcrunch", "BiqCrunch_second_release", "BiqCrunch", "src")
	subprocess.check_call("make", cwd=src_dir, shell=True)

class PreInstall(install):
    def run(self):
        compile()
        super().run()

setup(
	name="pybiqcrunch",
	description="A Python wrapper for the BiqCrunch quadratic solver",
    version="1.0",
    packages=find_packages(),
    install_requires=["numpy"],
    package_data={
    	"pybiqcrunch": [str(Path(*path.parts[1:])) for path in Path("pybiqcrunch/biqcrunch").rglob("*")]
    },
    cmdclass={"install": PreInstall}
)
