# Smina-plugin

<img width="346" alt="image" src="https://user-images.githubusercontent.com/102952395/162983719-287b957d-8cc1-4d77-9ef2-ba20833cfbb3.png"><img width="346" alt="image" src="https://user-images.githubusercontent.com/102952395/162983965-896ca987-7ba6-4b50-a326-d8f0aed61af0.png">


This is a plugin for PyMol 2.x or 3.x to use the vina/Autodock fork SMINA (https://sourceforge.net/projects/smina/) under Windows 10 or 11 with a linux subsystem installed. It is coded in Python 3 with PyQT5. It allows to dock single ligands or lists, with rigid or flexible sidechains and to post-refine the results. Import of residue selections or the display of results are intreractive beween the plugin and PyMol.   
# Why using smina under wsl2 in a Windows environnement ?
I like to work in a windows environnemt. Unfortunately there is no smina executable compiled for windows available. Therefore I wrote this plugin for PyMol under windows which executes the compiled smina static for linux available on the smina home page in a wsl2 environnement. The user has the impression to work completely under windows. I will also release a modified version of this interface for people who like to work completly under linux in the future. 
# How to install ?
Install the wsl2 environnement according to the well documented procedure : https://docs.microsoft.com/fr-fr/windows/wsl/install

You also need to install Openbabel (https://github.com/openbabel/openbabel/releases/tag/openbabel-3-1-1) under windows to let the plugin create pdbqt files.

Hints:

During wsL installation you have to execute the command line interface in the admin mode, else the error "the requested operation requires elevation" appears.
You can directly choose which linux distribution to install. I installed Ubuntu :> wsl -- install -d Ubuntu-18.04

Don't forget to change the subsystem from wsl to wsl2 : >wsl --set-default-version 2

If you use an open source version of PyMol from Christophe Gohlkes site (https://github.com/cgohlke/pymol-open-source-wheels) install the latest (!) Microsoft Visual C++ Redistributable packages for Visual Studio 2015, 2017, 2019, and 2022 files to make it work. You also have to install Python (versions 3.8 to 3.12 were tested) and either PySide2 (https://pypi.org/project/PySide2/) for Python versions up to 3.10 or PyQt5 for Python versions > 3.10 (using the pip installer : 'pip install -U PyQt5'). I installed the shiboken2-5.15.2 and PySide2-5.15.2 .whl files for my Python version 3.10 with pip (if the pip path is not in your PATH variable, open a command prompt in the pip.exe containing directory of your Python installation and enter : pip install PySide2) or PyQt5 directly with pip for Python versions > 3.10 (3.11 and 3.12 worked fine for me).
