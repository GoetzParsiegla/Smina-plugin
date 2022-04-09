# Smina-plugin
This is a plugin for PyMol 2.x tu use the vina/Autodock fork SMINA (https://sourceforge.net/projects/smina/) under Windows 10 or 11 with a linux subsystem installed. It is coded in Python 3 with PyQT5 (PySide2). It allows to dock single ligands or lists, with rigid or flexible sidechains and to post-refine the results. Residues selections from PyMol can be  imported and results are displayed in PyMol.   
# Why using smina under wsl2 in a Windows environnement ?
I like to work in a windows environnemt. Unfortunately there is no smina executable compiled for windows available. Therfore I wrote this plugin which works in PyMol under windows and executes the static smina version for linux available on the smina home page in a wsl2 environnement. The user has the impression to work completely under windows. I will also release a modified version of this interface for people who like to work completly under linux in the future. 
# How to install ?
Install the wsl2 environnement according to the well documented procedure : https://docs.microsoft.com/fr-fr/windows/wsl/install

Hints:

You have to execute the command line interface in the admin mode, else the error "the requested operation requires elevation" appears.
You can directly choose which linux distribution to install. I took :> wsl -- install -d Ubuntu-18.04

D'ont forget to change the subsystem from wsl to wsl2 : >wsl --set-default-version 2

If you use an open source version of PyMol from Christophe Gohlkes site, you have to installe Python (3.8 and 3.9 will work fine) and PySide2 (https://pypi.org/project/PySide2/) under Windows according to your Python version. I installed the shiboken2-5.12.2 and PySide2-5.12.2 .whl files for my Python version manually with pip. 
