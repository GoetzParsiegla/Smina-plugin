# Smina-plugin
This is a plugin for PyMol 2.x tu use the vina/Autodock fork SMINA (https://sourceforge.net/projects/smina/) under Windows 10 with a linux subsystem installed. It is coded in Python 3 with PyQT5 (PySide2). It allows to dock single ligands or lists, with rigid or flexible sidechains and to post-refine the results. Residues selections from PyMol can be  imported and results are displayed in PyMol.   
# Why using smina under wsl2 in a Windows environnement ?
I like to work in a windows environnemt. Unfortunately there is no smina executable compiled for windows available. Therfore I wrote this plugin which works in PyMol under windows and executes the static smina version for linux available on the smina home page in a wsl2 environnement. The user has the impression to work completely under windows. I will also release a modified version of this interface for people who like to work completly under linux in the future. 
# How to install ?
Install the wsl2 environnement as well documented in the microsaft help page : https://docs.microsoft.com/fr-fr/windows/wsl/install

Remarques:
Il faut executer l'invité de commandes en mode administrateur, sinon l'erreur "'opération demandée nécessite une élévation" s'affiche.
On peut directement choisir quel installation linux. J'ai choisi :> wsl -- install -d Ubuntu-18.04
Il ne faut pas oublier de changer le subsystem de wsl à wsl2 : >wsl --set-default-version 2
