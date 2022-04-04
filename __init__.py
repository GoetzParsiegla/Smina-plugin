# This Python 3.x file uses the following encoding: utf-8
# Smina-plugin for PyMol 3.x (Windows version avec WSL2) Copyright Notice.
# ====================================================================
#
# The smina-plugin source code is copyrighted, but you can freely
# use and copy it as long as you don't change or remove any of the
# copyright notices.
#
# ----------------------------------------------------------------------
# Smina plugin is Copyright (C) 2020 by Goetz Parsiegla
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name of Goetz Parsiegla not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# GOETZ PARSIEGLA DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
# SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS.  IN NO EVENT SHALL GOETZ PARSIEGLA BE LIABLE FOR ANY
# SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
# RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
# CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------
#
# Executes smina for molecular docking.
# If you find bugs or have any suggestions for future versions contact me:
# goetz.parsiegla@imm.cnrs.fr

import sys
import os
import time

# pymol.Qt is a wrapper which provides the PySide2/Qt5/Qt4 interface
# if it has been installed in python before !
from pymol.Qt import QtWidgets, QtCore
from pymol import cmd
from glob import glob
from pymol.cgo import *
from pymol.vfont import plain

__version__ = "0.9.8.0"

def __init_plugin__(app=None):
    '''
    Add an entry to the PyMOL "Plugin" menu
    '''
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('Smina', run_plugin_gui)


# global reference to avoid garbage collection of our dialog
dialog = None

def run_plugin_gui():
    global dialog
    if dialog is None:
        dialog = QtWidgets.QDialog() # now global dialog holds a Qtobject

    # filename of our UI file
    uifile = os.path.join(os.path.dirname(__file__), 'form.ui')

    # load the UI file into our dialog
    from pymol.Qt.utils import loadUi
    form = loadUi(uifile, dialog)
    Smina(form) # call the plugin class and pass the form as an argument
    dialog.show()

# --------------- Plugin code starts here --------------------

class Smina:
    def __init__(self, form):
        self.form = form     

        # get a temporary file directory
        if not sys.platform.startswith('win'):
            home = os.environ.get('HOME')
        else:
            home = os.environ.get('HOMEPATH')

        tmp_dir = os.path.join(home,'.PyMol_plugin')
        if not os.path.isdir(tmp_dir):
            os.mkdir(tmp_dir)
            print("Created temporary files directory:  %s" % tmp_dir)

        # Assign variables to lineEdits and spinBoxes  
        statusline = self.form.lineEdit
        self.smina_location = self.form.lineEdit_2
        self.openbabel_location = self.form.lineEdit_3
        self.ligand_dir_location = self.form.lineEdit_4
        self.scoring_table_dir_location = self.form.lineEdit_8
        self.center_selection = self.form.lineEdit_5
        self.form.doubleSpinBox.setMinimum(-1000)
        self.form.doubleSpinBox.setMaximum(1000)
        self.form.doubleSpinBox_2.setMinimum(-1000)
        self.form.doubleSpinBox_2.setMaximum(1000)
        self.form.doubleSpinBox_3.setMinimum(-1000)
        self.form.doubleSpinBox_3.setMaximum(1000)
        self.form.doubleSpinBox_4.setMaximum(14)
        self.form.doubleSpinBox_4.setSingleStep(0.1)
        self.form.doubleSpinBox_5.setMaximum(0)
        self.form.doubleSpinBox_5.setMinimum(-1000000)
        self.form.doubleSpinBox_5.setSingleStep(1)

        # make Buttongroups
        self.Buttongroup_1 = QtWidgets.QButtonGroup()
        self.Buttongroup_1.addButton(self.form.radioButton)
        self.Buttongroup_1.addButton(self.form.radioButton_2)
        self.Buttongroup_2 = QtWidgets.QButtonGroup()
        self.Buttongroup_2.addButton(self.form.radioButton_3)
        self.Buttongroup_2.addButton(self.form.radioButton_4)
        
        # defaults
        self.config_settings = {}
        self.smina_config = {}
        self.ligand_dir_path = ""
        self.scoring_table_dir_path = ""
        self.current_flexibles = []
        self.form.doubleSpinBox_4.setValue(7.00)
        self.current_ligands = []
        self.loaded_poses_list = []
        self.firstrun = True
        self.form.spinBox.setValue(8) # exhaustiveness
        self.form.spinBox_2.setValue(9) # maxposes
        self.form.spinBox_3.setValue(4) # autobox_buf
        self.form.doubleSpinBox_5.setValue(0) # seed
        self.scoring_table_file = ""
        vinascr = [[-0.035579,'gauss(o=0,_w=0.5,_c=8)'],[-0.005156,'gauss(o=3,_w=2,_c=8)'],\
                    [0.840245,'repulsion(o=0,_c=8)'],[-0.035069,'hydrophobic(g=0.5,_b=1.5,_c=8)'],\
                    [-0.587439,'non_dir_h_bond(g=-0.7,_b=0,_c=8)'],[1.923000,'num_tors_div']]
        approx_methods_list = ['linear', 'spline', 'exact']        
        self.scoring_list = vinascr
        self.current_poses_list = []
        #-----------------------------------------------------------

        # Config page

        install_text = """
        <html><body>This is the PyMol 2.x plugin for smina. It needs at least
        Windows 10 build 2004 (The may 2020 update) with the wsl2 (Windows
        subsystem for linux version 2) component activated.<br> An Ubuntu >14.04
        installation is required to work with smina (18.04. confirmed).<br>You find
        instructions to install wsl2 at: <a href="https://docs.microsoft.com/fr-fr/windows/wsl/install-win10">
        docs.microsoft.com/fr-fr/windows/wsl/install-win10</a>.<br>
        The Openbabel executive (obabel.exe) as well as the smina.static
        executive for linux must be located in a Windows directory. You will find
        the latest stable version of Openbabel at :<br> <a href="https://github.com/openbabel/openbabel/releases">
        github.com/openbabel/openbabel/releases</a>
        and smina.static at : <a href="https://sourceforge.net/projects/smina/files/">
        sourceforge.net/projects/smina/files/</a>.
        </body></html>
        """

        def set_statusline(text):
            statusline.clear()
            statusline.insert(text)

        # Config functions

        def get_smina_location():
            filedialog = QtWidgets.QFileDialog()
            filename = filedialog.getOpenFileName(None, "smina.static location", os.getcwd(), 'linux.static files (*.static)')
            filename = filename[0]
            set_smina_location(filename)

        def get_openbabel_location():
            filedialog = QtWidgets.QFileDialog()
            filename = filedialog.getOpenFileName(None, "obabel.exe location", 'C:\\Program Files', 'exe files (*.exe)')
            filename = filename[0]
            set_openbabel_location(filename)

        def get_ligand_dir_path():
            dirdialog = QtWidgets.QFileDialog()
            dirname = dirdialog.getExistingDirectory(None, "Ligand directory", os.getcwd())+'/'
            set_ligand_dir_path(dirname)

        def get_scoring_table_dir_path():
            dirdialog = QtWidgets.QFileDialog()
            dirname = dirdialog.getExistingDirectory(None, "Scoring-table directory", os.getcwd())+'/'
            set_scoring_table_dir_path(dirname)

        def set_smina_location(filename):
            self.smina_location.clear()
            self.smina_location.insert(filename)
            self.smina_exe = filename
            self.config_settings['smina_exe'] = filename

        def set_openbabel_location(filename):
            self.openbabel_location.clear()
            self.openbabel_location.insert(filename)
            self.openbabel_exe = filename
            self.config_settings['openbabel_exe'] = filename

        def set_ligand_dir_path(dirname):
            self.ligand_dir_location.clear()
            self.ligand_dir_location.insert(dirname)            
            self.ligand_dir_path = dirname
            self.config_settings['ligand_dir_path'] = dirname
            
        def set_scoring_table_dir_path(dirname):
            self.scoring_table_dir_location.clear()
            self.scoring_table_dir_location.insert(dirname)
            self.scoring_table_dir_path = dirname
            self.config_settings['scoring_table_dir_path'] = dirname

        def read_plugin_config_file():
            config_file_name = os.path.join(tmp_dir,"smina_plugin.conf")
            self.config_settings = {}
            self.config_settings['smina_exe'] = ''
            self.config_settings['openbabel_exe'] = ''
            self.config_settings['ligand_dir_path'] = ''
            self.config_settings['scoring_table_dir_path'] = ''
            if os.path.isfile(config_file_name):
                set_statusline('Reading configuration file: %s' % config_file_name)
                with open(config_file_name,'r') as f:
                    lst = f.readlines()
                    for line in lst:
                        if line[0]!='#':
                            entr = line.split('=')
                            self.config_settings[entr[0].strip()] = entr[1].strip()
                    set_smina_location(self.config_settings['smina_exe'])
                    set_openbabel_location(self.config_settings['openbabel_exe'])
                    set_ligand_dir_path(self.config_settings['ligand_dir_path'])
                    set_scoring_table_dir_path(self.config_settings['scoring_table_dir_path'])
            else:
                set_statusline('ERROR : Plugin configuration file not found')
            return self.config_settings

        def save_plugin_config_file():
            config_file_name = os.path.join(tmp_dir,"smina_plugin.conf")
            with open(config_file_name,'w') as fp:
                print('#========================================', file=fp)
                print('# Smina Plugin configuration file', file=fp)
                self.config_settings['smina_exe'] = self.smina_location.text()
                self.config_settings['openbabel_exe'] = self.openbabel_location.text()
                self.config_settings['ligand_dir_path'] = self.ligand_dir_location.text()
                self.config_settings['scoring_table_dir_path'] = self.scoring_table_dir_location.text()
                for key, val in self.config_settings.items():
                    print(key, '=', val, file=fp)
            set_statusline('Wrote smina-plugin configuration file %s' % config_file_name)

        # Run page buildup
        self.form.textBrowser.setHtml(install_text)
        self.form.textBrowser.setOpenExternalLinks(True)
        self.config_settings = read_plugin_config_file()
        
        # Box page

        Box_text = """
        <html><body>This is the page for the creation of a Box parameter file ('receptor'_config.txt),
        containing the coordinates of the Box defining the search space for
        docking.<br> * By default a standard box is centered around the COM of the selected receptor.
        The box can be centered around a selection of a residue imported from Pymol or entered manually
        in the coordinate input spaces.<br> * Existing parameter files for selected receptor objects are
        loaded automatically.<br> * An existing parameter file may also be loaded from eleswhere using
        'Load config file'.<br>Don't forget to write the current parameter file in the current directory
        with 'Write config file' beforre procceding. </body></html>
        """

        # Box functions        

        def read_config_file(filename):
            Box_settings = {}
            with open(filename, 'r') as f:
                lst = f.readlines()
                for line in lst:
                    entr = line.split('=')
                    Box_settings[entr[0].strip()] = entr[1].strip()
                self.form.doubleSpinBox.setValue(float(Box_settings['center_x']))
                self.form.doubleSpinBox_2.setValue(float(Box_settings['center_y']))
                self.form.doubleSpinBox_3.setValue(float(Box_settings['center_z']))
                self.form.spinBox_4.setValue(int(Box_settings['size_x']))
                self.form.spinBox_5.setValue(int(Box_settings['size_y']))
                self.form.spinBox_6.setValue(int(Box_settings['size_z']))

        def load_smina_config_file():
            prot = self.form.comboBox_3.currentText()
            filename = prot+"_config.txt"
            read_config_file(filename)

        def load_custom_config_file():
            filedialog = QtWidgets.QFileDialog()
            filename = filedialog.getOpenFileName(None, "load custom_config file", os.getcwd(), 'config files (*.txt)')
            filename = filename[0]
            read_config_file(filename)

        def get_Boxcenter_selection():
            prot = self.form.comboBox_3.currentText()
            if cmd.get_names("selections"):
                cmd.create("box_center", "sele")
                myspace ={'lst': []}
                cmd.iterate("box_center", "lst.append(resn+resi)", space=myspace)
                sel=[(myspace["lst"][0])]
                for i in range(len(myspace["lst"])-1):
                    if (myspace["lst"][i]) == (myspace["lst"][i+1]):
                        continue
                    else:
                        sel.append(myspace["lst"][i+1])
                self.form.lineEdit_5.clear()
                self.form.lineEdit_5.insert(",".join(sel))
                cmd.delete("box_center")
                Boxcenter_on_selection()
                
        def Boxcenter_on_selection():            
            res_center = self.form.lineEdit_5.text()
            res_i = res_center.strip('ACEGHILMNOPRSTUVY')
            cmd.select("sele",("resi "+res_i))
            cmd.create("box_center","sele")
            cmd.delete("sele")
            myspace ={'lst': []}
            cmd.iterate_state(0,"box_center", "lst.append((x,y,z))", space=myspace)
            self.form.doubleSpinBox.setValue(float(myspace["lst"][1][0]))
            self.form.doubleSpinBox_2.setValue(float(myspace["lst"][1][1]))
            self.form.doubleSpinBox_3.setValue(float(myspace["lst"][1][2]))
            cmd.delete("box_center")

        def set_Receptor():
            prot = self.form.comboBox_3.currentText()
            if prot != "":
                if os.path.isfile(prot+"_config.txt"):
                    set_statusline('Reading configuration file: %s' % (prot+"_config.txt"))
                    load_smina_config_file()
                else:
                    com = cmd.centerofmass(prot)
                    self.form.doubleSpinBox.setValue(com[0])
                    self.form.doubleSpinBox_2.setValue(com[1])
                    self.form.doubleSpinBox_3.setValue(com[2])
            else:
                None
            
        def set_Boxsize_defaults():
            self.form.spinBox_4.setValue(20)
            self.form.spinBox_5.setValue(20)
            self.form.spinBox_6.setValue(20)


# ---- The following code of the box calculation and th CGO object is extracted
# ---- from the Autodock/Vina plugin copyrighted by Daniel Seeliger :
#
# Autodock/Vina plugin  Copyright Notice
# ============================
#
# The Autodock/Vina plugin source code is copyrighted, but you can freely use and
# copy it as long as you don't change or remove any of the copyright
# notices.
#
# ----------------------------------------------------------------------
# Autodock/Vina plugin is Copyright (C) 2009 by Daniel Seeliger
#
#                        All Rights Reserved
#
# Permission to use, copy, modify, distribute, and distribute modified
# versions of this software and its documentation for any purpose and
# without fee is hereby granted, provided that the above copyright
# notice appear in all copies and that both the copyright notice and
# this permission notice appear in supporting documentation, and that
# the name of Daniel Seeliger not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# DANIEL SEELIGER DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS
# SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS.  IN NO EVENT SHALL DANIEL SEELIGER BE LIABLE FOR ANY
# SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER
# RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF
# CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# ----------------------------------------------------------------------

        def calculate_Box():
            x = float(self.form.doubleSpinBox.value())
            y = float(self.form.doubleSpinBox_2.value())
            z = float(self.form.doubleSpinBox_3.value())
            xpts = int(self.form.spinBox_4.value())
            ypts = int(self.form.spinBox_5.value())
            zpts = int(self.form.spinBox_6.value())
            size = [xpts, ypts, zpts]
            xmax = x + size[0]/2.
            xmin = x - size[0]/2.
            ymax = y + size[1]/2.
            ymin = y - size[1]/2.
            zmax = z + size[2]/2.
            zmin = z - size[2]/2.
            box_edge_x = [xmin,xmax]
            box_edge_y = [ymin,ymax]
            box_edge_z = [zmin,zmax]
            box_coords  = [box_edge_x,box_edge_y,box_edge_z]
            cmd.delete('box')
            display_Box(box_coords)

        def display_Box(box):
            view = cmd.get_view()
            name = "box"
            obj = []
            # build cgo object
            color = [1.,1.,1.]
            cylinder_size = 0.2
            for i in range(2):
                for k in range (2):
                    for j in range(2):
                        if i != 1:
                            obj.append(CYLINDER)
                            obj.extend([box[0][i],box[1][j],box[2][k]])
                            obj.extend([box[0][i+1],box[1][j],box[2][k]])
                            obj.append(cylinder_size)
                            obj.extend(color)
                            obj.extend(color)
                            obj.append(COLOR)
                            obj.extend(color)
                            obj.append(SPHERE)
                            obj.extend([box[0][i],box[1][j],box[2][k],cylinder_size])
                            
                        if j != 1:
                            obj.append(CYLINDER)
                            obj.extend([box[0][i],box[1][j],box[2][k]])
                            obj.extend([box[0][i],box[1][j+1],box[2][k]])
                            obj.append(cylinder_size)
                            obj.extend(color)
                            obj.extend(color)
                            obj.append(COLOR)
                            obj.extend(color)
                            obj.append(SPHERE)
                            obj.extend([box[0][i],box[1][j+1],box[2][k],cylinder_size])
                        if k != 1:
                            obj.append(CYLINDER)
                            obj.extend([box[0][i],box[1][j],box[2][k]])
                            obj.extend([box[0][i],box[1][j],box[2][k+1]])
                            obj.append(cylinder_size)
                            obj.extend(color)
                            obj.extend(color)
                            obj.append(COLOR)
                            obj.extend(color)
                            obj.append(SPHERE)
                            obj.extend([box[0][i],box[1][j],box[2][k+1],cylinder_size])
            axes = [[2.0,0.0,0.0],[0.0,2.0,0.0],[0.0,0.0,2.0]]
            xpos = [box[0][1]+(box[0][1]-box[0][0])/5.,box[1][0],box[2][0]]
            cyl_text(obj,plain,xpos,'X',0.10,axes=axes)
            ypos = [box[0][0],box[1][1]+(box[1][1]-box[1][0])/5,box[2][0]]
            cyl_text(obj,plain,ypos,'Y',0.10,axes=axes)
            zpos = [box[0][0],box[1][0],box[2][1]+(box[2][1]-box[2][0])/5]
            cyl_text(obj,plain,zpos,'Z',0.10,axes=axes)
            cmd.load_cgo(obj,name)
            cmd.set_view(view)            

        def hide_Box():
            cmd.delete("box")          

# Autodock/Vina code ends here
# ---------------------------------------------------------------------------------

        def save_smina_config_file():
            prot = self.form.comboBox_3.currentText()
            if prot == "":
                set_statusline("No structure selected")
            else:
                smina_config_file_name = prot+"_config.txt" 
                with open(smina_config_file_name,'w') as fp:
                    self.smina_config['center_x'] = self.form.doubleSpinBox.value()
                    self.smina_config['center_y'] = self.form.doubleSpinBox_2.value()
                    self.smina_config['center_z'] = self.form.doubleSpinBox_3.value()
                    self.smina_config['size_x'] = self.form.spinBox_4.value()
                    self.smina_config['size_y'] = self.form.spinBox_5.value()
                    self.smina_config['size_z'] = self.form.spinBox_6.value()
                    for key, val in self.smina_config.items():
                        print(key, '=', val, file=fp)
                set_statusline('Wrote smina configuration file %s' % smina_config_file_name)

        # Run page buildup
        self.form.textBrowser_3.setHtml(Box_text)
        set_Boxsize_defaults()

        # Receptor page

        receptor_text = """
        <html><body> * Objects already loaded in PyMol are imported in the smina-plugin when opening it.<br>
         * You may import other receptors from Pymol by using 'Import object'.<br>
         * The selected receptor in the receptor Window will be used in all other plugin pages.<br>
         * If no 'receptor'.pdbqt file of the 'receptor' is present in the working directory, it will be automatically
        created with openbabel.<br> * If no flexible residues are entered, 'Rigid Sidechains' are selected for docking.<br>
         * To use 'Flexible sidechains' for docking, select flexible residues of the receptor in PyMol and import them
         in the flexible residues window. You can still delete residues in the window.<br> * Actually only
        'gasteiger' charges are supported in the .pdbqt file, but you might select a pH for the protonation state
        of the receptor and the ligand when automatically creating their '.pdbqt' files. </body></html>
        """

        # Receptor functions

        def import_objects():
            self.form.comboBox.clear()
            self.form.comboBox_3.clear()
            lst = cmd.get_names()
            if 'sele' in lst:
                lst.remove('sele')
            if 'cgo' in lst:
                lst.remove('cgo')
            object_list = lst
            self.form.comboBox.addItems(object_list)
            self.form.comboBox_3.addItems(object_list)

        def import_charge_models():
            charge_model_list= ["gasteiger"]
            self.form.comboBox_5.addItems(charge_model_list)

        def select_current_flexibles():
            flex_id = []
            for i in range(len(self.current_flexibles)):
                flex_id.append("resi "+self.current_flexibles[i].strip(':ABCDEFGHIJKLMNOPQRSTUVWXYZ'))
            sel = str(flex_id).strip('[]').translate(str.maketrans('','', "'"))
            cmd.select("flexres", (sel))
            cmd.disable("flexres")
            return sel

        def import_flexibles():
            prot = self.form.comboBox_3.currentText()
            if cmd.get_names("selections"):
                cmd.create("flexibles", "sele")
                cmd.delete("sele")
                myspace ={'lst': []}
                cmd.iterate("flexibles", "lst.append(chain+':'+resn+resi)", space=myspace)
                cmd.delete("flexibles")
                flex = [(myspace["lst"][0])]
                for i in range(len(myspace["lst"])-1):
                    if (myspace["lst"][i]) == (myspace["lst"][i+1]):
                        continue
                    else:
                        flex.append(myspace["lst"][i+1])
                for i in range(len(flex)):
                    self.form.listWidget.addItem(flex[i])
                    self.current_flexibles.append(flex[i])
                select_current_flexibles()
                self.form.radioButton_2.setChecked(True)
                self.form.radioButton_4.setChecked(True)
            else:
                set_statusline('ERROR : Could not find selection in PyMol')

        def clear_flexibles():
            cmd.delete("flexres")
            self.form.listWidget.clear()
            self.current_flexibles = []

        def delete_flexible():
            sel = self.form.listWidget.selectedItems()
            if not sel:
                return
            for item in sel:
                del_flex = self.form.listWidget.currentItem().text()
                self.current_flexibles.remove(del_flex)
                print("deleted: "+del_flex)
                self.form.listWidget.takeItem(self.form.listWidget.row(item))            
            select_current_flexibles()
                
        def synchronize_Radiobuttons_2():
            self.Buttongroup_2.setExclusive(False)
            self.form.radioButton_3.setChecked(self.form.radioButton.isChecked())
            self.form.radioButton_4.setChecked(self.form.radioButton_2.isChecked())
            self.Buttongroup_2.setExclusive(True)

        # Run page buildup
        self.form.comboBox.setCurrentText("")        
        self.form.comboBox_2.setCurrentText("")
        self.form.textBrowser_2.setHtml(receptor_text)
        import_charge_models()

        # Ligand page

        ligand_text = """
        <html><body> * Ligands are loaded from and are saved into the Ligand directory.<br>
         * If ligands are selected in .pdb format and their file in .pdbqt format is missing, the latter
        is automatically generated in the Ligand directory : either before docking or when adding them to the
        multirun list.<br> * When the multirun list is checked, ligands of the list are consecutively docked to the
        receptor in the selected docking configuration (rigid or flexible SC).<br>
         * Ligands have no more limits in torsion angles like it was in vina.
        </body></html>
        """

        def File_type_selector_buildup():
            File_type_list = ['pdbqt','pdb','sdf']
            self.form.comboBox_4.addItems(File_type_list)            

        def import_ligands():
            self.form.comboBox_2.clear()
            ligand_list =[]
            list_raw = glob(os.path.join(self.ligand_dir_path,"*."+self.form.comboBox_4.currentText()))
            for item in list_raw:
                ligand_name_raw = item.split("\\")[-1]
                ligand_name = ligand_name_raw.split(".")[0]
                ligand_list.append(ligand_name) 
            self.form.comboBox_2.addItems(ligand_list)

        def make_ligand_pdbqt(ligand):
            charge_model = self.form.comboBox_5.currentText()
            ligand_pdbqt = ligand.split(".")[0]+".pdbqt"
            print("Running openbable to create ligand.pdbqt :")
            command = 'call "%s" %s -O %s --partialcharge %s' % (self.openbabel_exe, ligand, ligand_pdbqt, charge_model)
            if self.form.groupBox_17.isChecked() == True:
                command = command+' -p %s' % (float(self.form.doubleSpinBox_4.value()))
            else :
                command = command+' -h'
            print(command)
            os.system(command)
            if os.path.isfile(ligand_pdbqt):
                set_statusline("Created %s" % ligand_pdbqt)
            else :
                set_statusline("ERROR when trying to create %s" % ligand_pdbqt)

        def clear_multirun_list():
            self.form.listWidget_2.clear()
            self.current_ligands = []

        def delete_ligands_from_list():
            sel = self.form.listWidget_2.selectedItems()
            if not sel: return
            for item in sel:
                del_lig = self.form.listWidget_2.currentItem().text()
                self.current_ligands.remove(del_lig)
                print("deleted: "+del_lig)
                self.form.listWidget_2.takeItem(self.form.listWidget_2.row(item))            

        def make_multirun_list():
            self.form.groupBox_14.setChecked(True)
            ligand = self.ligand_dir_path+self.form.comboBox_2.currentText()+"."+self.form.comboBox_4.currentText()
            if not (os.path.isfile(ligand.split(".")[0]+".pdbqt") or self.form.comboBox_4.currentText() == "pdbqt"):
                make_ligand_pdbqt(ligand)
            self.form.listWidget_2.addItem(self.form.comboBox_2.currentText())
            self.current_ligands.append(self.form.comboBox_2.currentText())

        def set_lig_directory_indicator(path):
            self.form.lineEdit_6.clear()
            self.form.lineEdit_6.insert(path)
            self.ligand_dir_path = path # refresh
            import_ligands() # new list when path changed
            
        # Run page buildup
        self.form.textBrowser_4.setHtml(ligand_text)
        File_type_selector_buildup()
        set_lig_directory_indicator(self.ligand_dir_path)
        import_ligands()

        # Smina page

        docking_text = """
        <html><body>Smina is a fork of vina including advanced scoring and minimization :<br>
         * Docking with Smina allows the use of user customed scoring tables. 
        The actual scoring table might be edited or new scoring tables can be created or loaded.<br>
         * A log-file containing all docking parameters is created if this option is checked.<br>
         * Results are selectable and ancient results can be loaded on the results page.<br>
         * All results can be post refined.
         </body></html>
        """

        def make_receptor_pdbqt(prot):
            cmd.save(prot+".pdb", prot)
            receptor = os.getcwd()+"\\"+prot+".pdbqt"
            charge_model = self.form.comboBox_5.currentText()
            print("running openbable to create %s.pdbqt :" % prot)
            command = 'call "%s" %s.pdb -O %s --partialcharge %s' % (self.openbabel_exe, prot, receptor, charge_model)
            if self.form.groupBox_17.isChecked() == True:
                command = command+' -p %s' % (float(self.form.doubleSpinBox_4.value()))
            else :
                command = command+' -h'
            print(command)
            os.system(command)
            if os.path.isfile(receptor):
                set_statusline("Created %s" % receptor)
            else :
                set_statusline("ERROR when trying to create %s" % receptor)            

        def use_custom_scoring():
            if self.form.groupBox_26.isChecked() == True:
                scoring_table = self.form.tableView
                scoring_model = ScoringTableModel(self.scoring_list)
                scoring_table.setModel(scoring_model)
                scoring_table.resizeColumnsToContents()

        def show_scoring_tables(path):
            self.scoring_table_dir_path = str(path) # refresh
            self.form.comboBox_6.clear()
            self.form.comboBox_7.clear()
            scoring_table_list = ["vina"] # current table is default
            list_raw = glob(os.path.join(self.scoring_table_dir_path,"*.scr"))
            for item in list_raw:
                scoring_table_name = item.split(".")[0].split("\\")[-1]
                scoring_table_list.append(scoring_table_name) 
            self.form.comboBox_6.addItems(scoring_table_list)
            self.form.comboBox_7.addItems(scoring_table_list)
            self.form.comboBox_7.addItems(["new"])
            change_current_score_table()
            
        def change_current_score_table():
            scoring_table_name = self.form.comboBox_6.currentText()
            if scoring_table_name == 'vina':
                self.scoring_list = vinascr
                use_custom_scoring()
                set_statusline('Default scoring table: vina')
            else :
                if os.path.isfile(os.path.join(self.scoring_table_dir_path,scoring_table_name+".scr")):
                    set_statusline('Reading scoring table: %s' % scoring_table_name)
                    self.scoring_list = []
                    filename = os.path.join(self.scoring_table_dir_path,scoring_table_name+".scr")
                    with open(filename,'r') as f:
                        lst = f.readlines()
                        for line in lst:
                            weight = line.split(' ')[0].translate(str.maketrans('','', "\n"))
                            potential = line.split(' ')[-1].translate(str.maketrans('','', "\n"))
                            self.scoring_list.append([weight, potential])
                    self.scoring_table_file = filename
                    use_custom_scoring()
                else:
                  set_statusline('ERROR : Could not find scoring table: %s' % scoring_table_name)  

        def select_scoring_table():
            selected_scoring_table = self.form.comboBox_7.currentText()
            if selected_scoring_table == 'vina':
                edit_scoring_list = vinascr
            elif selected_scoring_table == 'new':
                edit_scoring_list = ['','']
            else :
                edit_scoring_list = []
                if os.path.isfile(os.path.join(self.scoring_table_dir_path,selected_scoring_table+".scr")):                    
                    filename = os.path.join(self.scoring_table_dir_path,selected_scoring_table+".scr")
                    with open(filename,'r') as f:
                        lst = f.readlines()
                        for line in lst:
                            weight = line.split(' ')[0].translate(str.maketrans('','', "\n"))
                            potential = line.split(' ')[-1].translate(str.maketrans('','', "\n"))
                            edit_scoring_list.append([weight, potential]) 
                else:                    
                    set_statusline('ERROR : Could not find scoring table %s for editing' % selected_scoring_table)
            return edit_scoring_list

        def edit_score_table(): 
            data = select_scoring_table()
            selected_scoring_table = self.form.comboBox_7.currentText()
            editor = Editor(data, self.scoring_table_dir_path, selected_scoring_table)
            editor.show()
            show_scoring_tables(self.scoring_table_dir_path)
       
        def load_docked(outfile):
            if ( not os.path.isfile(outfile)):
                set_statusline('ERROR : Could not find %s in current directory' % outfile)
            cmd.load(outfile)

        def fill_score_list(outfile):
            filename = ""
            if outfile == "": # load _docked file
                filedialog = QtWidgets.QFileDialog()
                if filedialog.exec_():
                    filename = ''.join(map(str, filedialog.selectedFiles()))
                cmd.load(filename)
            else:
                filename = outfile
            lst = open(filename, 'r').readlines()
            score_list = []
            for line in lst:
                if 'MODEL' in line:
                    modnum = line.strip(" MODEL\n")
                if 'minimizedAffinity' in line:
                    modE = line.strip('REMAK minzedAfty\n') 
                    score_list.append(modnum+";"+modE)
            if filename.split('.')[0].rsplit('_', 1)[-1] != "docked": # not a '_docked.pdbqt' file 
                ligand_name = filename.translate(str.maketrans('\\','/','')).rsplit('/', 1)[-1].split('.')[0]
            else :
                ligand_name = filename.translate(str.maketrans('\\','/','')).rsplit('_', 1)[0].rsplit('/', 1)[-1] # normal outfile
            if self.firstrun == True: # only at first cyle of fill to get rid of tab1
                self.form.tabWidget_2.clear()
                self.firstrun = False
            score_table = QtWidgets.QTableWidget(parent = self.form.tabWidget_2)
            self.form.tabWidget_2.addTab(score_table,ligand_name)
            score_table.setGeometry(QtCore.QRect(0, 0, 450, 310))
            score_table.setColumnCount(3)
            score_table.setRowCount(len(score_list))
            score_table.setHorizontalHeaderLabels(["Model","Affinty","Select"])
            score_table.setMouseTracking(True)
            score_table.setEditTriggers(QtWidgets.QTableWidget.NoEditTriggers) # no doubleclick editing
            for i in range(len(score_list)):
                mod = QtWidgets.QTableWidgetItem(score_list[i].split(";")[0])
                aff = QtWidgets.QTableWidgetItem(score_list[i].split(";")[-1])
                check = QtWidgets.QTableWidgetItem("[ ]")
                mod.setTextAlignment(QtCore.Qt.AlignCenter)
                aff.setTextAlignment(QtCore.Qt.AlignCenter)
                check.setTextAlignment(QtCore.Qt.AlignCenter)
                score_table.setItem(i,0,mod)
                score_table.setItem(i,1,aff)
                score_table.setItem(i,2,check)
            score_table.cellClicked.connect(check_button_select) # this gets the current cell for selection

        def combine_flexres(flexout): # merges all residues of the same model into a multi-model file starting with model 1
            flexres_merged_name = "" # just in case no flexres is selected
            if os.path.isfile(flexout.split(".")[0]+"_merged.pdbqt"):
                os.remove(flexout.split(".")[0]+"_merged.pdbqt")
            if os.path.isfile(flexout):
                with open(flexout, 'r') as file:
                    oldfile = file.readlines()
                flexres_merged_name = flexout.split(".")[0]+"_merged.pdbqt"
                for line in oldfile:
                    if 'MODEL' in line:
                        actual_model = line
                        break
                first_model = 1
                with open(flexres_merged_name, 'a') as newfile:
                    newfile.write("MODEL "+str(first_model)+"\n")
                    for line in oldfile:
                        if 'MODEL' in line and line != actual_model:
                            newfile.write('ENDMDL\n')
                            first_model = first_model + 1
                            actual_model = line
                            newfile.write("MODEL "+str(first_model)+"\n")
                            continue
                        if 'MODEL' in line and line == actual_model:
                            continue        
                        if 'ENDMDL' in line:
                            continue
                        else :
                            newfile.write(line)
                    newfile.write('ENDMDL\n')
            return flexres_merged_name
           
        def smina_loop(ligand):
            prot = self.form.comboBox.currentText()
            if prot == "":
                set_statusline("ERROR : No structure selected")
                return
            else:
                receptor = os.getcwd()+"\\"+prot+".pdbqt"
                if ( not os.path.isfile(receptor)):
                    make_receptor_pdbqt(prot)
            ligand_name = ligand.split("/")[-1]
            config = os.getcwd()+"\\"+prot+"_config.txt"
            if self.form.groupBox_15.isChecked() == True: # add suffix to outfile
                ligand_name = ligand_name+"_"+self.form.lineEdit_7.text()
            outfile = os.getcwd()+"\\"+ligand_name+"_docked.pdbqt"
            if ( not os.path.isfile(config)): # check presence of smina config file for receptor 
                set_statusline('ERROR : Could not find %s_config.txt in current directory' % prot)
                return
            else:
                print("running smina :")
                receptor_wsl = "/mnt/c"+receptor.split(":")[-1].translate(str.maketrans('\\','/','')) # format receptor path for wsl
                config_wsl = "/mnt/c"+config.split(":")[-1].translate(str.maketrans('\\','/','')) # format config path for wsl
                outfile_wsl = "/mnt/c"+outfile.split(":")[-1].translate(str.maketrans('\\','/','')) # format outfile path for wsl
                smina = "/mnt/c"+self.smina_exe.split(":")[-1] # format smina_path for wsl
                ligand = "/mnt/c"+ligand.split(":")[-1] # format ligand_dir_path for wsl
                exhaust = self.form.spinBox.value()
                maxmodes = self.form.spinBox_2.value()
                command = 'call wsl %s -r %s -l %s.pdbqt --config %s -o %s --flex_hydrogens --exhaustiveness %s --num_modes %s ' % (smina,
                    receptor_wsl, ligand, config_wsl, outfile_wsl, exhaust, maxmodes)                 
                if self.form.radioButton_2.isChecked() == True: # run smina with flexibles
                    if self.current_flexibles != []:
                        flex_raw = []
                        for i in range(len(self.current_flexibles)):
                            A = self.current_flexibles[i].split(":")[0]
                            B = self.current_flexibles[i].split(":")[-1].strip(':ABCDEFGHIJKLMNOPQRSTUVWXYZ')
                            flex_raw.append(A+":"+B)
                        flexibles = str(flex_raw).strip('[]').translate(str.maketrans('','', "'")).translate(str.maketrans('','', " "))
                        flexout = os.getcwd()+"\\"+ligand_name+"_flexres.pdbqt"
                        flexout_wsl = "/mnt/c"+flexout.split(":")[-1].translate(str.maketrans('\\','/','')) # format outfile path for wsl
                        command = command+' --flexres %s --out_flex %s' % (flexibles, flexout_wsl)
                    else:
                        set_statusline("ERROR : No flexibles selected -> using rigid side chains")
                        self.Buttongroup_1.setExclusive(False)
                        self.form.radioButton.setChecked(True)
                        self.form.radioButton_2.setChecked(False)
                        self.Buttongroup_1.setExclusive(True)
                        print("Missing flexible residues switched to rigid side chains")
                        return
                if self.form.groupBox_25.isChecked() == True: # set seed
                    seed = int(self.form.doubleSpinBox_5.value())
                    command = command+' --seed %s' % (seed)
                if self.form.groupBox_26.isChecked() == True: # custom score table
                    if self.form.comboBox_6.currentText() == "vina":
                        None
                    else:
                        custscrtbl = "/mnt/c"+self.scoring_table_file.split(":")[-1] # format score_table_path for wsl
                        command = command+' --custom_scoring %s' % (custscrtbl)
                if self.form.checkBox.isChecked() == True: # write smina Logfile
                    if self.form.groupBox_15.isChecked() == True: # add suffix to logfile too !
                        logfile = receptor.split(".")[0]+'_'+ligand.split('/')[-1].split('.')[0]+"_"+self.form.lineEdit_7.text()+'.log'
                    else :
                        logfile = receptor.split(".")[0]+'_'+ligand.split('/')[-1].split('.')[0]+'.log'
                logfile_wsl = "/mnt/c"+logfile.split(":")[-1].translate(str.maketrans('\\','/','')) # format receptor path for wsl
                command = command+' --log %s' % (logfile_wsl)
                # print(command)
                os.system(command)
                if self.form.checkBox.isChecked() == True:
                    if os.path.isfile(logfile):
                        lst = []
                        with open(config, 'r') as f:
                            lst = f.readlines()
                        with open (logfile, 'a') as log_file:
                            log_file.write("Parameters : \n"+command+"\n"+str(lst))
                    else:
                        set_statusline("ERROR : Could not find "+logfile)
                load_docked(outfile)
                fill_score_list(outfile)
                if self.form.radioButton_2.isChecked() == True: 
                    flexres_merged = combine_flexres(flexout) # still has "_merged" in title !
                    if os.path.isfile(flexres_merged):
                        flexres_name = flexres_merged.rsplit("_",maxsplit=1)[0].split('\\')[-1]
                        if self.form.groupBox_15.isChecked() == True: # add suffix
                            flexres_name = flexres_name+'_'+self.form.lineEdit_7.text()
                        cmd.load(flexres_merged, flexres_name) 
                    else :
                        set_statusline('ERROR : Could not find %s in current directory' % flexres_merged)  

        def run_smina():                   
            if self.form.groupBox_14.isChecked() == True: # Multirun
                for i in range(self.form.listWidget_2.count()):
                    ligand = self.ligand_dir_path+self.form.listWidget_2.item(i).text()
                    smina_loop(ligand)
            else:
                if ( not os.path.isfile(self.ligand_dir_path+self.form.comboBox_2.currentText()+".pdbqt")):
                    ligand = self.ligand_dir_path+self.form.comboBox_2.currentText()+"."+self.form.comboBox_4.currentText()
                    make_ligand_pdbqt(ligand)
                ligand = self.ligand_dir_path+self.form.comboBox_2.currentText()
                smina_loop(ligand)

        def synchronize_Radiobuttons_1():
            self.Buttongroup_1.setExclusive(False)
            self.form.radioButton.setChecked(self.form.radioButton_3.isChecked())
            self.form.radioButton_2.setChecked(self.form.radioButton_4.isChecked())
            self.Buttongroup_1.setExclusive(True)            

        # Run page buildup
        self.form.textBrowser_5.setHtml(docking_text)
        show_scoring_tables(self.scoring_table_dir_path)

        # Results Page
            
        def set_checked(row):
            current_widget = self.form.tabWidget_2.currentWidget()
            score_table = current_widget
            checked = QtWidgets.QTableWidgetItem("[X]")
            checked.setTextAlignment(QtCore.Qt.AlignCenter)
            score_table.setItem(row,2,checked)
            modnr = score_table.item(row,0).text()
            ligand_name = self.form.tabWidget_2.tabText(self.form.tabWidget_2.currentIndex())
            if ligand_name.rsplit("_", 1)[-1] == ("minimized" or "localdocked" or "randomized"): # select a post-refined ligand
                None
            else :
                ligand = ligand_name+"_docked"
                cmd.select("selection", ligand, 1, 1, 0, modnr)    
                cmd.disable(ligand)
                if self.form.radioButton_2.isChecked() == True: # result has flexres
                    residues = self.form.tabWidget_2.tabText(self.form.tabWidget_2.currentIndex())+"_flexres"
                    cmd.select("selection", residues, 1, 1, 1, modnr)
                    cmd.disable(residues)
                    cmd.disable(ligand)
                cmd.create ("pose"+modnr+"_"+ligand_name, "selection", modnr, 1)
                cmd.delete("selection")
            
        def set_unchecked(row):
            check = QtWidgets.QTableWidgetItem("[ ]")
            check.setTextAlignment(QtCore.Qt.AlignCenter)
            current_widget = self.form.tabWidget_2.currentWidget()
            score_table = current_widget
            score_table.setItem(row,2,check)
            modnr = score_table.item(row,0).text()
            ligand_name = self.form.tabWidget_2.tabText(self.form.tabWidget_2.currentIndex())
            if ligand_name.rsplit("_", 1)[-1] == ("minimized" or "localdocked" or "randomized"): # select a post-refined ligand
                None
            else :
                cmd.delete ("pose"+modnr+"_"+ligand_name)            

        def check_button_select(row,column):
            current_widget = self.form.tabWidget_2.currentWidget() # score_table
            check_item = current_widget.item(row,2)
            if check_item.text() == "[X]":
                set_unchecked(row)
            else :
                set_checked(row)
            if self.form.checkBox_2.isChecked() == True: # update Postrefinement list if checked
                show_current_poses() 

        def get_docked_file():
            outfile = ""
            fill_score_list(outfile)

        def export_current_results():
            score_list = []
            sel = self.form.tabWidget_2.currentIndex()
            for i in range(self.form.tabWidget_2.widget(sel).rowCount()):
                score_list.append(self.form.tabWidget_2.widget(sel).item(i, 0).text()
                    +","+self.form.tabWidget_2.widget(sel).item(i, 1).text())
            fp = open(self.form.tabWidget_2.tabText(sel)+".csv",'w')
            print("model,Affinity", file=fp)
            for i in range(len(score_list)):
                print(score_list[i], file=fp)
            fp.close()
            set_statusline('Wrote results for docking as %s.csv file' % self.form.tabWidget_2.tabText(sel))

        def get_selected_poses():
            current_widget = self.form.tabWidget_2.currentWidget() # score_table
            ligand_name = self.form.tabWidget_2.tabText(self.form.tabWidget_2.currentIndex())
            ligand_poses_selected_list =[]
            for i in range(current_widget.rowCount()):
                check_item = current_widget.item(i,2)
                if check_item.text() == "[X]":
                    modnr = current_widget.item(i,0).text()
                    ligand_poses_selected_list.append([ligand_name,modnr])
            return ligand_poses_selected_list

        def export_pdbs():
            ligand_poses_selected_list = get_selected_poses()
            for i in range(len(ligand_poses_selected_list)):
                ligand = ligand_poses_selected_list[i][0]
                if ligand.rsplit("_", 1)[-1] == ("minimized" or "localdocked" or "randomized"): # select a post-refined ligand
                    print("writing: "+ligand+" to "+ligand+".pdb")
                    cmd.save(ligand+".pdb", ligand)
                else :
                    modnr = ligand_poses_selected_list[i][1]
                    print("writing: "+"pose"+modnr+"_"+ligand+" to "+"pose"+modnr+"_"+ligand+".pdb")
                    cmd.save("pose"+modnr+"_"+ligand+".pdb", "pose"+modnr+"_"+ligand)            

        # Post refinement Page

        refinement_text = """
        <html><body>Ligands selected for post-refinement show up in the 'Current Poses' Window:<br>
         * Multiple minimization options might be selected simultainously.<br>
         * Only minimization with rigid sidechains is possible up to now.
         </body></html>
        """

        def get_current_pose_selections():
            if self.form.checkBox_2.isChecked() == True:
                selected_poses_list = []
                for i in range(self.form.tabWidget_2.count()):
                    self.form.tabWidget_2.setCurrentIndex(i)
                    ligand_poses_selected_list = get_selected_poses()
                    for i in range(len(ligand_poses_selected_list)):
                        ligand = ligand_poses_selected_list[i][0]
                        if ligand.rsplit("_", 1)[-1] == ("minimized" or "localdocked" or "randomized"): # select a post-refined ligand
                            selected_poses_list.append(ligand)
                        else :
                            modnr = ligand_poses_selected_list[i][1]
                            selected_poses_list.append("pose"+modnr+"_"+ligand) 
            else:
                selected_poses_list = []
            return selected_poses_list

        def show_current_poses(): # calculates an updated list every time it is called
            self.form.listWidget_3.clear()
            self.current_poses_list = []
            if self.form.tabWidget_2.findChildren(QtWidgets.QTableWidget): # Only if results exist
                selected_poses_list = get_current_pose_selections()
                for i in range(len(selected_poses_list)):
                    self.current_poses_list.append(selected_poses_list[i])
            for i in range(len(self.loaded_poses_list)):
                self.current_poses_list.append(self.loaded_poses_list[i])
            for i in range(len(self.current_poses_list)):
                self.form.listWidget_3.addItem(self.current_poses_list[i])            

        def delete_ligand_flexres(ligand_pdbqt): #delete flexible residues from file
            with open(ligand_pdbqt, 'r') as file:
                ligand_file = file.readlines()
            os.remove(ligand_pdbqt)
            with open(ligand_pdbqt, 'a') as newfile:
                for line in ligand_file:
                    if 'TORSDOF' in line:
                        newfile.write(line)
                        break
                    else:
                        newfile.write(line)
            return
                    
        def minimize():
            if self.current_poses_list == []:
                set_statusline("ERROR : No poses selected")
                return
            else:
                if self.form.groupBox_30.isChecked() == False and self.form.checkBox_4.isChecked() == False and self.form.checkBox_5.isChecked() == False :
                    set_statusline("ERROR : Select post refinement method first")
                    return
                else :
                    for i in range(len(self.current_poses_list)):
                        ligand_name = self.current_poses_list[i]
                        cmd.save(ligand_name+".pdb", ligand_name) # create file just before use
                        ligand = os.getcwd()+"\\"+ligand_name
                        ligand_pdb = ligand+".pdb"
                        make_ligand_pdbqt(ligand_pdb)
                        os.remove(ligand_pdb) # cleanup
                        # check if ligand_file has flexres                    
                        ligand_pdbqt = ligand+".pdbqt"
                        with open(ligand_pdbqt, 'r') as file: 
                            ligand_file = file.readlines()
                            ligand_flexres = 0
                            for line in ligand_file:
                                if 'TORSDOF' in line:
                                    ligand_flexres = ligand_flexres +1
                        if ligand_flexres > 1:
                            delete_ligand_flexres(ligand_pdbqt)
                        prot = self.form.comboBox.currentText()
                        if prot == "":
                            set_statusline("ERROR : No structure selected")
                            return
                        else:
                            receptor = os.getcwd()+"\\"+prot+".pdbqt"
                            if ( not os.path.isfile(receptor)):
                                make_receptor_pdbqt(prot)
                            config = os.getcwd()+"\\"+prot+"_config.txt"
                            if ( not os.path.isfile(config)): # check presence of smina config file for receptor 
                                set_statusline('ERROR : Could not find %s_config.txt in current directory' % (prot))
                                return
                            config_wsl = "/mnt/c"+config.split(":")[-1].translate(str.maketrans('\\','/','')) # format config path for wsl
                            receptor_wsl = "/mnt/c"+receptor.split(":")[-1].translate(str.maketrans('\\','/','')) # format receptor path for wsl
                            ligand_wsl = "/mnt/c"+ligand.split(":")[-1].translate(str.maketrans('\\','/','')) # format ligand path for wsl
                            smina = "/mnt/c"+self.smina_exe.split(":")[-1] # format smina path for wsl
                            command = 'call wsl %s -r %s -l %s.pdbqt --config %s ' % (smina, receptor_wsl, ligand_wsl, config_wsl)
                            if self.form.groupBox_26.isChecked() == True: # custom score table
                                if self.form.comboBox_6.currentText() == "vina":
                                    None
                                else:
                                    custscrtbl = "/mnt/c"+self.scoring_table_file.split(":")[-1] # format score_table_path for wsl
                                    command = command+' --custom_scoring %s' % (custscrtbl)
                            if self.form.groupBox_30.isChecked() == True:
                                print("running smina for local docking around "+ligand_name+" :")
                                outfile = ligand+"_localdocked.pdbqt"
                                outfile_wsl = "/mnt/c"+outfile.split(":")[-1].translate(str.maketrans('\\','/','')) # format outfile path for wsl
                                autobox_buf = str(self.form.spinBox_3.value())
                                command = 'call wsl %s -r %s -l %s.pdbqt -o %s --local_only --autobox_ligand %s.pdbqt --autobox_add %s' % (smina, receptor_wsl, ligand_wsl, outfile_wsl, ligand_wsl, autobox_buf) 
                            if self.form.checkBox_5.isChecked() == True:
                                print("running smina to randomize pose: "+ligand_name)
                                outfile = ligand+"_randomize.pdbqt"
                                outfile_wsl = "/mnt/c"+outfile.split(":")[-1].translate(str.maketrans('\\','/','')) # format outfile path for wsl
                                command = command+' --randomized_only -o %s ' % (outfile_wsl)
                            if self.form.checkBox_4.isChecked() == True:
                                if self.form.checkBox_5.isChecked() == True or self.form.groupBox_30.isChecked() == True: # d'ont change their -out filename
                                    command = command+' --minimize '
                                else:
                                    print("running smina to minimize pose: "+ligand_name)
                                    outfile = ligand+"_minimized.pdbqt"
                                    outfile_wsl = "/mnt/c"+outfile.split(":")[-1].translate(str.maketrans('\\','/','')) # format outfile path for wsl
                                    command = command+' --minimize -o %s ' % (outfile_wsl)
                                approx_method = self.form.comboBox_8.currentText()
                                command = command+' --approximation %s ' % (approx_method)
                                if self.form.checkBox_3.isChecked() == True:
                                    command = command+' --accurate_line '
                                if self.form.checkBox_6.isChecked() == True:
                                    command = command+' --minimize_early_term '
                                if self.form.groupBox_31.isChecked() == True:
                                    minim_steps = str(self.form.spinBox_7.value())
                                    command = command+' --minimize_iters %s ' % (minim_steps)
                                if self.form.groupBox_33.isChecked() == True:
                                    aprox_factor = str(self.form.spinBox_8.value())
                                    command = command+' --factor %s ' % (aprox_factor)
                                if self.form.groupBox_32.isChecked() == True:
                                    force_cap = str(self.form.spinBox_9.value())
                                    command = command+' --force_cap %s ' % (force_cap)
                            if self.form.checkBox.isChecked() == True: # write smina Logfile
                                logfile = outfile.split(".")[0]+'.log'
                                logfile_wsl = "/mnt/c"+logfile.split(":")[-1].translate(str.maketrans('\\','/','')) # format outfile path for wsl
                                command = command+' --log %s' % (logfile_wsl)
                        # print(command)
                        os.system(command)
                        if self.form.checkBox_7.isChecked() == True:
                            if os.path.isfile(logfile):
                                with open (logfile, 'a') as log_file:
                                    log_file.write(command)
                            else:
                                set_statusline("ERROR : Could not find "+logfile)
                        load_docked(outfile)
                        fill_minimized_list(outfile)

        def score_poses():
            for i in range(len(self.current_poses_list)):
                if not os.path.isfile(os.getcwd()+"\\"+self.current_poses_list[i]+".pdbqt"):
                    ligand = os.getcwd()+"\\"+self.current_poses_list[i]+".pdb"
                    make_ligand_pdbqt(ligand)
                prot = self.form.comboBox.currentText()
                if prot == "":
                    set_statusline("ERROR : No structure selected")
                else:
                    receptor = os.getcwd()+"\\"+prot+".pdbqt"
                    if ( not os.path.isfile(receptor)):
                        make_receptor_pdbqt(prot)
                    outfile = os.getcwd()+"\\"+self.current_poses_list[i]+"_scored.pdbqt"
                    ligand = os.getcwd()+"\\"+self.current_poses_list[i]
                    ligand_name = self.current_poses_list[i]
                    if self.form.groupBox_15.isChecked() == True: # add suffix to outfile
                        outfile = outfile.split(".")[0]+"_"+self.form.lineEdit_7.text()+".pdbqt"
                    print("running smina to score pose "+ligand_name+" :")
                    smina = "/mnt/c"+self.smina_exe.split(":")[-1] # format smina path for wsl
                    receptor_wsl = "/mnt/c"+receptor.split(":")[-1].translate(str.maketrans('\\','/','')) # format receptor path for wsl
                    ligand_wsl = "/mnt/c"+ligand.split(":")[-1].translate(str.maketrans('\\','/','')) # format ligand path for wsl
                    outfile_wsl = "/mnt/c"+outfile.split(":")[-1].translate(str.maketrans('\\','/','')) # format outfile path for wsl
                    command = 'call wsl %s -r %s -l %s.pdbqt --score_only -o %s' % (smina, receptor_wsl, ligand_wsl, outfile_wsl)
                    if self.form.groupBox_26.isChecked() == True: # custom score table
                            if self.form.comboBox_6.currentText() == "vina":
                                None
                            else:
                                custscrtbl = "/mnt/c"+self.scoring_table_file.split(":")[-1] # format score_table_path for wsl
                                command = command+' --custom_scoring %s' % (custscrtbl)
                    # print(command)
                    os.system(command)
                    lst = []
                    with open(outfile, 'r') as f:
                        lst = f.readlines()
                        for line in lst:
                            if 'minimizedAffinity' in line:
                                modE = line.strip('REMAK minzedAfty\n') 
                                print("Affinty for "+ligand_name+" : "+modE)
                    time.sleep(2) # to let the results be printed
            
        def load_pose(): # (from anywhere)
            loaded_poses_list=[]
            filedialog = QtWidgets.QFileDialog()
            filename = filedialog.getOpenFileName(None, "load pose", os.getcwd(), 'docking_pose.pdb files (*.pdb)')
            filename = filename[0]
            cmd.load(filename)
            loaded_pose = filename.split(".")[0].split("/")[-1]
            self.loaded_poses_list.append(loaded_pose)
            if not os.path.isfile(loaded_pose+".pdb"): #  create pose.pdb file in cwd if loaded eleswhere
                cmd.save(loaded_pose+".pdb", loaded_pose)
            show_current_poses()

        # Run page buildup
        self.form.comboBox_8.addItems(approx_methods_list)
        self.form.textBrowser_6.setHtml(refinement_text)

        # Results Page

        def set_PRtable_checked(row):
            checked = QtWidgets.QTableWidgetItem("[X]")
            checked.setTextAlignment(QtCore.Qt.AlignCenter)
            self.form.tableWidget.setItem(row,2,checked)
            ligand_name = self.form.tableWidget.item(row,0).text()
            cmd.enable(ligand_name)
            
        def set_PRtable_unchecked(row):
            check = QtWidgets.QTableWidgetItem("[ ]")
            check.setTextAlignment(QtCore.Qt.AlignCenter)
            self.form.tableWidget.setItem(row,2,check)
            ligand_name = self.form.tableWidget.item(row,0).text()
            cmd.disable(ligand_name)

        def PRtable_check_button_select(row,column):
            check_item = self.form.tableWidget.item(row,2)
            if check_item.text() == "[X]":
                set_PRtable_unchecked(row)
            else:
                set_PRtable_checked(row) 

        def format_minimized_list():
            self.form.tableWidget.setColumnCount(3)
            self.form.tableWidget.setColumnWidth(1,100)
            self.form.tableWidget.setColumnWidth(2,80)
            self.form.tableWidget.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
            self.form.tableWidget.setHorizontalHeaderLabels(["Model","Affinty","Select"])
            self.form.tableWidget.setMouseTracking(True)
            self.form.tableWidget.setEditTriggers(QtWidgets.QTableWidget.NoEditTriggers) # no doubleclick editing
            self.form.tableWidget.cellClicked.connect(PRtable_check_button_select) # this gets the current cell for selection
            

        def fill_minimized_list(outfile):
            filename = ""
            if outfile == "": # load _minimized file
                filedialog = QtWidgets.QFileDialog()
                if filedialog.exec_():
                    filename = ''.join(map(str, filedialog.selectedFiles()))
                cmd.load(filename)
            else:
                filename = outfile
            lst = open(filename, 'r').readlines()
            ligand_name = filename.translate(str.maketrans('\\','/','')).rsplit('/', 1)[-1].split('.')[0] #  complete name
            for line in lst:
                if 'minimizedAffinity' in line:
                    modE = line.strip('REMAK minzedAfty\n') 
            new_row_number = self.form.tableWidget.rowCount()
            self.form.tableWidget.insertRow(new_row_number)
            check = QtWidgets.QTableWidgetItem("[X]")
            mod = QtWidgets.QTableWidgetItem(ligand_name)
            aff = QtWidgets.QTableWidgetItem(modE)
            mod.setTextAlignment(QtCore.Qt.AlignCenter)
            aff.setTextAlignment(QtCore.Qt.AlignCenter)
            check.setTextAlignment(QtCore.Qt.AlignCenter)
            self.form.tableWidget.setItem(new_row_number,0,mod)
            self.form.tableWidget.setItem(new_row_number,1,aff)
            self.form.tableWidget.setItem(new_row_number,2,check)
                    

        def get_minimized_file():
            outfile = ""
            fill_minimized_list(outfile)

        def export_current_minimization_results():
            score_list = []
            for i in range(self.form.tableWidget.rowCount()):
                score_list.append(self.form.tableWidget.item(i, 0).text()
                    +","+self.form.tableWidget.item(i, 1).text())
            fp = open("minimization.csv",'w')
            print("model,Affinity", file=fp)
            for i in range(len(score_list)):
                print(score_list[i], file=fp)
            fp.close()
            set_statusline('Wrote results of minimization as minimization.csv file') 

        def get_selected_minimization_poses():
            minimized_poses_selected_list =[]
            for i in range(self.form.tableWidget.rowCount()):
                check_item = self.form.tableWidget.item(i,2)
                if check_item.text() == "[X]":
                    ligand_name = self.form.tableWidget.item(i,0).text()
                    minimized_poses_selected_list.append([ligand_name])
            return minimized_poses_selected_list

        def export_minimization_pdbs():
            ligand_poses_selected_list = get_selected_minimization_poses()
            for i in range(len(ligand_poses_selected_list)):
                ligand = str(ligand_poses_selected_list[i]).strip("[]'")
                print("writing: "+ligand+" to "+ligand+".pdb")
                cmd.save(ligand+".pdb", ligand)

        # Run page buildup
        format_minimized_list()
                
        # callback bindings
        self.form.doubleSpinBox.valueChanged.connect(calculate_Box)
        self.form.doubleSpinBox_2.valueChanged.connect(calculate_Box)
        self.form.doubleSpinBox_3.valueChanged.connect(calculate_Box)
        self.form.spinBox_4.valueChanged.connect(calculate_Box)
        self.form.spinBox_5.valueChanged.connect(calculate_Box)
        self.form.spinBox_6.valueChanged.connect(calculate_Box)
        self.form.comboBox_3.currentIndexChanged.connect(set_Receptor)
        self.form.lineEdit_5.returnPressed.connect(Boxcenter_on_selection)
        self.form.comboBox_4.currentIndexChanged.connect(import_ligands)
        self.ligand_dir_location.textChanged.connect(set_lig_directory_indicator)
        self.form.groupBox_26.clicked.connect(use_custom_scoring)
        self.form.comboBox_6.currentIndexChanged.connect(change_current_score_table)
        self.form.comboBox_7.currentIndexChanged.connect(select_scoring_table)
        self.scoring_table_dir_location.textChanged.connect(show_scoring_tables)
        self.ligand_dir_location.textChanged.connect(import_ligands)
        self.form.checkBox_2.stateChanged.connect(show_current_poses)
        self.Buttongroup_1.buttonClicked.connect(synchronize_Radiobuttons_2)
        self.Buttongroup_2.buttonClicked.connect(synchronize_Radiobuttons_1)
         
        # launch on startup :
        import_objects()


        # ---------------------------------------------
        # All Button bindings :

        self.form.pushButton.clicked.connect(run_smina)
        self.form.pushButton_2.clicked.connect(get_smina_location)
        self.form.pushButton_3.clicked.connect(get_openbabel_location)
        self.form.pushButton_4.clicked.connect(get_ligand_dir_path)
        self.form.pushButton_5.clicked.connect(save_plugin_config_file)
        self.form.pushButton_6.clicked.connect(import_objects)
        self.form.pushButton_7.clicked.connect(import_flexibles)
        self.form.pushButton_8.clicked.connect(calculate_Box)
        self.form.pushButton_9.clicked.connect(hide_Box)
        self.form.pushButton_10.clicked.connect(save_smina_config_file)
        self.form.pushButton_11.clicked.connect(get_Boxcenter_selection)
        self.form.pushButton_12.clicked.connect(clear_flexibles)
        self.form.pushButton_13.clicked.connect(delete_flexible)
        self.form.pushButton_14.clicked.connect(make_multirun_list)
        self.form.pushButton_15.clicked.connect(clear_multirun_list)
        self.form.pushButton_16.clicked.connect(delete_ligands_from_list)
        self.form.pushButton_17.clicked.connect(get_docked_file)
        self.form.pushButton_18.clicked.connect(export_pdbs)
        self.form.pushButton_19.clicked.connect(export_current_results) 
        self.form.pushButton_20.clicked.connect(edit_score_table)
        self.form.pushButton_21.clicked.connect(get_scoring_table_dir_path)
        self.form.pushButton_22.clicked.connect(load_pose)
        self.form.pushButton_23.clicked.connect(score_poses)
        self.form.pushButton_24.clicked.connect(minimize)
        self.form.pushButton_25.clicked.connect(load_custom_config_file)
        self.form.pushButton_26.clicked.connect(get_minimized_file)
        self.form.pushButton_27.clicked.connect(export_current_minimization_results)
        self.form.pushButton_28.clicked.connect(export_minimization_pdbs)
        # ----------------------------------------------

class ScoringTableModel(QtCore.QAbstractTableModel):
    
    def __init__(self, datain, parent=None):
        QtCore.QAbstractTableModel.__init__(self, parent)
        self.arraydata = datain
        self._header = {0: 'Weight', 1: 'Potential',}
    
    def rowCount(self, parent):
        return len(self.arraydata)
    
    def columnCount(self, parent):
        return len(self.arraydata[0])
    
    def data(self, index, role):
        if not index.isValid():
            return None
        elif role != QtCore.Qt.DisplayRole:
            return None
        return self.arraydata[index.row()][index.column()]

    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
        if role == QtCore.Qt.DisplayRole:
            if orientation == QtCore.Qt.Horizontal:
                return self._header[section]
            else:
                return str(section)

class Editor(QtWidgets.QDialog):

    def __init__(self, data, scoring_table_dir_path, selected_scoring_table, parent=None):
        super(Editor, self).__init__(parent)
        self.setWindowTitle("Scoring Table Editor")
        self.setGeometry(50,50,500,500)
        Layout = QtWidgets.QVBoxLayout()
        self.scoring_table = QtWidgets.QTableWidget(parent=None)
        self.add_row_button = QtWidgets.QPushButton("Add Row")
        self.add_row_button.clicked.connect(self.edit_insert_row)
        self.delete_row_button = QtWidgets.QPushButton("Delete Row")
        self.delete_row_button.clicked.connect(self.edit_delete_row)
        self.button_box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Save
                                                | QtWidgets.QDialogButtonBox.Cancel)
        self.button_box.addButton(self.add_row_button, QtWidgets.QDialogButtonBox.ActionRole)
        self.button_box.addButton(self.delete_row_button, QtWidgets.QDialogButtonBox.ActionRole)
        self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject) 
        self.scoring_table.setGeometry(QtCore.QRect(10, 10, 450, 310))
        self.scoring_table.setColumnCount(2)
        self.scoring_table.horizontalHeader().setStretchLastSection(True)
        self.scoring_table.setHorizontalHeaderLabels(["Weight","Potential"])
        self.scoring_table.setRowCount(len(data))
        for i in range(len(data)):
            wei = QtWidgets.QTableWidgetItem(str(data[i][0]))
            pot = QtWidgets.QTableWidgetItem(str(data[i][1]))
            self.scoring_table.setItem(i,0,wei)
            self.scoring_table.setItem(i,1,pot)
        Layout.addWidget(self.scoring_table)
        Layout.addWidget(self.button_box)
        self.setLayout(Layout)
        if self.exec():
            print("saving file")
            rows = self.scoring_table.rowCount()
            new_filename = os.path.join(scoring_table_dir_path,selected_scoring_table+"_new.scr")
            with open(new_filename,'w') as f:
                for i in range(rows):
                    print(f'{self.scoring_table.item(i,0).text():13}{self.scoring_table.item(i,1).text()}', file=f)

    def edit_insert_row(self):
        row = self.scoring_table.rowCount()
        self.scoring_table.insertRow(row)
        print("Editor added row: ",row)

    def edit_delete_row(self):
        row = self.scoring_table.currentRow()
        self.scoring_table.removeRow(row)
        print("Editor deleted row: ", row)    
        
