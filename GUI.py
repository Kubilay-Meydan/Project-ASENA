from Bio.SeqUtils.ProtParam import ProteinAnalysis
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from ipywidgets import widgets
from PyQt5.QtWidgets import QCheckBox, QGroupBox, QRadioButton, QMessageBox, QFormLayout, QLineEdit,QInputDialog, QDialog, QApplication, QMainWindow, QTextEdit, QAction, QFileDialog, QFontDialog, QSplitter, QWidget, QVBoxLayout, QLabel, QToolBar, QPushButton, QHBoxLayout, QMenu
from PyQt5.QtCore import Qt,QTimer
from PyQt5.QtGui import QPixmap, QFont, QTextImageFormat, QFontDatabase,QPalette,QBrush,QTextDocument
from PyQt5.QtCore import Qt, QEvent
from IPython.display import display
import sys
import os
import random
import subprocess
import platform
from src.main import *

class PopupWindow(QDialog):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("ASENA")
        layout = QVBoxLayout()
        self.label = QLabel("A.S.E.N.A", self)
        self.label.setAlignment(Qt.AlignCenter)
        font = self.label.font()
        font.setPointSize(53)
        self.label.setFont(font)
        self.label.setStyleSheet('color : black')
        layout.addWidget(self.label)
        self.sub_label = QLabel("Automatic Sequence Editing and Nomography Assistant", self)
        self.sub_label.setAlignment(Qt.AlignCenter)
        font.setPointSize(7)
        self.sub_label.setStyleSheet('color : black')
        self.sub_label.setFont(font)
        layout.addWidget(self.sub_label)
        self.random_label = QLabel(random.choice([str("🐭🐭🐭🐭🐭🐭🐭"),str("Someday, we'll all be free"),str("5AM in Versailles"),str("No one man should have all that power"), str("Work well, you're important"), str("Up from the ashes, into the sky"), str("Make it all come to life")]), self)
        self.random_label.setAlignment(Qt.AlignCenter)
        font = QFont()
        font.setPointSize(10)
        font_id = QFontDatabase.addApplicationFont("src/font.ttf")
        font_family = QFontDatabase.applicationFontFamilies(font_id)[0]
        font.setFamily(font_family)
        self.random_label.setStyleSheet('color : black')
        self.random_label.setFont(font)
        layout.addWidget(self.random_label)
        self.setLayout(layout)
        self.setFixedSize(400, 200)
        self.setWindowFlag(Qt.FramelessWindowHint)
        def read_settings():
            global tutorial
            global background 
            global games
            global glitch_mode
            tutorial = ''
            background = ''
            games = ''
            glitch_mode = ''
            with open('src/settings.txt') as f:
                contents = f.readlines()
                for line in contents:
                    if line.strip().startswith('tutorial'):
                        parsed = line.strip().split('= ')
                        if parsed[1] == 'True':
                            tutorial = True
                        else:
                            tutorial = False
                    if line.strip().startswith('email'):
                        parsed = line.strip().split('= ')
                        email = parsed[1]
                    if line.strip().startswith('background'):
                        parsed = line.strip().split('= ')
                        background = parsed[1]
                    if line.strip().startswith('games'):
                        parsed = line.strip().split('= ')
                        if parsed[1] == 'True':
                            games = True
                        else:
                            games = False
                    if line.strip().startswith('glitch_mode'):
                        parsed = line.strip().split('= ')
                        if parsed[1] == 'True':
                            glitch_mode = True
                        else:
                            glitch_mode = False    

        read_settings()
        QTimer.singleShot(2500, self.close)

class TextEditor(QMainWindow):
    # Get the current working directory
    def __init__(self):
        super().__init__()
        # Create a QPixmap object with the path to the image file
        pixmap = QPixmap('src/'+background)
        
        # Create a QPalette object and set its brush to the QPixmap object
        palette = self.palette()
        palette.setBrush(QPalette.Window, QBrush(pixmap))

        # Set the QPalette object as the background of the TextEditor
        self.setPalette(palette)
        font_id = QFontDatabase.addApplicationFont("src/font.ttf")
        font_family = QFontDatabase.applicationFontFamilies(font_id)[0]
        font = QFont(font_family, 10)
        QApplication.setFont(font)
        # Create main toolbar with buttons A, B, C, and Text Editor
        main_toolbar = self.addToolBar('Main Toolbar')
        main_toolbar.setMovable(False) # Disable toolbar movement

        
        tools_menu = QMenu('Tools', self)
        dna_to_rna_action = QAction('DNA to RNA', self)
        rna_to_dna_action = QAction('RNA to DNA', self)
        dna_to_dnac_action = QAction('DNA to DNAc', self)
        seq_find_action = QAction('Sequence find', self)
        patern_frequence_action = QAction('Patern frequence', self)
        prot_stats_action = QAction('Protein statistics', self)
        tools_menu.addAction(dna_to_rna_action)
        tools_menu.addAction(rna_to_dna_action)
        tools_menu.addAction(dna_to_dnac_action)
        tools_menu.addAction(seq_find_action)
        tools_menu.addAction(patern_frequence_action)
        tools_menu.addAction(prot_stats_action)
        tools_button = QPushButton('Tools')
        tools_button.setMenu(tools_menu)
        main_toolbar.addWidget(tools_button)

        prot_stats_action.triggered.connect(self.display_protein_stats)
        seq_find_action.triggered.connect(self.sequence_find_dialog)
        patern_frequence_action.triggered.connect(self.pattern_frequency_dialog)
        dna_to_rna_action.triggered.connect(self.dna_to_rna)
        rna_to_dna_action.triggered.connect(self.rna_to_dna)
        dna_to_dnac_action.triggered.connect(self.dna_to_dnac)

        #Create Widget for Prot with sub-buttons
        prot_menu = QMenu('Prot', self)
        blast = QAction('Blast (Soon)', self)
        prot_menu.addAction(blast)
        uni_prot = QAction('UniProt', self)
        prot_menu.addAction(uni_prot)
        prot_button = QPushButton('Prot')
        prot_button.setMenu(prot_menu)
        main_toolbar.addWidget(prot_button)

        #Create Widget for Gene with sub-buttons
        gene_menu = QMenu('Gene', self)
        gene_bank = QAction('Genbank Info', self)
        gene_menu.addAction(gene_bank)
        RNA_align = QAction('Rna score and alignment',self)
        gene_menu.addAction(RNA_align)
        plasmid_editor = QAction('Plasmid Editor (Soon)', self)
        gene_menu.addAction(plasmid_editor)
        gene_button = QPushButton('Gene')
        gene_button.setMenu(gene_menu)
        main_toolbar.addWidget(gene_button)

        #Create Widget for Phylo with sub-buttons
        phylo_menu = QMenu('Phylogeny', self)
        one_click = QAction('One Click', self)
        phylo_menu.addAction(one_click)
        create_alignment = QAction('Create Alignment', self)
        phylo_menu.addAction(create_alignment)
        alignement_curation = QAction('BGME Curation', self)
        phylo_menu.addAction(alignement_curation)
        multiple_alignement = QAction('From Alignment File', self)
        phylo_menu.addAction(multiple_alignement)
        phylo_button = QPushButton('Phylogeny')
        phylo_button.setMenu(phylo_menu)
        main_toolbar.addWidget(phylo_button)

        #connect button
        one_click.triggered.connect(self.one_click)
        create_alignment.triggered.connect(self.create_alignement)
        alignement_curation.triggered.connect(self.BMGE_curation)
        multiple_alignement.triggered.connect(self.from_alignement)
        gene_bank.triggered.connect(self.gene_bank)
        uni_prot.triggered.connect(self.uniprot)
        RNA_align.triggered.connect(self.display_rna_stats)
        settings_button = QPushButton('Settings', self)
        main_toolbar.addWidget(settings_button)
        settings_button.clicked.connect(self.show_settings_window)
        # Create widget for text editor
        editor_widget = QWidget()
        editor_layout = QVBoxLayout()

        # Create toolbar for text editor with buttons for File, Save, Font, Insert Image, and Close
        editor_toolbar = QToolBar()
        editor_toolbar.setObjectName("Editor Toolbar")
        editor_toolbar.setMovable(False) # Disable toolbar movement

        # Create buttons for text editor toolbar
        file_button = QPushButton('File')
        save_button = QPushButton('Save')
        font_button = QPushButton('Font')
        image_button = QPushButton('Insert Image')
        one_click_button = QPushButton('One_click')
        create_alignment = QPushButton('Create_Alignment')
        alignement_curation_button = QPushButton('Curate alignement with BMGE')
        gene_bank = QPushButton("get genebank info")
        uni_prot = QPushButton("get prot sequences")

        # Connect buttons to their respective functions
        file_button.clicked.connect(self.open_file)
        save_button.clicked.connect(self.save_file)
        font_button.clicked.connect(self.change_font)
        image_button.clicked.connect(self.insert_image)
        one_click_button.clicked.connect(self.one_click)
        create_alignment.clicked.connect(self.create_alignement)
        alignement_curation_button.clicked.connect(self.BMGE_curation)
        gene_bank.clicked.connect(self.gene_bank)
        uni_prot.clicked.connect(self.uniprot)

        # Add buttons to the text editor toolbar
        editor_toolbar.addWidget(file_button)
        editor_toolbar.addWidget(save_button)
        editor_toolbar.addWidget(font_button)
        editor_toolbar.addWidget(image_button)

        # Create text edit widget
        self.text_edit = QTextEdit()
        self.text_edit.setMinimumWidth(400)
        self.text_edit.setReadOnly(False)
        editor_layout.addWidget(editor_toolbar)
        editor_layout.addWidget(self.text_edit)

        # Set the editor widget layout and add to the main window
        editor_widget.setLayout(editor_layout)

        # Create a triangle logo widget
        triangle_widget = QWidget()
        triangle_layout = QVBoxLayout()

        triangle_pixmap = QPixmap('src/triangle.png')
        triangle_label = QLabel()
        triangle_label.setPixmap(triangle_pixmap)

        triangle_layout.addWidget(triangle_label, alignment=Qt.AlignRight)
        triangle_widget.setLayout(triangle_layout)

        # Create a splitter widget and add it to a vertical box layout with the text editor widget
        splitter_widget = QSplitter(Qt.Horizontal)
        splitter_widget.addWidget(triangle_widget)
        splitter_widget.addWidget(editor_widget)
        splitter_widget.setHandleWidth(8)
        splitter_widget.setStyleSheet("QSplitter::handle { background-color: gray; } ")
        splitter_widget.show()
        if background == 'Darkmode' or background == 'Dark_mode.png':
            main_toolbar.setStyleSheet('background-color: darkgray')
            tools_button.setStyleSheet('background-color: grey')
            prot_button.setStyleSheet('background-color: grey')
            gene_button.setStyleSheet('background-color: grey')
            phylo_button.setStyleSheet('background-color: grey')
            settings_button.setStyleSheet('background-color: grey')
            file_button.setStyleSheet('background-color: grey')
            save_button.setStyleSheet('background-color: grey')
            font_button.setStyleSheet('background-color: grey')
            image_button.setStyleSheet('background-color: grey')
            self.text_edit.setStyleSheet("background-color: darkgray;") 
        # Create a central widget to hold the splitter widget
        central_widget = QWidget()
        central_layout = QVBoxLayout()

        central_layout.addWidget(splitter_widget)
        central_widget.setLayout(central_layout)
        splitter_widget.moveSplitter(1,0)

        self.setCentralWidget(central_widget)

        # Set window title and dimensions
        self.setWindowTitle('Project ASENA')
        self.showMaximized()  # Set the window to take up the full screen on first open
    
    def show_settings_window(self):
        settings_window = QDialog()
        settings_window.setWindowTitle('Settings')
        settings_window.setGeometry(100, 100, 400, 300)

        tutorial_group = QGroupBox('Tutorial', settings_window)
        tutorial_layout = QHBoxLayout()
        tutorial_checkbox = QCheckBox()
        if tutorial is True:
            tutorial_checkbox.setChecked(True)
        else:
            tutorial_checkbox.setChecked(False)
        tutorial_checkbox.toggled.connect(self.update_tutorial)
        tutorial_layout.addWidget(tutorial_checkbox)
        tutorial_group.setLayout(tutorial_layout)

        color_scheme_group = QGroupBox('Color Scheme', settings_window)
        color_scheme_layout = QVBoxLayout()
        color_scheme_radio1 = QRadioButton('Light_mode')
        if background == 'Light_mode' or background == "Light_mode.png":
            color_scheme_radio1.setChecked(True)
        color_scheme_radio2 = QRadioButton('Dark_mode')
        if background == "Dark_mode" or background == "Dark_mode.png":
            color_scheme_radio2.setChecked(True)
        color_scheme_radio1.toggled.connect(self.update_color_scheme)
        color_scheme_radio2.toggled.connect(self.update_color_scheme)
        color_scheme_layout.addWidget(color_scheme_radio1)
        color_scheme_layout.addWidget(color_scheme_radio2)
        color_scheme_group.setLayout(color_scheme_layout)

        games_group = QGroupBox('Games', settings_window)
        games_layout = QHBoxLayout()
        games_checkbox = QCheckBox()
        if games is True:
            games_checkbox.setChecked(True)
        else:
            games_checkbox.setChecked(False)
        games_checkbox.toggled.connect(self.update_game)
        games_layout.addWidget(games_checkbox)
        games_group.setLayout(games_layout)


        main_layout = QVBoxLayout()
        main_layout.addWidget(tutorial_group)
        main_layout.addWidget(color_scheme_group)
        main_layout.addWidget(games_group)
        settings_window.setLayout(main_layout)
        # Add 'Save' button
        save_button = QPushButton('Save')
        save_button.clicked.connect(lambda: self.save_and_restart(settings_window))
        main_layout.addWidget(save_button)
        settings_window.setLayout(main_layout)
        settings_window.exec_()
    def update_tutorial(self, checked):
        self.tutorial = checked
        global tutorial
        self.tutorial = checked
        if tutorial is False:
            print('tutorial True now')
            tutorial = True
            input_file = "src/settings.txt"
            output_file = "src/text_temp.txt"

            with open(input_file, 'r') as file:
                lines = file.readlines()

            with open(output_file, 'w') as file:
                for line in lines:
                    if line.strip().startswith('tutorial'):
                        file.write('tutorial = True\n')
                    else:
                        file.write(line)

            # Replace the original file with the modified one
            import os
            os.remove(input_file)
            os.rename(output_file, input_file)
            return True
        if tutorial is True:
            print('tutorial False now')
            tutorial = False
            input_file = "src/settings.txt"
            output_file = "src/text_temp.txt"

            with open(input_file, 'r') as file:
                lines = file.readlines()

            with open(output_file, 'w') as file:
                for line in lines:
                    if line.strip().startswith('tutorial'):
                        file.write('tutorial = False\n')
                    else:
                        file.write(line)

            # Replace the original file with the modified one
            import os
            os.remove(input_file)
            os.rename(output_file, input_file)
            return False
    def save_and_restart(self, settings_window):
        # Save settings here if necessary
        # (Assuming you have already saved the settings in their respective methods)
        # Close the settings window
        settings_window.close()
        self.close()
        # Set the restart flag
        self.should_restart = True

    def update_color_scheme(self, checked):
        if checked:
            global background
            sender = self.sender()
            self.color_scheme = sender.text()
            option = self.color_scheme
            if background is not str(option):
                background = str(option)
                input_file = "src/settings.txt"
                output_file = "src/text_temp.txt"       
                with open(input_file, 'r') as file:
                    lines = file.readlines()        
                with open(output_file, 'w') as file:
                    for line in lines:
                        if line.strip().startswith('background'):
                            file.write('background = '+ str(option) + '.png\n')
                        else:
                            file.write(line)        
                # Replace the original file with the modified one
                import os
                os.remove(input_file)
                os.rename(output_file, input_file)
                return False

    def update_game(self, checked):
        global games
        self.game = checked
        if games is False:
            print('game True now')
            games = True
            input_file = "src/settings.txt"
            output_file = "src/text_temp.txt"

            with open(input_file, 'r') as file:
                lines = file.readlines()

            with open(output_file, 'w') as file:
                for line in lines:
                    if line.strip().startswith('games'):
                        file.write('games = True\n')
                    else:
                        file.write(line)

            # Replace the original file with the modified one
            import os
            os.remove(input_file)
            os.rename(output_file, input_file)
            return True
        if games is True:
            print('game False now')
            games = False
            input_file = "src/settings.txt"
            output_file = "src/text_temp.txt"

            with open(input_file, 'r') as file:
                lines = file.readlines()

            with open(output_file, 'w') as file:
                for line in lines:
                    if line.strip().startswith('games'):
                        file.write('games = False\n')
                    else:
                        file.write(line)

            # Replace the original file with the modified one
            import os
            os.remove(input_file)
            os.rename(output_file, input_file)
            return False
    def dna_to_rna(button):
        # Ouvre une fenêtre de dialogue pour entrer la séquence d'ADN
        dna_sequence, ok = QInputDialog.getText(None, 'DNA to RNA', 'Enter DNA sequence:')
        if ok:
            # Convertit la séquence d'ADN en ARN
            dna_sequence_upper = dna_sequence.upper()
            if not is_valid_enter_DNA(dna_sequence_upper):
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Critical)
                msg.setText("Invalid sequence")
                msg.setInformativeText("The sequence you entered contains invalid characters. Please enter a valid nucleic sequence")
                msg.setWindowTitle("Error")
                msg.exec_()
                return 
            rna_sequence = dna_sequence_upper.replace('T', 'U')
            rna_sequence = rna_sequence.replace('t', 'u')
            # Affiche la séquence d'ARN dans la fenêtre de dialogue
            msg = QMessageBox()
            msg.setWindowTitle("DNA to RNA")
            msg.setText(f"DNA sequence: {dna_sequence_upper}\nRNA sequence: {rna_sequence}")
            msg.exec_()

    def rna_to_dna(button):
        # Ouvre une fenêtre de dialogue pour entrer la séquence d'ARN
        rna_sequence, ok = QInputDialog.getText(None, 'RNA to DNA', 'Enter RNA sequence:')
        if ok:
            # Convertit la séquence d'ARN en ADN
            rna_sequence_upper = rna_sequence.upper()
            if not is_valid_enter_RNA(rna_sequence_upper):
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Critical)
                msg.setText("Invalid sequence")
                msg.setInformativeText("The sequence you entered contains invalid characters. Please enter a valid nucleic sequence")
                msg.setWindowTitle("Error")
                msg.exec_()
                return
            dna_sequence = rna_sequence_upper.replace('U', 'T')
            dna_sequence = dna_sequence.replace('u', 't')
            # Affiche la séquence d'ADN dans la fenêtre de dialogue
            msg = QMessageBox()
            msg.setWindowTitle("RNA to DNA")
            msg.setText(f"RNA sequence: {rna_sequence_upper}\nDNA sequence: {dna_sequence}")
            msg.exec_()

    def dna_to_dnac(button):
        # Ouvre une fenêtre de dialogue pour entrer la séquence d'ADN
        dna_sequence, ok = QInputDialog.getText(None, 'DNA to DNAc', 'Enter DNA sequence:')
        if ok:
            # Convertit la séquence d'ADN en ADN complémentaire
            dna_sequence_upper = dna_sequence.upper()
            if not is_valid_enter_DNA(dna_sequence_upper):
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Critical)
                msg.setText("Invalid sequence")
                msg.setInformativeText("The sequence you entered contains invalid characters. Please enter a valid nucleic sequence")
                msg.setWindowTitle("Error")
                msg.exec_()
                return          
            dnac_sequence = dna_sequence_upper.translate(str.maketrans("ATCG", "TAGC"))[::-1]
            dnac_sequence = dnac_sequence.translate(str.maketrans("atcg", "tagc"))[::-1]
            # Affiche la séquence d'ADN complémentaire dans la fenêtre de dialogue
            msg = QMessageBox()
            msg.setWindowTitle("DNA to DNAc")
            msg.setText(f"DNA sequence: {dna_sequence_upper}\nDNAc sequence: {dnac_sequence}")
            msg.exec_()

    def pattern_frequency_dialog(self):
            # Create a new dialog window
            dialog = QDialog(self)
            dialog.setWindowTitle('Pattern Frequency')
            dialog.resize(300, 100)

            # Create a form layout for the dialog
            layout = QFormLayout()

            # Add fields for the sequence and pattern length
            sequence_field = QLineEdit()
            layout.addRow('Sequence:', sequence_field)
            pattern_length_field = QLineEdit()
            layout.addRow('Pattern Length:', pattern_length_field)

            # Add a button to perform the search
            search_button = QPushButton('Recherche')
            layout.addRow(search_button)

            # Connect the button to the search function
            search_button.clicked.connect(lambda: self.perform_pattern_search(sequence_field.text(), int(pattern_length_field.text()), dialog))

            # Set the layout for the dialog and display it
            dialog.setLayout(layout)
            dialog.exec_()
    
    def perform_pattern_search(self, sequence, pattern_length, dialog):
        # Get all patterns of the specified length from the sequence
        patterns = [sequence[i:i+pattern_length] for i in range(len(sequence)-pattern_length+1)]

        # Count the occurrences of each pattern and store them in a dictionary
        counts = {}
        for i, pattern in enumerate(patterns):
            if pattern in counts:
                counts[pattern]['count'] += 1
                counts[pattern]['positions'].append(i+1)
            else:
                counts[pattern] = {'count': 1, 'positions': [i+1]}

        # Find the most frequent pattern and its occurrences
        most_frequent_pattern = max(counts, key=lambda x: counts[x]['count'])
        occurrences = counts[most_frequent_pattern]['positions']

        # Create a new dialog window to display the results
        result_dialog = QDialog(dialog)
        result_dialog.setWindowTitle('Résultats')
        result_dialog.resize(300, 100)

        # Create a layout for the dialog
        layout = QVBoxLayout()

        # Display the input sequence
        sequence_label = QLabel(f"Séquence d'entrée : {sequence}")
        layout.addWidget(sequence_label)

        # Display the most frequent pattern and its occurrences
        pattern_label = QLabel(f'Le motif le plus fréquent de longueur {pattern_length} est "{most_frequent_pattern}"')
        layout.addWidget(pattern_label)
        occurrences_label = QLabel(f'Il apparaît {len(occurrences)} fois aux positions suivantes : {occurrences}')
        layout.addWidget(occurrences_label)

        # Add a button to close the dialog
        close_button = QPushButton('Fermer')
        layout.addWidget(close_button)

        # Connect the button to close the dialog
        close_button.clicked.connect(result_dialog.close)

        # Set the layout for the dialog and display it
        result_dialog.setLayout(layout)
        result_dialog.exec_()

    def sequence_find_dialog(self):
        #Create a dialog box for sequence search
        dialog = QDialog(self)
        dialog.setWindowTitle('Sequence Find')
        layout = QFormLayout()

        # Create a text box for user to input sequence and pattern
        sequence_box = QLineEdit()
        layout.addRow('Sequence:', sequence_box)
        pattern_box = QLineEdit()
        layout.addRow('Pattern:', pattern_box)

        # Create a button to initiate pattern search
        search_button = QPushButton('Search')
        layout.addRow(search_button)

        # Create a label to display search results
        result_label = QLabel()
        layout.addRow('Result:', result_label)

        # Define function to search for pattern in sequence
        def search_pattern():
            # Get the sequence and pattern from user input
            sequence = sequence_box.text().upper()
            pattern = pattern_box.text().upper()

            # Search for the pattern in the sequence
            pat = search_a_pattern(sequence, pattern)[0]
            count = search_a_pattern(sequence,pattern)[3]
            positions= search_a_pattern(sequence,pattern)[2]

            # Display the number of occurrences and positions of the pattern in the sequence
            result_label.setText('Pattern "{}" found {} times at positions: {}'.format(pat, count, positions))

        # Connect the search button to the search_pattern function
        search_button.clicked.connect(search_pattern)

        dialog.setLayout(layout)
        dialog.exec_()

    def display_protein_stats(button):
    # Open a dialog window to get user input
        seq, ok_pressed = QInputDialog.getText(None, "Protein Statistics", "Enter a protein sequence of amino acids in capital letters:")
        # Only continue if the user clicked OK
        if ok_pressed:
            if seq == '':
                return False
            # Check if the input sequence is valid
            if not is_valid_sequence(seq):
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Critical)
                msg.setText("Invalid sequence")
                msg.setInformativeText("The sequence you entered contains invalid characters. Please enter a valid protein sequence.")
                msg.setWindowTitle("Error")
                msg.exec_()
                return

            # Calculate protein statistics
            protein_analysis = ProteinAnalysis(seq)
            sec_struct = protein_analysis.secondary_structure_fraction()
            stats = {
                "Molecular weight": f"{protein_analysis.molecular_weight():.2f} Da",
                "Aromaticity": f"{protein_analysis.aromaticity():.2f}",
                "Instability index (%)": f"{protein_analysis.instability_index():.2f}",
                "Isoelectric point (pH)": f"{protein_analysis.isoelectric_point():.2f}",
                "Secondary structure fraction (%) (helix)": f"{sec_struct[0]:.2f}",
                "Secondary structure fraction  (%) (sheet)": f"{sec_struct[1]:.2f}",
                "Secondary structure fraction (%) (coil)": f"{sec_struct[2]:.2f}"
            }

            # Display protein statistics in a new window
            stats_window = QDialog()
            stats_layout = QVBoxLayout()

            for stat, value in stats.items():
                label = QLabel(f"{stat}: {value}")
                stats_layout.addWidget(label)

            close_button = QPushButton("Close")
            close_button.clicked.connect(stats_window.close)
            stats_layout.addWidget(close_button)

            stats_window.setLayout(stats_layout)
            stats_window.setWindowTitle("Protein Statistics")
            stats_window.exec_()
    
    def display_rna_stats(button):
    # Create a custom QDialog to get user input
        input_dialog = QDialog()
        input_dialog.setWindowTitle("RNA Statistics")

        layout = QVBoxLayout()

        # Create input fields for RNA sequence 1, sequence 1 structure, and RNA sequence 2
        rna_seq1_label = QLabel("Enter RNA sequence 1 (in capital letters):")
        rna_seq1_input = QLineEdit()
        rna_struct1_label = QLabel("Enter sequence 1 structure: ex: 1-6")
        rna_struct1_input = QLineEdit()
        rna_seq2_label = QLabel("Enter RNA sequence 2 (in capital letters):")
        rna_seq2_input = QLineEdit()

        # Add input fields to the layout
        layout.addWidget(rna_seq1_label)
        layout.addWidget(rna_seq1_input)
        layout.addWidget(rna_struct1_label)
        layout.addWidget(rna_struct1_input)
        layout.addWidget(rna_seq2_label)
        layout.addWidget(rna_seq2_input)

        # Create buttons for submitting and canceling
        buttons_layout = QHBoxLayout()
        submit_button = QPushButton("Submit")
        cancel_button = QPushButton("Cancel")

        buttons_layout.addWidget(submit_button)
        buttons_layout.addWidget(cancel_button)
        layout.addLayout(buttons_layout)

        input_dialog.setLayout(layout)

        # Define button actions
        def on_submit():
            input_dialog.seq1 = rna_seq1_input.text()
            input_dialog.struct1 = rna_struct1_input.text()
            input_dialog.seq2 = rna_seq2_input.text()
            input_dialog.accept()
            
        def on_cancel():
            input_dialog.reject()

        submit_button.clicked.connect(on_submit)
        cancel_button.clicked.connect(on_cancel)
        def parse_ranges(s):
            ranges = s.split(", ")
            result = []
            for r in ranges:
                start, end = r.split("-")
                result.append((int(start), int(end)))
            return result
        result = input_dialog.exec_()
        if result == QDialog.Accepted:
            if not is_valid_enter_RNA(input_dialog.seq1) or not is_valid_enter_RNA(input_dialog.seq2):
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Critical)
                    msg.setText("Invalid sequence")
                    msg.setInformativeText("The sequence you entered contains invalid characters. Please enter a valid RNA sequence.")
                    msg.setWindowTitle("Error")
                    msg.exec_()
                    return
            if len(input_dialog.seq1) != len(input_dialog.seq2):
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Critical)
                    msg.setText("Invalid sequence")
                    msg.setInformativeText("The sequences you entered are not the same length.")
                    msg.setWindowTitle("Error")
                    msg.exec_()
                    return
            if len(input_dialog.struct1) == 0:
                msg = QMessageBox()
                msg.setIcon(QMessageBox.Critical)
                msg.setText("Invalid sequence")
                msg.setInformativeText("Please enter a structure for the first sequence")
                msg.setWindowTitle("Error")
                msg.exec_()
                return
            seq1 = input_dialog.seq1
            strstruct1 = input_dialog.struct1
            seq2 = input_dialog.seq2
            struct1 = parse_ranges(strstruct1)
            # Compute the edit distance and alignment
            output = edit_nested_plain(str(seq1), struct1, str(seq2))
            edit_distance = output[1, len(seq1), 1, len(seq2)]
            alignment = backtrace(str(seq1), str(seq2), output)

            # Display edit distance and alignment in a new window
            stats_window = QDialog()
            stats_layout = QVBoxLayout()

            edit_distance_label = QLabel(f"Edit distance: {edit_distance}")
            stats_layout.addWidget(edit_distance_label)

            alignment_label = QLabel(f"Alignment:\n{alignment[0]}")
            stats_layout.addWidget(alignment_label)
            alignment_label2 = QLabel(f"{alignment[1]}")
            stats_layout.addWidget(alignment_label2)

            close_button = QPushButton("Close")
            close_button.clicked.connect(stats_window.close)
            stats_layout.addWidget(close_button)

            stats_window.setLayout(stats_layout)
            stats_window.setWindowTitle("RNA Statistics")
            stats_window.exec_()

    def uniprot(button):
        # Create a QDialog with two input fields
        dialog = QDialog()
        dialog.setWindowTitle("get fasta")
        dialog.setModal(True)
        layout = QFormLayout()

        seq_input = QLineEdit()
        layout.addRow("Accession number(s) (comma separated no space):", seq_input)

        names_input = QLineEdit()
        layout.addRow("Protein name(s) (comma separated no space):", names_input)

        ok_button = QPushButton("OK")
        layout.addWidget(ok_button)
        dialog.setLayout(layout)

        # When the OK button is clicked, get the input values and close the dialog
        def ok_clicked():
            global acc
            global names
            acc = seq_input.text().split(",")
            names = names_input.text().split(",")
            dialog.close()
            # protein info and names
            sequences = all_sequences(acc)
            print('searching')
            write_fasta(sequences,names)
        # Open the dialog
        ok_button.clicked.connect(ok_clicked)
        dialog.exec_()
        print(seq_input.text())
        if seq_input.text() != '':
            with open('results/output.fasta', "r") as f:
                aligned_text = f.read()

            # Create a new window to display the aligned text
            aligned_window = QDialog(button)
            aligned_layout = QVBoxLayout()

            # Create a text edit widget and add the aligned text to it
            aligned_edit = QTextEdit()
            aligned_edit.setPlainText(aligned_text)
            aligned_edit.setReadOnly(True)
            aligned_layout.addWidget(aligned_edit)

            # Add a close button to the layout
            close_button = QPushButton("Close")
            close_button.clicked.connect(aligned_window.close)
            aligned_layout.addWidget(close_button)

            # Set the layout of the window and show it
            aligned_window.setLayout(aligned_layout)
            aligned_window.setWindowTitle("Protein Sequences (results/output.fasta)")
            aligned_window.exec_()

    def gene_bank(button):
        seq, ok_pressed = QInputDialog.getText(None, "Id", "Enter Gene Id:")
        if seq == '':
            return "no entry"
        if ok_pressed:
            # get data
            protein_analysis = get_genbank_info(seq)
        with open('results/'+str(seq), "r") as f:
            aligned_text = f.read()

        # Create a new window to display the aligned text
        aligned_window = QDialog(button)
        aligned_layout = QVBoxLayout()

        # Create a text edit widget and add the aligned text to it
        aligned_edit = QTextEdit()
        aligned_edit.setPlainText(aligned_text)
        aligned_edit.setReadOnly(True)
        aligned_layout.addWidget(aligned_edit)

        # Add a close button to the layout
        close_button = QPushButton("Close")
        close_button.clicked.connect(aligned_window.close)
        aligned_layout.addWidget(close_button)

        # Set the layout of the window and show it
        aligned_window.setLayout(aligned_layout)
        aligned_window.setWindowTitle("Genbank Info")
        aligned_window.exec_()
   
    def from_alignement(self):
        file_path, _ = QFileDialog.getOpenFileName(None, "Open File", "", "Text Files (*.txt); Fasta Files (*.fasta)")
        # Get the path to the selected file
        if file_path == '':
            return "no file selected"
        if not is_fasta(file_path):
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Critical)
                    msg.setText("Invalid sequence")
                    msg.setInformativeText("The file is not in fasta format")
                    msg.setWindowTitle("Error")
                    msg.exec_()
                    return 'error'
        run_bmge_on_alignment(file_path, 'results/curated.fasta')
        make_phylo_tree_newick('results/curated.fasta')

    def BMGE_curation(button):
        file_path, _ = QFileDialog.getOpenFileName(None, "Open File", "", "Text Files (*.txt); Fasta Files (*.fasta)")
        if file_path != '':

            if not is_fasta(file_path):
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Critical)
                    msg.setText("Invalid sequence")
                    msg.setInformativeText("The file is not in fasta format")
                    msg.setWindowTitle("Error")
                    msg.exec_()
                    return 'error'
            # Create a new window to display the aligned text
            run_bmge_on_alignment(file_path, 'results/curated.fasta')
            with open('results/curated.fasta', "r") as f:
                aligned_text = f.read()
            aligned_window = QDialog(button)
            aligned_layout = QVBoxLayout()

            # Create a text edit widget and add the aligned text to it
            aligned_edit = QTextEdit()
            aligned_edit.setPlainText(aligned_text)
            aligned_edit.setReadOnly(True)
            aligned_layout.addWidget(aligned_edit)

            # Add a close button to the layout
            close_button = QPushButton("Close")
            close_button.clicked.connect(aligned_window.close)
            aligned_layout.addWidget(close_button)

            # Set the layout of the window and show it
            aligned_window.setLayout(aligned_layout)
            aligned_window.setWindowTitle("Curated alignement")
            aligned_window.exec_()

    def create_alignement(button):
        file_path, _ = QFileDialog.getOpenFileName(None, "Open File", "", "Text Files (*.txt); Fasta Files (*.fasta)")
        if file_path != '':
            if not is_fasta(file_path):
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Critical)
                    msg.setText("Invalid sequence")
                    msg.setInformativeText("The file is not in fasta format")
                    msg.setWindowTitle("Error")
                    msg.exec_()
                    return 'error'
            Align_muscle(file_path, 'results/aligned')
            with open('results/aligned', "r") as f:
                aligned_text = f.read()

            # Create a new window to display the aligned text
            aligned_window = QDialog(button)
            aligned_layout = QVBoxLayout()

            # Create a text edit widget and add the aligned text to it
            aligned_edit = QTextEdit()
            aligned_edit.setPlainText(aligned_text)
            aligned_edit.setReadOnly(True)
            aligned_layout.addWidget(aligned_edit)

            # Add a close button to the layout
            close_button = QPushButton("Close")
            close_button.clicked.connect(aligned_window.close)
            aligned_layout.addWidget(close_button)

            # Set the layout of the window and show it
            aligned_window.setLayout(aligned_layout)
            aligned_window.setWindowTitle("Aligned Text")
            aligned_window.exec_()

    def open_file(self):
    # Open file dialog to select file
        file_name, _ = QFileDialog.getOpenFileName(self, 'Open File', '', 'Rich Text Files (*.rtf *.docx);;Text Files (*.txt);;All Files (*)')
          # Print the selected file name to the console
        if file_name:
            # Read file and set text in text edit widget
            with open(file_name, 'r') as file:
                self.text_edit.setPlainText(file.read())
    
    def one_click(button):   
        file_path, _ = QFileDialog.getOpenFileName(None, "Open File", "", "Text Files (*.txt); Fasta Files (*.fasta)")
        # Get the path to the selected file
        if file_path != '':
            if not is_fasta(file_path):
                    msg = QMessageBox()
                    msg.setIcon(QMessageBox.Critical)
                    msg.setText("Invalid sequence")
                    msg.setInformativeText("The file is not in fasta format")
                    msg.setWindowTitle("Error")
                    msg.exec_()
                    return 'error'
            Align_muscle(file_path,'results/aligned')
            run_bmge_on_alignment('results/aligned.fasta', 'results/curated.fasta')
            make_phylo_tree_newick('results/curated.fasta')

    def save_file(self):
        # Open file dialog to select file to save to
        file_name, _ = QFileDialog.getSaveFileName(self, 'Save File', '', 'Rich Text Files (*.rtf);;Word Document Files (*.docx);;Fasta File (*.fasta);;Text Files (*.txt)')
        if file_name:
            # Write text in text edit widget to file
            with open(file_name, 'w') as file:
                file.write(self.text_edit.toPlainText())

    def change_font(self):
        # Open font dialog to change font
        font, ok = QFontDialog.getFont()
        if ok:
            self.text_edit.setFont(font)

    def new_file(self):
        # Clear text in text edit widget and show the text edit widget and its toolbar
        self.text_edit.clear()
        self.text_edit.show()
        self.findChild(QToolBar, "Editor Toolbar").show()
        # Show the splitter widget
        self.centralWidget().show()
    
    def insert_image(self):
        fileName, _ = QFileDialog.getOpenFileName(self, "Open Image", "", "Image Files (*.png *.jpg *.bmp)")
        if fileName:
            imageFormat = QTextImageFormat()
            imageFormat.setName(fileName)
            self.text_edit.textCursor().insertImage(imageFormat)

def open_tutorial_pdf():
    pdf_file = "tutorial.pdf"
    try:
        if platform.system() == "Windows":
            os.startfile(pdf_file)
        elif platform.system() == "Darwin":
            subprocess.Popen(["open", pdf_file])
        else:
            subprocess.Popen(["xdg-open", pdf_file])
    except FileNotFoundError:
        print("The specified file was not found.")

def run_gui():
    global tutorial
    should_restart = True
    while should_restart:
        app = QApplication(sys.argv)
        popup = PopupWindow()  # create an instance of the PopupWindow class
        popup.exec_()  # show the popup window and wait for it to be closed

        editor = TextEditor()
        editor.show()
        if tutorial == True:
            open_tutorial_pdf()
            print('tutorial False now')
            tutorial = False
            input_file = "src/settings.txt"
            output_file = "src/text_temp.txt"

            with open(input_file, 'r') as file:
                lines = file.readlines()

            with open(output_file, 'w') as file:
                for line in lines:
                    if line.strip().startswith('tutorial'):
                        file.write('tutorial = False\n')
                    else:
                        file.write(line)

            # Replace the original file with the modified one
            import os
            os.remove(input_file)
            os.rename(output_file, input_file)
        exit_code = app.exec_()

        # Check if the application should restart
        should_restart = editor.should_restart if hasattr(editor, 'should_restart') else False

        # Clean up QApplication instance
        del app

    sys.exit(exit_code)

run_gui()
