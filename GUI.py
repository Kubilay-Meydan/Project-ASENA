from Bio import Entrez, SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import AlignIO, Phylo
from Bio.Align.Applications import MuscleCommandline
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from ipywidgets import Box, widgets
from PyQt5.QtWidgets import QDialogButtonBox, QMessageBox, QFormLayout, QLineEdit,QInputDialog, QDialog, QApplication, QMainWindow, QTextEdit, QAction, QFileDialog, QFontDialog, QSplitter, QWidget, QVBoxLayout, QLabel, QToolBar, QPushButton, QHBoxLayout, QMenu
from PyQt5.QtCore import Qt, QMimeData,QTimer
from PyQt5.QtGui import QPixmap, QFont, QDrag, QTextImageFormat, QFontDatabase
from IPython.display import display
import sys
import os
import ipywidgets as widgets
import ipyfilechooser as filechooser
from IPython.display import display, FileLinks
from matplotlib.backends.backend_agg import FigureCanvasAgg
import requests, random
from main import *

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
        self.random_label = QLabel(random.choice([str("Someday, we'll all be free"),str("5AM in Versailles"),str("No one man should have all that power"), str("Work well, you're important"), str("Up from the ashes, into the sky"), str("Make it all come to life")]), self)
        self.random_label.setAlignment(Qt.AlignCenter)
        font = QFont()
        font.setPointSize(15)
        font_id = QFontDatabase.addApplicationFont("font2.ttf")
        font_family = QFontDatabase.applicationFontFamilies(font_id)[0]
        font.setFamily(font_family)
        self.random_label.setStyleSheet('color : black')
        self.random_label.setFont(font)
        layout.addWidget(self.random_label)
        self.setLayout(layout)
        self.setFixedSize(400, 200)
        self.setWindowFlag(Qt.FramelessWindowHint)
        QTimer.singleShot(3000, self.close)


class TextEditor(QMainWindow):
    # Get the current working directory
    def __init__(self):
        super().__init__()
        font_id = QFontDatabase.addApplicationFont("font.ttf")
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
        blast = QAction('Blast', self)
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
        plasmid_editor = QAction('Plasmid Editor', self)
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
        multiple_alignment = QAction('From Alignment File', self)
        phylo_menu.addAction(multiple_alignment)
        phylo_button = QPushButton('Phylogeny')
        phylo_button.setMenu(phylo_menu)
        main_toolbar.addWidget(phylo_button)

        #connect button
        one_click.triggered.connect(self.one_click)
        create_alignment.triggered.connect(self.create_alignement)
        multiple_alignment.triggered.connect(self.from_alignement)
        gene_bank.triggered.connect(self.gene_bank)
        uni_prot.triggered.connect(self.uniprot)

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
        multiple_alignment = QPushButton('From Alignement File')
        gene_bank = QPushButton("get genebank info")
        uni_prot = QPushButton("get prot sequences")

        # Connect buttons to their respective functions
        file_button.clicked.connect(self.open_file)
        save_button.clicked.connect(self.save_file)
        font_button.clicked.connect(self.change_font)
        image_button.clicked.connect(self.insert_image)
        one_click_button.clicked.connect(self.one_click)
        create_alignment.clicked.connect(self.create_alignement)
        multiple_alignment.clicked.connect(self.from_alignement)
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

        triangle_pixmap = QPixmap('triangle.png')
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
        seq, ok_pressed = QInputDialog.getText(None, "Protein Statistics", "Enter protein sequence in capital letters:")
        # Only continue if the user clicked OK
        if ok_pressed:
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
            with open('output.fasta', "r") as f:
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
            aligned_window.setWindowTitle("Protein Sequences (output.fasta)")
            aligned_window.exec_()
        
    def gene_bank(button):
        seq, ok_pressed = QInputDialog.getText(None, "Id", "Enter Gene Id:")
        if seq == '':
            return "no file selected"
        if ok_pressed:
            # get data
            protein_analysis = get_genbank_info(seq)
        with open(str(seq), "r") as f:
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
        file_path, _ = QFileDialog.getOpenFileName(None, "Open File", "", "All Files (*);;Text Files (*.txt)")
        # Get the path to the selected file
        if file_path == '':
            return "no file selected"
        run_bmge_on_alignment(file_path, 'curated.fasta')
        fig = make_phylogenetic_tree_bof(file_path)
        # create a matplotlib figure canvas
        canvas = FigureCanvasAgg(fig)
        # create a widget box to hold the canvas and add it to the app
        canvas_widget = widgets.Output()
        with canvas_widget:
            display(canvas)
        box = widgets.Box(children=[canvas_widget])
        display(box)

    def create_alignement(button):
        file_path, _ = QFileDialog.getOpenFileName(None, "Open File", "", "All Files (*);;Text Files (*.txt)")
        if file_path != '':
            Align_muscle(file_path, 'aligned')
            with open('aligned', "r") as f:
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
        file_path, _ = QFileDialog.getOpenFileName(None, "Open File", "", "All Files (*);;Text Files (*.txt)")
        # Get the path to the selected file
        if file_path == '':
            return "no file selected"
        Align_muscle('output.fasta','aligned')
        run_bmge_on_alignment('aligned.fasta', 'curated.fasta')
        fig = make_phylogenetic_tree_bof('curated.fasta')
        # create a matplotlib figure canvas
        canvas = FigureCanvasAgg(fig)
        # create a widget box to hold the canvas and add it to the app
        canvas_widget = widgets.Output()
        with canvas_widget:
            display(canvas)
        box = widgets.Box(children=[canvas_widget])
        display(box)

    def save_file(self):
        # Open file dialog to select file to save to
        file_name, _ = QFileDialog.getSaveFileName(self, 'Save File', '', 'Rich Text Files (*.rtf *.docx);;Text Files (*.txt);;All Files (*)')
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

def run_gui():
    app = QApplication(sys.argv)
    popup = PopupWindow()  # create an instance of the PopupWindow class
    popup.exec_()  # show the popup window and wait for it to be closed
    editor = TextEditor()
    editor.show()
    sys.exit(app.exec_())
run_gui()