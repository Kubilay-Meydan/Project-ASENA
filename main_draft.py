from Bio import Entrez, SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import AlignIO, Phylo
from Bio.Align.Applications import MuscleCommandline
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from ipywidgets import Box, widgets
from PyQt5.QtWidgets import QInputDialog, QDialog, QApplication, QMainWindow, QTextEdit, QAction, QFileDialog, QFontDialog, QSplitter, QWidget, QVBoxLayout, QLabel, QToolBar, QPushButton, QHBoxLayout, QMenu
from PyQt5.QtCore import Qt, QMimeData
from PyQt5.QtGui import QPixmap, QFont, QDrag, QTextImageFormat, QFontDatabase
from IPython.display import display
import sys
import os
import ipywidgets as widgets
import ipyfilechooser as filechooser
from IPython.display import display, clear_output
from matplotlib.backends.backend_agg import FigureCanvasAgg

entry = ''
def DNA_to_RNA(entry):
    ans = ''
    for char in entry.lower():
        if char == 'a':
            ans+=('a')
        if char == 'c':
            ans+=('c')
        if char == 'g':
            ans+=('g') 
        if char == 't':
            ans+=('u')
    return ans

Entrez.email = "your@email.com"  # replace with your email address

def get_genbank_info(gene_id):
    with open(str(gene_id)+'.txt', 'w') as f:
        handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        cds_sequence = []
        numberexon = 1
        numbercds = 1
        numberintron = 1
        f.write(record.name)
        f.write('\n')
        f.write('\n')
        f.write(record.description)
        f.write('\n')
        f.write('\n')
        for feature in record.features:
            if feature.type == "exon":
                f.write('Exon ')
                f.write(str(numberexon))
                f.write('\n')
                f.write(str(feature.location))                
                f.write('\n')
                f.write(str(feature.extract(record).seq))
                f.write('\n')
                f.write('\n')
                numberexon += 1
            if feature.type == "CDS":
                f.write('CDS ')
                f.write(str(numbercds))
                f.write('\n')
                f.write(str(feature.location))
                f.write('\n')
                f.write(str(feature.extract(record).seq))
                f.write('\n')
                f.write('\n')
                numbercds += 1
            if feature.type == "intron":
                f.write('Intron ')
                f.write(str(numberintron))
                f.write('\n')
                f.write(str(feature.location))
                f.write('\n')
                f.write(str(feature.extract(record).seq))
                f.write('\n')
                f.write('\n')
                numberintron += 1



def get_sequence(acc, seq_type='protein'):
    Entrez.email = 'your_email@example.com' # required by NCBI
    handle = Entrez.efetch(db='protein', id=acc, rettype='gb', retmode='text')
    record = SeqIO.read(handle, 'genbank')
    if seq_type == 'protein':
        return record.seq
    elif seq_type == 'DNA':
        return record.seq.transcribe()
    else:
        return None

def all_sequences(accessions):
    sequences = []
    for e in accessions:
        sequences.append(str(get_sequence(e)))
    return sequences

def write_fasta(sequences, names):
    with open("output.fasta", "w") as fasta:
        for seq, name in zip(sequences, names):
            fasta.write(">" + name + "\n")
            fasta.write(seq + "\n")



def Align_muscle(sequences,output):
# Align the sequences using muscle
    muscle_exe = "muscle.exe"
    cline = MuscleCommandline(muscle_exe,input=sequences, out=output, clw=True)
    cline()
    alignment = AlignIO.read(output, "clustal")
    return alignment



def make_phylogenetic_tree_bof(alignment_file):
    # Load the multiple sequence alignment from file
    alignment = AlignIO.read(alignment_file, "clustal")
    # Calculate the pairwise distances between sequences
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)
    # Construct the phylogenetic tree using the UPGMA method
    constructor = DistanceTreeConstructor(calculator, 'upgma')
    tree = constructor.build_tree(alignment)
    # Draw and show the tree
    Phylo.draw(tree)




# List of protein accession numbers
accessions = ['P66928','P14949','P10599','P34723','P0A4L3','P0AA25','P08629','P10639','P42115']
# List of protein names 
names = ['THIO_HELPY', 'THIO_BACSU', 'THIO_HUMAN', 'THIO_PENCH', 'THIO_LISMO','THIO_ECOLI', 'THIO_CHICK', 'THIO_MOUSE', 'THIO_NEUCR']


#
#       Gene id for human insulin, gets Exon and Intron details:
#

gene_id = 'L15440.1'
#get_genbank_info(gene_id)

#
#       Makes Phylo tree from accession:
#

#   Writes protein sequences to a fasta file using the uniprot accessions above and the names list given by user:
#sequences = all_sequences(accessions)
#write_fasta(sequences,names)

#   Uses MUSCLE alignement to make a Multiple Sequence Alignement File
#Align_muscle('output.fasta','aligned')
#   Makes a simple Phylo tree, no bootstrap yet
#make_phylogenetic_tree_bof('aligned')

def display_protein_stats():
    # Open a dialog window to get user input
    seq, ok_pressed = QInputDialog.getText(None, "Protein Statistics", "Enter protein sequence:")

    # Only continue if the user clicked OK
    if ok_pressed:
        # Calculate protein statistics
        protein_analysis = ProteinAnalysis(seq)
        stats = {
            "Molecular weight": protein_analysis.molecular_weight(),
            "Aromaticity": protein_analysis.aromaticity(),
            "Instability index": protein_analysis.instability_index(),
            "Isoelectric point": protein_analysis.isoelectric_point(),
            "Secondary structure fraction": protein_analysis.secondary_structure_fraction(),
        }

        # Display protein statistics in a new window
        stats_window = QDialog()
        stats_layout = QVBoxLayout()

        for stat, value in stats.items():
            label = QLabel(f"{stat}: {value}")
            stats_layout.addWidget(label)

        close_button = QPushButton("Fermer")
        close_button.clicked.connect(stats_window.close)
        stats_layout.addWidget(close_button)

        stats_window.setLayout(stats_layout)
        stats_window.setWindowTitle("Protein Statistics")
        stats_window.exec_()

class TextEditor(QMainWindow):
    # Get the current working directory
    def __init__(self):
        current_dir = os.getcwd()
        print(current_dir)
        super().__init__()
        font_id = QFontDatabase.addApplicationFont("font.ttf")
        font_family = QFontDatabase.applicationFontFamilies(font_id)[0]
        font = QFont(font_family, 10)
        QApplication.setFont(font)
        # Create main toolbar with buttons A, B, C, and Text Editor
        main_toolbar = self.addToolBar('Main Toolbar')
        main_toolbar.setMovable(False) # Disable toolbar movement

        #Create Widget File with sub-buttons
        file_menu = QMenu('File', self)
        new_action = QAction('New', self)
        new_action.setShortcut('Ctrl+N')
        open_action = QAction('Open', self)
        open_action.setShortcut('Ctrl+O')
        save_action = QAction('Save', self)
        save_action.setShortcut('Ctrl+S')
        save_as_action = QAction('Save As', self)
        file_menu.addAction(new_action)
        file_menu.addAction(open_action)
        file_menu.addAction(save_action)
        file_menu.addAction(save_as_action)
        file_button = QPushButton('File')
        file_button.setMenu(file_menu)
        main_toolbar.addWidget(file_button)
        # Connect buttons to their respective functions
        new_action.triggered.connect(self.new_file)
        open_action.triggered.connect(lambda: self.open_save_file_dialog('Open File', 'r'))
        save_action.triggered.connect(lambda: self.open_save_file_dialog('Save File', 'w'))
        save_as_action.triggered.connect(lambda: self.open_save_file_dialog('Save File As', 'w'))
        # Create dropdown menu for Tools with sub-buttons
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
        prot_stats_action.triggered.connect(display_protein_stats)

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
        multiple_alignment = QAction('Multiple Alignment from File', self)
        phylo_menu.addAction(multiple_alignment)
        phylo_button = QPushButton('Phylogeny')
        phylo_button.setMenu(phylo_menu)
        main_toolbar.addWidget(phylo_button)
        one_click.triggered.connect(self.one_click)
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

        # Connect buttons to their respective functions
        file_button.clicked.connect(self.open_file)
        save_button.clicked.connect(self.save_file)
        font_button.clicked.connect(self.change_font)
        image_button.clicked.connect(self.insert_image)
        one_click_button.clicked.connect(self.one_click)

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
        self.setWindowTitle('KN Gui')
        self.showMaximized()  # Set the window to take up the full screen on first open

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
        Align_muscle(file_path, 'aligned')
        fig = make_phylogenetic_tree_bof('aligned')
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
    editor = TextEditor()
    editor.show()
    sys.exit(app.exec_())
run_gui()
'''
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo.Consensus import bootstrap_consensus
import matplotlib.pyplot as plt
from Bio import Phylo

def make_phylogenetic_tree_bof(alignment_file, bootstraps=100):
    # Load the multiple sequence alignment from file
    alignment = AlignIO.read(alignment_file, "clustal")
    # Calculate the pairwise distances between sequences
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)
    # Construct the phylogenetic tree using the UPGMA method
    constructor = DistanceTreeConstructor(calculator, 'upgma')
    tree = constructor.build_tree(alignment)
    # Perform bootstrapping to estimate the confidence in the branches of the tree
    bootstrapped_trees = bootstrap_consensus(alignment_file, bootstraps,tree, 'strict_consensus')
    # Draw and show the tree
    Phylo.draw(bootstrapped_trees)
    plt.show()

make_phylogenetic_tree_bof('aligned')
'''