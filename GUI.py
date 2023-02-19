from PyQt5.QtWidgets import QApplication, QMainWindow, QTextEdit, QAction, QFileDialog, QFontDialog, QSplitter, QWidget, QVBoxLayout, QLabel, QToolBar, QPushButton, QHBoxLayout, QMenu
from PyQt5.QtCore import Qt, QMimeData
from PyQt5.QtGui import QPixmap, QFont, QDrag, QTextImageFormat, QFontDatabase
import sys
import os

class TextEditor(QMainWindow):
    def __init__(self):
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

        # Connect buttons to their respective functions
        file_button.clicked.connect(self.open_file)
        save_button.clicked.connect(self.save_file)
        font_button.clicked.connect(self.change_font)
        image_button.clicked.connect(self.insert_image)

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
        print(file_name)  # Print the selected file name to the console
        if file_name:
            # Read file and set text in text edit widget
            with open(file_name, 'r') as file:
                self.text_edit.setPlainText(file.read())

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

if __name__ == '__main__':
    app = QApplication(sys.argv)
    editor = TextEditor()
    editor.show()
    sys.exit(app.exec_())
