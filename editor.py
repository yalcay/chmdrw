from PySide6.QtWidgets import (QApplication, QMainWindow, QGraphicsView, 
                             QGraphicsScene, QVBoxLayout, QHBoxLayout, 
                             QPushButton, QWidget, QToolBar, QComboBox, QToolButton,
                             QFileDialog, QMessageBox, QFrame)
from PySide6.QtGui import QPen, QColor, QKeySequence, QAction, QIcon, QPainter, QPixmap
from PySide6.QtCore import Qt, QPointF, QSize, QPoint  # QPoint eklendi
import sys
import json
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import os
try:
    from rdkit.Chem import CDXMLWriter
    from rdkit.Chem import CDXMLReader
    HAS_CDX = True
except ImportError:
    HAS_CDX = False

class ChemDrawToolbar(QToolBar):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setOrientation(Qt.Vertical)
        self.setIconSize(QSize(24, 24))
        self.setFixedWidth(40)
        
        self.setStyleSheet("""
            QToolBar {
                background: #f0f0f0;
                border: 1px solid #c0c0c0;
                spacing: 2px;
            }
            QToolButton {
                border: 1px solid transparent;
                padding: 2px;
            }
            QToolButton:hover {
                background: #e0e0e0;
                border: 1px solid #a0a0a0;
            }
            QToolButton:checked {
                background: #d0d0d0;
                border: 1px solid #808080;
            }
        """)
        
        self.setup_tools()

    def setup_tools(self):
        # Seçim araçları grubu
        tools = [
            ("select", "Seçim Aracı", self.create_icon("select")),
            ("lasso", "Lasso Seçim", self.create_icon("lasso")),
            ("marquee", "Dikdörtgen Seçim", self.create_icon("marquee")),
            ("eraser", "Silgi", self.create_icon("eraser"))
        ]
        self.add_tool_group(tools)
        self.addSeparator()
        
        # Bağ araçları grubu
        tools = [
            ("single_bond", "Tek Bağ", self.create_icon("single_bond")),
            ("double_bond", "Çift Bağ", self.create_icon("double_bond")),
            ("triple_bond", "Üçlü Bağ", self.create_icon("triple_bond")),
            ("wedge_bond", "Kama Bağı", self.create_icon("wedge_bond")),
            ("hash_bond", "Kesik Bağ", self.create_icon("hash_bond"))
        ]
        self.add_tool_group(tools)
        self.addSeparator()
        
        # Atom araçları grubu
        tools = [
            ("C", "Karbon", self.create_text_icon("C")),
            ("O", "Oksijen", self.create_text_icon("O")),
            ("N", "Azot", self.create_text_icon("N")),
            ("H", "Hidrojen", self.create_text_icon("H")),
            ("S", "Kükürt", self.create_text_icon("S")),
            ("P", "Fosfor", self.create_text_icon("P"))
        ]
        self.add_tool_group(tools)
        self.addSeparator()
        
        # Şekil araçları grubu
        tools = [
            ("benzene", "Benzen", self.create_icon("benzene")),
            ("cyclopentane", "Siklopentan", self.create_icon("cyclopentane")),
            ("cyclohexane", "Siklohekzan", self.create_icon("cyclohexane"))
        ]
        self.add_tool_group(tools)

    def add_tool_group(self, tools):
        # Her araç için buton oluştur
        for tool_id, tooltip, icon in tools:
            button = QToolButton()
            button.setIcon(icon)
            button.setToolTip(tooltip)
            button.setCheckable(True)
            button.setIconSize(QSize(20, 20))
            button.clicked.connect(lambda checked, t=tool_id: self.tool_clicked(t))
            self.addWidget(button)

    def create_icon(self, name):
        # Temel simge oluşturma (örnek olarak basit şekiller)
        pixmap = QPixmap(24, 24)
        pixmap.fill(Qt.transparent)
        painter = QPainter(pixmap)
        painter.setRenderHint(QPainter.Antialiasing)
        
        try:
            if name == "select":
                painter.drawPolyline([QPoint(5, 5), QPoint(15, 15), QPoint(5, 20)])
            elif name == "lasso":
                points = [QPoint(5, 15), QPoint(10, 5), QPoint(15, 5), 
                         QPoint(20, 10), QPoint(15, 15), QPoint(10, 15)]
                for i in range(len(points)-1):
                    painter.drawLine(points[i], points[i+1])
            elif name == "marquee":
                painter.drawRect(5, 5, 14, 14)
            elif name == "eraser":
                painter.drawRect(5, 5, 14, 14)
                painter.drawLine(5, 5, 19, 19)
            elif name == "single_bond":
                painter.drawLine(5, 5, 19, 19)
            elif name == "double_bond":
                painter.drawLine(7, 7, 17, 17)
                painter.drawLine(5, 17, 17, 5)
            elif name == "triple_bond":
                painter.drawLine(7, 7, 17, 17)
                painter.drawLine(5, 12, 19, 12)
                painter.drawLine(5, 17, 17, 5)
            elif name == "wedge_bond":
                # Düzeltilmiş drawPolygon kullanımı
                points = [QPoint(5, 5), QPoint(19, 19), QPoint(15, 19)]
                painter.drawPolygon(points)  # Nokta listesini tek argüman olarak ver
            elif name == "hash_bond":
                for i in range(5):
                    painter.drawLine(5 + i*3, 5 + i*3, 8 + i*3, 8 + i*3)
            elif name == "benzene":
                painter.drawEllipse(4, 4, 16, 16)
            elif name == "cyclopentane":
                points = []
                from math import cos, sin, pi
                for i in range(5):
                    angle = i * 2 * pi / 5
                    x = 12 + 8 * cos(angle)
                    y = 12 + 8 * sin(angle)
                    points.append(QPoint(int(x), int(y)))
                # Düzeltilmiş drawPolygon kullanımı
                painter.drawPolygon(points)
            elif name == "cyclohexane":
                points = []
                from math import cos, sin, pi
                for i in range(6):
                    angle = i * 2 * pi / 6
                    x = 12 + 8 * cos(angle)
                    y = 12 + 8 * sin(angle)
                    points.append(QPoint(int(x), int(y)))
                # Düzeltilmiş drawPolygon kullanımı
                painter.drawPolygon(points)
        finally:
            painter.end()
        
        return QIcon(pixmap)

    def create_text_icon(self, text):
        pixmap = QPixmap(24, 24)
        pixmap.fill(Qt.transparent)
        painter = QPainter(pixmap)
        painter.setRenderHint(QPainter.Antialiasing)
        
        # Text ayarları
        font = painter.font()
        font.setPointSize(12)
        font.setBold(True)
        painter.setFont(font)
        
        # Metni merkeze çiz
        painter.drawText(pixmap.rect(), Qt.AlignCenter, text)
        painter.end()
        return QIcon(pixmap)

    def tool_clicked(self, tool_id):
        # Diğer butonların seçimini kaldır
        for button in self.findChildren(QToolButton):
            if button.isChecked() and button.toolTip() != tool_id:
                button.setChecked(False)
        
        # Ana pencereye sinyal gönder
        if isinstance(self.parent(), MoleculeEditor):
            if tool_id in ['C', 'O', 'N', 'H', 'S', 'P']:
                self.parent().set_current_element(tool_id)
            elif tool_id in ['single_bond', 'double_bond', 'triple_bond']:
                bond_type = {'single_bond': 'Tek Bağ',
                           'double_bond': 'Çift Bağ',
                           'triple_bond': 'Üçlü Bağ'}[tool_id]
                self.parent().set_bond_type(bond_type)
            elif tool_id in ['select', 'lasso', 'marquee', 'eraser']:
                self.parent().scene.current_tool = tool_id

class MoleculeDrawingScene(QGraphicsScene):
    def __init__(self):
        super().__init__()
        self.atoms = []
        self.bonds = []
        self.drawing = False
        self.last_point = None
        self.current_element = 'C'
        self.current_bond_type = 'single'
        self.selected_item = None
        self.undo_stack = []
        self.redo_stack = []
        self.current_tool = 'select'
        self.selection_path = None
        self.selection_items = []

    def add_to_undo_stack(self, action):
        self.undo_stack.append(action)
        self.redo_stack.clear()  # Yeni bir eylem yapıldığında redo stack'i temizle

    def undo(self):
        if not self.undo_stack:
            return
        
        action = self.undo_stack.pop()
        if action['type'] == 'add_atom':
            self.remove_atom(action['data'])
            self.redo_stack.append(action)
        elif action['type'] == 'add_bond':
            self.remove_bond(action['data'])
            self.redo_stack.append(action)
        elif action['type'] == 'remove_atom':
            atom = action['data']
            self.atoms.append(atom)
            self.draw_atom(atom)
            self.redo_stack.append({'type': 'add_atom', 'data': atom})
        elif action['type'] == 'remove_bond':
            bond = action['data']
            self.bonds.append(bond)
            self.draw_bond(bond)
            self.redo_stack.append({'type': 'add_bond', 'data': bond})

    def redo(self):
        if not self.redo_stack:
            return
        
        action = self.redo_stack.pop()
        if action['type'] == 'add_atom':
            atom = action['data']
            self.atoms.append(atom)
            self.draw_atom(atom)
            self.undo_stack.append(action)
        elif action['type'] == 'add_bond':
            bond = action['data']
            self.bonds.append(bond)
            self.draw_bond(bond)
            self.undo_stack.append(action)
        elif action['type'] == 'remove_atom':
            self.remove_atom(action['data'])
            self.undo_stack.append(action)
        elif action['type'] == 'remove_bond':
            self.remove_bond(action['data'])
            self.undo_stack.append(action)

    def remove_atom(self, atom):
        if atom in self.atoms:
            # Atomla ilişkili görsel öğeleri kaldır
            for item in atom.get('visual_items', []):
                self.removeItem(item)
            # Atomu listeden çıkar
            self.atoms.remove(atom)
            # Atomla ilişkili bağları kaldır
            bonds_to_remove = []
            for bond in self.bonds:
                if bond['start'] == atom['position'] or bond['end'] == atom['position']:
                    bonds_to_remove.append(bond)
            for bond in bonds_to_remove:
                self.remove_bond(bond)

    def remove_bond(self, bond):
        if bond in self.bonds:
            # Bağla ilişkili görsel öğeleri kaldır
            for item in bond.get('visual_items', []):
                self.removeItem(item)
            # Bağı listeden çıkar
            self.bonds.remove(bond)

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            pos = event.scenePos()
            
            if self.current_tool == 'select':
                self.select_item_at_position(pos)
            elif self.current_tool == 'lasso':
                self.start_lasso_selection(pos)
            elif self.current_tool == 'marquee':
                self.start_marquee_selection(pos)
            elif self.current_tool == 'eraser':
                self.erase_at_position(pos)
            else:  # Atom ve bağ çizimi için
                self.drawing = True
                atom = {
                    'position': pos,
                    'element': self.current_element,
                    'id': len(self.atoms)
                }
                self.atoms.append(atom)
                self.draw_atom(atom)
                self.add_to_undo_stack({'type': 'add_atom', 'data': atom})
                self.last_point = pos

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.LeftButton:
            if self.drawing and self.last_point:
                pos = event.scenePos()
                if (pos - self.last_point).manhattanLength() > 10:
                    bond = {
                        'start': self.last_point,
                        'end': pos,
                        'type': self.current_bond_type
                    }
                    self.bonds.append(bond)
                    self.draw_bond(bond)
                    self.add_to_undo_stack({'type': 'add_bond', 'data': bond})
            
            self.drawing = False
            self.last_point = None

    def mouseMoveEvent(self, event):
        pos = event.scenePos()
        
        if self.current_tool == 'lasso' and self.drawing:
            self.update_lasso_selection(pos)
        elif self.current_tool == 'marquee' and self.drawing:
            self.update_marquee_selection(pos)
        elif self.current_tool == 'eraser':
            self.erase_at_position(pos)
        elif self.drawing and self.last_point:
            self.update_drawing(pos)

    def start_lasso_selection(self, pos):
        self.drawing = True
        self.selection_path = [pos]
        # Önceki seçimi temizle
        self.clearSelection()

    def update_lasso_selection(self, pos):
        if self.selection_path:
            self.selection_path.append(pos)
            # Geçici çizgiyi güncelle
            self.update()

    def start_marquee_selection(self, pos):
        self.drawing = True
        self.selection_start = pos
        self.selection_rect = None
        # Önceki seçimi temizle
        self.clearSelection()

    def update_marquee_selection(self, pos):
        if self.selection_start:
            # Seçim dikdörtgenini güncelle
            rect = QRectF(self.selection_start, pos).normalized()
            if self.selection_rect:
                self.removeItem(self.selection_rect)
            self.selection_rect = self.addRect(rect, QPen(Qt.DashLine))

    def erase_at_position(self, pos):
        items = self.items(pos)
        for item in items:
            # Atom veya bağ silme işlemi
            self.removeItem(item)
            # İlgili veri yapılarından da sil
            self.remove_atom_or_bond(item)

    def draw_atom(self, atom):
        """Atomu çiz"""
        pos = atom['position']
        # Atom çemberi
        circle = self.addEllipse(pos.x()-3, pos.y()-3, 6, 6, 
                               QPen(QColor('black')), QColor('black'))
        # Element etiketi
        text = self.addText(atom['element'])
        text.setPos(pos.x()+5, pos.y()-10)
        atom['visual_items'] = [circle, text]

    def draw_bond(self, bond):
        """Bağı çiz"""
        start, end = bond['start'], bond['end']
        if bond['type'] == 'single':
            line = self.addLine(start.x(), start.y(), end.x(), end.y(), 
                              QPen(QColor('black')))
            bond['visual_items'] = [line]
        elif bond['type'] == 'double':
            # Çift bağ çizimi
            offset = 2
            dx = end.x() - start.x()
            dy = end.y() - start.y()
            length = (dx*dx + dy*dy)**0.5
            if length == 0:
                return
            nx = -dy/length * offset
            ny = dx/length * offset
            
            line1 = self.addLine(start.x()+nx, start.y()+ny, 
                               end.x()+nx, end.y()+ny, 
                               QPen(QColor('black')))
            line2 = self.addLine(start.x()-nx, start.y()-ny, 
                               end.x()-nx, end.y()-ny, 
                               QPen(QColor('black')))
            bond['visual_items'] = [line1, line2]
        elif bond['type'] == 'triple':
            # Üçlü bağ çizimi
            offset = 3
            dx = end.x() - start.x()
            dy = end.y() - start.y()
            length = (dx*dx + dy*dy)**0.5
            if length == 0:
                return
            nx = -dy/length * offset
            ny = dx/length * offset
            
            line1 = self.addLine(start.x(), start.y(), 
                               end.x(), end.y(), 
                               QPen(QColor('black')))
            line2 = self.addLine(start.x()+nx, start.y()+ny, 
                               end.x()+nx, end.y()+ny, 
                               QPen(QColor('black')))
            line3 = self.addLine(start.x()-nx, start.y()-ny, 
                               end.x()-nx, end.y()-ny, 
                               QPen(QColor('black')))
            bond['visual_items'] = [line1, line2, line3]

    def update_drawing(self, pos):
        """Çizim sırasında geçici bağ göster"""
        # Önceki geçici bağları temizle
        items = self.items()
        for item in items:
            if hasattr(item, 'isTemporary') and item.isTemporary:
                self.removeItem(item)
        
        # Geçici bağ çiz
        if self.last_point:
            temp_line = self.addLine(self.last_point.x(), self.last_point.y(),
                                   pos.x(), pos.y(), QPen(QColor('gray')))
            temp_line.isTemporary = True

    def select_item_at_position(self, pos):
        """Tıklanan pozisyondaki öğeyi seç"""
        items = self.items(pos)
        if items:
            if self.selected_item:
                # Önceki seçimi temizle
                if isinstance(self.selected_item, list):
                    for item in self.selected_item:
                        if hasattr(item, 'defaultPen'):
                            item.setPen(item.defaultPen)
                else:
                    if hasattr(self.selected_item, 'defaultPen'):
                        self.selected_item.setPen(self.selected_item.defaultPen)
            
            # Yeni seçimi yap
            selected = items[0]
            self.selected_item = selected
            if not hasattr(selected, 'defaultPen'):
                selected.defaultPen = selected.pen()
            selected.setPen(QPen(QColor('red')))

    def remove_atom_or_bond(self, item):
        """Atom veya bağı sil"""
        # Atomları kontrol et
        for atom in self.atoms:
            if item in atom.get('visual_items', []):
                self.remove_atom(atom)
                self.add_to_undo_stack({'type': 'remove_atom', 'data': atom})
                return

        # Bağları kontrol et
        for bond in self.bonds:
            if item in bond.get('visual_items', []):
                self.remove_bond(bond)
                self.add_to_undo_stack({'type': 'remove_bond', 'data': bond})
                return

class MoleculeEditor(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Molekül Editörü")
        self.setGeometry(100, 100, 1000, 600)

        # Ana widget ve layout
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QHBoxLayout(main_widget)

        # Sol taraftaki araç çubuğu
        self.chem_toolbar = ChemDrawToolbar(self)
        main_layout.addWidget(self.chem_toolbar)

        # Sağ taraf için konteyner
        right_container = QWidget()
        right_layout = QVBoxLayout(right_container)

        # Scene ve view
        self.scene = MoleculeDrawingScene()
        self.view = QGraphicsView(self.scene)
        right_layout.addWidget(self.view)

        # Sağ konteyner ekle
        main_layout.addWidget(right_container)

        # Menü oluştur
        self.create_menu()

        # Çizim alanı ayarları
        self.scene.setSceneRect(0, 0, 800, 500)

        # Kısayol tuşları
        self.setup_shortcuts()

        # Durum çubuğu
        self.statusBar().showMessage('Hazır')

    def create_menu(self):
        menubar = self.menuBar()
        
        # Dosya menüsü
        file_menu = menubar.addMenu('Dosya')
        
        new_action = QAction('Yeni', self)
        new_action.setShortcut('Ctrl+N')
        new_action.triggered.connect(self.clear_scene)
        
        save_action = QAction('Kaydet', self)
        save_action.setShortcut('Ctrl+S')
        save_action.triggered.connect(self.save_molecule)
        
        load_action = QAction('Aç', self)
        load_action.setShortcut('Ctrl+O')
        load_action.triggered.connect(self.load_molecule)
        
        file_menu.addAction(new_action)
        file_menu.addAction(save_action)
        file_menu.addAction(load_action)

        # Düzenle menüsü
        edit_menu = menubar.addMenu('Düzenle')
        
        undo_action = QAction('Geri Al', self)
        undo_action.setShortcut('Ctrl+Z')
        undo_action.triggered.connect(self.scene.undo)
        
        redo_action = QAction('İleri Al', self)
        redo_action.setShortcut('Ctrl+Y')
        redo_action.triggered.connect(self.scene.redo)
        
        edit_menu.addAction(undo_action)
        edit_menu.addAction(redo_action)

    def setup_shortcuts(self):
        delete_shortcut = QKeySequence(Qt.Key_Delete)
        delete_action = QAction(self)
        delete_action.setShortcut(delete_shortcut)
        delete_action.triggered.connect(self.delete_selected)
        self.addAction(delete_action)

    def set_current_element(self, element):
        self.scene.current_element = element
        self.statusBar().showMessage(f'Seçili element: {element}')

    def set_bond_type(self, bond_type):
        if bond_type == 'Tek Bağ':
            self.scene.current_bond_type = 'single'
        elif bond_type == 'Çift Bağ':
            self.scene.current_bond_type = 'double'
        elif bond_type == 'Üçlü Bağ':
            self.scene.current_bond_type = 'triple'
        self.statusBar().showMessage(f'Seçili bağ tipi: {bond_type}')

    def clear_scene(self):
        reply = QMessageBox.question(self, 'Temizle', 
                                   'Tüm çizimi temizlemek istediğinizden emin misiniz?',
                                   QMessageBox.Yes | QMessageBox.No)
        if reply == QMessageBox.Yes:
            self.scene.clear()
            self.scene.atoms = []
            self.scene.bonds = []
            self.scene.undo_stack.clear()
            self.scene.redo_stack.clear()
            self.statusBar().showMessage('Çizim temizlendi')

    def delete_selected(self):
        if self.scene.selected_item:
            self.scene.removeItem(self.scene.selected_item)
            self.scene.selected_item = None
            self.statusBar().showMessage('Seçili öğe silindi')

    def save_molecule(self):
        filename, selected_filter = QFileDialog.getSaveFileName(
            self, 
            "Molekülü Kaydet",
            "",
            "JSON Dosyası (*.json);;MOL Dosyası (*.mol);;ChemDraw Dosyası (*.cdx)" if HAS_CDX else "JSON Dosyası (*.json);;MOL Dosyası (*.mol)"
        )
        
        if filename:
            try:
                if filename.endswith('.json'):
                    self.save_as_json(filename)
                elif filename.endswith('.mol'):
                    self.save_as_mol(filename)
                elif filename.endswith('.cdx') and HAS_CDX:
                    self.save_as_cdx(filename)
                
                self.statusBar().showMessage(f'Molekül kaydedildi: {filename}')
            except Exception as e:
                QMessageBox.warning(self, "Hata", f"Kaydetme hatası: {str(e)}")

    def save_as_mol(self, filename):
        # Önce RDKit molekülüne dönüştür
        mol = self.create_rdkit_mol()
        if mol:
            # MOL dosyası olarak kaydet
            writer = Chem.SDWriter(filename)
            writer.write(mol)
            writer.close()

    def save_as_cdx(self, filename):
        if not HAS_CDX:
            QMessageBox.warning(self, "Uyarı", "ChemDraw dosya desteği mevcut değil.")
            return
        
        # Önce RDKit molekülüne dönüştür
        mol = self.create_rdkit_mol()
        if mol:
            # CDX dosyası olarak kaydet
            writer = CDXMLWriter.CDXMLWriter(filename)
            writer.write(mol)
            writer.close()

    def create_rdkit_mol(self):
        """RDKit molekülü oluştur"""
        try:
            mol = Chem.RWMol()
            atom_indices = {}
            
            # Atomları ekle
            for i, atom in enumerate(self.scene.atoms):
                rd_atom = Chem.Atom(atom['element'])
                atom_idx = mol.AddAtom(rd_atom)
                atom_indices[atom['id']] = atom_idx
            
            # Bağları ekle
            for bond in self.scene.bonds:
                start_atom = next(a for a in self.scene.atoms if a['position'] == bond['start'])
                end_atom = next(a for a in self.scene.atoms if a['position'] == bond['end'])
                
                start_idx = atom_indices[start_atom['id']]
                end_idx = atom_indices[end_atom['id']]
                
                if bond['type'] == 'single':
                    bond_type = Chem.BondType.SINGLE
                elif bond['type'] == 'double':
                    bond_type = Chem.BondType.DOUBLE
                elif bond['type'] == 'triple':
                    bond_type = Chem.BondType.TRIPLE
                
                mol.AddBond(start_idx, end_idx, bond_type)
            
            return mol.GetMol()
        except Exception as e:
            QMessageBox.warning(self, "Hata", f"RDKit molekülü oluşturma hatası: {str(e)}")
            return None

    def load_molecule(self):
        filename, selected_filter = QFileDialog.getOpenFileName(
            self,
            "Molekül Yükle",
            "",
            "Tüm Desteklenen Formatlar (*.json *.mol *.cdx);;JSON Dosyası (*.json);;MOL Dosyası (*.mol);;ChemDraw Dosyası (*.cdx)" if HAS_CDX else "Tüm Desteklenen Formatlar (*.json *.mol);;JSON Dosyası (*.json);;MOL Dosyası (*.mol)"
        )
        
        if filename:
            try:
                extension = os.path.splitext(filename)[1].lower()
                if extension == '.json':
                    self.load_from_json(filename)
                elif extension == '.mol':
                    self.load_from_mol(filename)
                elif extension == '.cdx' and HAS_CDX:
                    self.load_from_cdx(filename)
                else:
                    QMessageBox.warning(self, "Uyarı", "Desteklenmeyen dosya formatı!")
                    return
                
                self.statusBar().showMessage(f'Molekül yüklendi: {filename}')
            except Exception as e:
                QMessageBox.warning(self, "Hata", f"Yükleme hatası: {str(e)}")

    def load_from_mol(self, filename):
        # MOL dosyasını oku
        mol = Chem.SDMolSupplier(filename)[0]
        if mol is None:
            raise ValueError("Geçersiz MOL dosyası")
        
        self.clear_scene()
        self.molecule_to_scene(mol)

    def load_from_cdx(self, filename):
        if not HAS_CDX:
            QMessageBox.warning(self, "Uyarı", "ChemDraw dosya desteği mevcut değil.")
            return
        
        # CDX dosyasını oku
        mol = CDXMLReader.readCDXML(filename)
        if mol is None:
            raise ValueError("Geçersiz CDX dosyası")
        
        self.clear_scene()
        self.molecule_to_scene(mol)

    def molecule_to_scene(self, mol):
        """RDKit molekülünü sahneye çiz"""
        # Atomları ekle
        conf = mol.GetConformer()
        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            atom_data = {
                'position': QPointF(pos.x * 50 + 350, pos.y * 50 + 250),  # Ekran merkezine göre ölçeklendir
                'element': atom.GetSymbol(),
                'id': len(self.scene.atoms)
            }
            self.scene.atoms.append(atom_data)
            self.scene.draw_atom(atom_data)

        # Bağları ekle
        for bond in mol.GetBonds():
            start_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            
            bond_type = str(bond.GetBondType()).lower()
            if bond_type == 'single':
                bond_type = 'single'
            elif bond_type == 'double':
                bond_type = 'double'
            elif bond_type == 'triple':
                bond_type = 'triple'
            else:
                bond_type = 'single'  # Varsayılan olarak tek bağ
            
            bond_data = {
                'start': self.scene.atoms[start_idx]['position'],
                'end': self.scene.atoms[end_idx]['position'],
                'type': bond_type
            }
            self.scene.bonds.append(bond_data)
            self.scene.draw_bond(bond_data)
    
    def save_as_json(self, filename):
        data = {
            'atoms': [{
                'element': atom['element'],
                'x': atom['position'].x(),
                'y': atom['position'].y()
            } for atom in self.scene.atoms],
            'bonds': [{
                'start_idx': self.scene.atoms.index(next(a for a in self.scene.atoms 
                                                       if a['position'] == bond['start'])),
                'end_idx': self.scene.atoms.index(next(a for a in self.scene.atoms 
                                                     if a['position'] == bond['end'])),
                'type': bond['type']
            } for bond in self.scene.bonds]
        }
        
        with open(filename, 'w') as f:
            json.dump(data, f)

    def load_from_json(self, filename):
        with open(filename, 'r') as f:
            data = json.load(f)
        
        self.clear_scene()
        
        # Atomları yükle
        for atom_data in data['atoms']:
            pos = QPointF(atom_data['x'], atom_data['y'])
            atom = {
                'position': pos,
                'element': atom_data['element'],
                'id': len(self.scene.atoms)
            }
            self.scene.atoms.append(atom)
            self.scene.draw_atom(atom)

        # Bağları yükle
        for bond_data in data['bonds']:
            start_atom = self.scene.atoms[bond_data['start_idx']]
            end_atom = self.scene.atoms[bond_data['end_idx']]
            bond = {
                'start': start_atom['position'],
                'end': end_atom['position'],
                'type': bond_data['type']
            }
            self.scene.bonds.append(bond)
            self.scene.draw_bond(bond)

    def convert_to_rdkit(self):
        if not self.scene.atoms:
            QMessageBox.warning(self, "Uyarı", "Önce bir molekül çizmelisiniz!")
            return

        # RDKit molekülü oluşturma
        mol = Chem.RWMol()
        
        # Atom indekslerini takip etmek için sözlük
        atom_indices = {}
        
        # Atomları ekle
        for i, atom in enumerate(self.scene.atoms):
            rd_atom = Chem.Atom(atom['element'])
            atom_idx = mol.AddAtom(rd_atom)
            atom_indices[atom['id']] = atom_idx
        
        # Bağları ekle
        for bond in self.scene.bonds:
            try:
                start_atom = next(a for a in self.scene.atoms if a['position'] == bond['start'])
                end_atom = next(a for a in self.scene.atoms if a['position'] == bond['end'])
                
                start_idx = atom_indices[start_atom['id']]
                end_idx = atom_indices[end_atom['id']]
                
                if bond['type'] == 'single':
                    bond_type = Chem.BondType.SINGLE
                elif bond['type'] == 'double':
                    bond_type = Chem.BondType.DOUBLE
                elif bond['type'] == 'triple':
                    bond_type = Chem.BondType.TRIPLE
                
                mol.AddBond(start_idx, end_idx, bond_type)
            except Exception as e:
                QMessageBox.warning(self, "Hata", f"Bağ eklenirken hata oluştu: {str(e)}")
                return
        
        # Molekülü tamamla
        try:
            mol = mol.GetMol()
            # SMILES formatına çevir
            smiles = Chem.MolToSmiles(mol)
            QMessageBox.information(self, "SMILES", f"Molekül SMILES: {smiles}")
            self.statusBar().showMessage(f'SMILES oluşturuldu: {smiles}')
        except Exception as e:
            QMessageBox.warning(self, "Hata", f"Molekül dönüştürme hatası: {str(e)}")

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MoleculeEditor()
    window.show()
    sys.exit(app.exec())