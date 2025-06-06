import os
import sqlite3
import json
import threading
import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.widgets import Cursor
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, rdMolDescriptors
import rdkit.RDLogger as RDLogger
from PIL import Image, ImageTk

# Suppress RDKit warnings
RDLogger.DisableLog('rdApp.*')

class EnhancedFTIRPredictor:
    def __init__(self):
        self.model = None
        self.wavenumbers = np.linspace(400, 4000, 1000)
        self.data_directory = "ftir_data"
        self.database_file = os.path.join(self.data_directory, "ftir_database.json")
        self.initialize_database()
        self.ir_table = self.load_ir_table()
        
    def load_ir_table(self):
        """IR Spectrum Table.pdf'den alınan bilgilerle oluşturulmuş genişletilmiş IR tablosu"""
        return {
            # O-H stretching
            'alcohol': [
                {'range': (3700, 3584), 'center': 3650, 'width': 50, 'height': 0.6, 'appearance': 'medium, sharp', 'type': 'free'},
                {'range': (3550, 3200), 'center': 3400, 'width': 150, 'height': 0.8, 'appearance': 'strong, broad', 'type': 'intermolecular bonded'},
                {'range': (3200, 2700), 'center': 3000, 'width': 200, 'height': 0.4, 'appearance': 'weak, broad', 'type': 'intramolecular bonded'}
            ],
            # N-H stretching
            'amine': [
                {'range': (3400, 3300), 'center': 3350, 'width': 80, 'height': 0.5, 'appearance': 'medium', 'type': 'primary amine'},
                {'range': (3330, 3250), 'center': 3300, 'width': 60, 'height': 0.5, 'appearance': 'medium', 'type': 'secondary amine'},
                {'range': (3000, 2800), 'center': 2900, 'width': 100, 'height': 0.7, 'appearance': 'strong, broad', 'type': 'amine salt'}
            ],
            # C-H stretching
            'alkyne': [
                {'range': (3333, 3267), 'center': 3300, 'width': 30, 'height': 0.7, 'appearance': 'strong, sharp', 'type': 'terminal alkyne'}
            ],
            'alkene': [
                {'range': (3100, 3000), 'center': 3050, 'width': 40, 'height': 0.5, 'appearance': 'medium', 'type': 'vinyl C-H'}
            ],
            'alkane': [
                {'range': (3000, 2840), 'center': 2920, 'width': 60, 'height': 0.5, 'appearance': 'medium', 'type': 'sp3 C-H'}
            ],
            'aldehyde': [
                {'range': (2830, 2695), 'center': 2720, 'width': 60, 'height': 0.4, 'appearance': 'medium', 'type': 'aldehyde C-H (doublet)'}
            ],
            # C=O stretching
            'carboxylic_acid': [
                {'range': (1720, 1706), 'center': 1715, 'width': 25, 'height': 0.9, 'appearance': 'strong', 'type': 'dimer'},
                {'range': (1760, 1760), 'center': 1760, 'width': 20, 'height': 0.8, 'appearance': 'strong', 'type': 'monomer'}
            ],
            'ester': [
                {'range': (1750, 1735), 'center': 1740, 'width': 25, 'height': 0.9, 'appearance': 'strong', 'type': 'normal ester'},
                {'range': (1770, 1780), 'center': 1775, 'width': 25, 'height': 0.9, 'appearance': 'strong', 'type': 'vinyl/phenyl ester'}
            ],
            'aldehyde': [
                {'range': (1740, 1720), 'center': 1730, 'width': 20, 'height': 0.9, 'appearance': 'strong', 'type': 'normal aldehyde'},
                {'range': (1710, 1685), 'center': 1700, 'width': 20, 'height': 0.9, 'appearance': 'strong', 'type': 'conjugated aldehyde'}
            ],
            'ketone': [
                {'range': (1725, 1705), 'center': 1715, 'width': 25, 'height': 0.9, 'appearance': 'strong', 'type': 'aliphatic ketone'},
                {'range': (1685, 1666), 'center': 1675, 'width': 25, 'height': 0.9, 'appearance': 'strong', 'type': 'conjugated ketone'}
            ],
            'amide': [
                {'range': (1690, 1690), 'center': 1690, 'width': 30, 'height': 0.8, 'appearance': 'strong', 'type': 'primary amide (free)'},
                {'range': (1680, 1680), 'center': 1680, 'width': 30, 'height': 0.8, 'appearance': 'strong', 'type': 'secondary amide (free)'}
            ],
            # C=C stretching
            'alkene': [
                {'range': (1678, 1668), 'center': 1675, 'width': 30, 'height': 0.5, 'appearance': 'weak', 'type': 'trans disubstituted'},
                {'range': (1662, 1626), 'center': 1640, 'width': 30, 'height': 0.6, 'appearance': 'medium', 'type': 'cis disubstituted'},
                {'range': (1648, 1638), 'center': 1645, 'width': 30, 'height': 0.7, 'appearance': 'strong', 'type': 'monosubstituted'}
            ],
            # C≡N stretching
            'nitrile': [
                {'range': (2260, 2222), 'center': 2240, 'width': 30, 'height': 0.5, 'appearance': 'weak', 'type': 'nitrile'}
            ],
            # Other important groups
            'thiol': [
                {'range': (2600, 2550), 'center': 2575, 'width': 40, 'height': 0.3, 'appearance': 'weak', 'type': 'S-H stretch'}
            ],
            'ether': [
                {'range': (1150, 1085), 'center': 1120, 'width': 50, 'height': 0.6, 'appearance': 'strong', 'type': 'aliphatic ether C-O'}
            ],
            'aromatic': [
                {'range': (1600, 1585), 'center': 1590, 'width': 20, 'height': 0.5, 'appearance': 'medium', 'type': 'aromatic C=C'},
                {'range': (1500, 1500), 'center': 1500, 'width': 20, 'height': 0.5, 'appearance': 'medium', 'type': 'aromatic C=C'},
                {'range': (3030, 3030), 'center': 3030, 'width': 30, 'height': 0.3, 'appearance': 'weak', 'type': 'aromatic C-H'}
            ]
        }

    def initialize_database(self):
        if not os.path.exists(self.data_directory):
            os.makedirs(self.data_directory)
        if not os.path.exists(self.database_file):
            self.create_empty_database()
            
    def create_empty_database(self):
        empty_db = {
            "compounds": {},
            "version": "1.0"
        }
        with open(self.database_file, 'w') as f:
            json.dump(empty_db, f, indent=4)

    def identify_functional_groups(self, mol):
        """Moleküldeki fonksiyonel grupları tespit et - Genişletilmiş versiyon"""
        functional_groups = {
            'alcohol': False,
            'carboxylic_acid': False,
            'ester': False,
            'aldehyde': False,
            'ketone': False,
            'ether': False,
            'alkene': False,
            'alkyne': False,
            'aromatic': False,
            'amine': False,
            'amide': False,
            'nitrile': False,
            'thiol': False
        }
        
        # SMARTS patterns for functional groups (genişletilmiş)
        patterns = {
            'alcohol': '[OH1][CX4]',
            'carboxylic_acid': '[CX3](=O)[OH]',
            'ester': '[CX3](=O)[OX2H0]',
            'aldehyde': '[CX3H1](=O)[#6]',
            'ketone': '[CX3](=O)[#6]',
            'ether': '[OX2]([#6])[#6]',
            'alkene': '[CX3]=[CX3]',
            'alkyne': '[CX2]#[CX2]',
            'aromatic': 'c1ccccc1',
            'amine': '[NX3;H2,H1;!$(NC=O)]',
            'amide': '[NX3][CX3](=[OX1])[#6]',
            'nitrile': '[NX1]#[CX2]',
            'thiol': '[SH]'
        }
        
        # Check each pattern
        for group, pattern in patterns.items():
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                functional_groups[group] = True
                
        return functional_groups

    def generate_realistic_spectrum(self, mol):
        wavenumbers = np.linspace(400, 4000, 1000)
        spectrum = np.zeros_like(wavenumbers)
        
        # Fonksiyonel grupları tespit et
        groups = self.identify_functional_groups(mol)
        
        # IR tablosundaki bilgilere göre pikler oluştur
        for group, present in groups.items():
            if present and group in self.ir_table:
                for peak_info in self.ir_table[group]:
                    spectrum += self.add_peak(
                        wavenumbers,
                        peak_info['center'],
                        peak_info['width'],
                        peak_info['height']
                    )
        
        # Aliphatic C-H stretching (if any carbon present)
        if any(atom.GetSymbol() == 'C' for atom in mol.GetAtoms()):
            spectrum += self.add_peak(wavenumbers, 2960, 30, 0.3)  # CH3 asymmetric
            spectrum += self.add_peak(wavenumbers, 2930, 30, 0.3)  # CH2 asymmetric
            spectrum += self.add_peak(wavenumbers, 2870, 30, 0.2)  # CH3 symmetric
            
        # Add minimal noise
        spectrum += np.random.normal(0, 0.01, len(wavenumbers))
        
        return {
            'wavenumbers': wavenumbers.tolist(),
            'absorbance': spectrum.tolist(),
            'functional_groups': groups
        }
        
    def add_peak(self, x, center, width, height):
        """Gaussian peak oluştur"""
        return height * np.exp(-(x - center)**2 / (2 * width**2))
        
    def predict_ftir(self, mol_input, input_format='mol'):
        try:
            mol = self.read_molecule(mol_input, input_format)
            if mol is None:
                return None
                
            spectrum_data = self.generate_realistic_spectrum(mol)
            
            if spectrum_data is None:
                return None
                
            # DataFrame'i doğru column isimleriyle oluştur
            df = pd.DataFrame({
                'Wavenumber': spectrum_data['wavenumbers'],
                'Absorbance': spectrum_data['absorbance']
            })
            
            # Fonksiyonel grup bilgilerini de ekleyelim
            df.attrs['functional_groups'] = spectrum_data['functional_groups']
            df.attrs['peak_info'] = self.get_peak_assignments(spectrum_data['functional_groups'])
                
            return df
                
        except Exception as e:
            print(f"Prediction error: {str(e)}")
            return None
            
    def get_peak_assignments(self, functional_groups):
        """Tespit edilen fonksiyonel gruplara göre pik atamalarını getir"""
        peak_assignments = []
        
        for group, present in functional_groups.items():
            if present and group in self.ir_table:
                for peak in self.ir_table[group]:
                    peak_assignments.append({
                        'group': group,
                        'center': peak['center'],
                        'range': peak['range'],
                        'appearance': peak['appearance'],
                        'type': peak['type']
                    })
        
        return peak_assignments
            
    def read_molecule(self, input_string, input_format='mol'):
        try:
            if input_format.lower() in ['mol', 'sdf']:
                if os.path.exists(input_string):
                    mol = Chem.MolFromMolFile(input_string)
                else:
                    return None
            else:
                raise ValueError(f"Unsupported format: {input_format}")
                
            if mol is None:
                raise ValueError("Failed to parse molecule")
                
            return mol
        except Exception as e:
            print(f"Molecule reading error: {str(e)}")
            return None

class FTIRDatabase:
    def __init__(self):
        self.db_path = os.path.join(os.path.dirname(__file__), "ftir_database.sqlite")
        self.initialize_database()

    def initialize_database(self):
        # Veritabanı bağlantısı oluştur
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        # Fonksiyonel gruplar tablosu
        cursor.execute('''
        CREATE TABLE IF NOT EXISTS functional_groups (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT NOT NULL,
            wave_min FLOAT NOT NULL,
            wave_max FLOAT NOT NULL,
            intensity TEXT,
            description TEXT,
            compound_type TEXT
        )
        ''')

        # Örnek veriler
        default_data = [
            # Karbonil grupları (C=O)
            ("Aldehit C=O", 1740, 1720, "Strong", "Aldehyde carbonyl stretch", "aldehyde"),
            ("Keton C=O", 1720, 1708, "Strong", "Ketone carbonyl stretch", "ketone"),
            ("Karboksilik asit C=O", 1725, 1700, "Strong", "Carboxylic acid carbonyl stretch", "carboxylic_acid"),
            ("Amit C=O", 1680, 1630, "Strong", "Amide carbonyl stretch", "amide"),
            ("Ester C=O", 1750, 1730, "Strong", "Ester carbonyl stretch", "ester"),
            
            # O-H grupları
            ("Alkol O-H", 3400, 3200, "Strong, broad", "Alcohol O-H stretch", "alcohol"),
            ("Karboksilik asit O-H", 3000, 2500, "Strong, very broad", "Carboxylic acid O-H stretch", "carboxylic_acid"),
            ("Fenol O-H", 3550, 3200, "Strong, broad", "Phenol O-H stretch", "phenol"),
            
            # N-H grupları
            ("Primer amin N-H", 3500, 3300, "Medium", "Primary amine N-H stretch", "primary_amine"),
            ("Sekonder amin N-H", 3350, 3310, "Medium", "Secondary amine N-H stretch", "secondary_amine"),
            ("Amit N-H", 3350, 3180, "Medium", "Amide N-H stretch", "amide"),
            
            # C-H grupları
            ("Alkan C-H", 2960, 2850, "Medium", "Alkane C-H stretch", "alkane"),
            ("Alken C-H", 3100, 3000, "Medium", "Alkene C-H stretch", "alkene"),
            ("Alkin C-H", 3300, 3270, "Strong, sharp", "Alkyne C-H stretch", "alkyne"),
            ("Aromatik C-H", 3100, 3000, "Weak", "Aromatic C-H stretch", "aromatic"),
            
            # C-O grupları
            ("Alkol C-O", 1260, 1050, "Strong", "Alcohol C-O stretch", "alcohol"),
            ("Eter C-O", 1150, 1070, "Strong", "Ether C-O stretch", "ether"),
            ("Ester C-O", 1300, 1000, "Two bands or more", "Ester C-O stretch", "ester"),
            
            # Diğer gruplar
            ("C=C", 1680, 1620, "Variable", "Alkene C=C stretch", "alkene"),
            ("C≡C", 2260, 2100, "Variable", "Alkyne C≡C stretch", "alkyne"),
            ("C≡N", 2260, 2220, "Medium", "Nitrile C≡N stretch", "nitrile")
        ]

        # Verileri ekle
        cursor.execute("SELECT COUNT(*) FROM functional_groups")
        if cursor.fetchone()[0] == 0:  # Tablo boşsa
            cursor.executemany('''
            INSERT INTO functional_groups 
            (name, wave_max, wave_min, intensity, description, compound_type)
            VALUES (?, ?, ?, ?, ?, ?)
            ''', default_data)

        conn.commit()
        conn.close()

    def get_functional_groups(self, wavenumber):
        """Belirli bir dalga sayısına karşılık gelen fonksiyonel grupları bul"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        cursor.execute('''
        SELECT name, description, intensity, compound_type
        FROM functional_groups
        WHERE ? BETWEEN wave_min AND wave_max
        ''', (wavenumber,))
        
        results = cursor.fetchall()
        conn.close()
        return results

    def get_compound_info(self, compound_type):
        """Belirli bir bileşik türüne ait tüm fonksiyonel grupları getir"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        cursor.execute('''
        SELECT name, wave_min, wave_max, intensity, description
        FROM functional_groups
        WHERE compound_type = ?
        ''', (compound_type,))
        
        results = cursor.fetchall()
        conn.close()
        return results

    def add_functional_group(self, name, wave_min, wave_max, intensity, description, compound_type):
        """Yeni fonksiyonel grup ekle"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        cursor.execute('''
        INSERT INTO functional_groups 
        (name, wave_min, wave_max, intensity, description, compound_type)
        VALUES (?, ?, ?, ?, ?, ?)
        ''', (name, wave_min, wave_max, intensity, description, compound_type))
        
        conn.commit()
        conn.close()

class FTIRPredictor:
    def __init__(self):
        self.model = None
        self.wavenumbers = np.linspace(400, 4000, 1000)
        self.data_directory = "ftir_data"
        self.database_file = os.path.join(self.data_directory, "ftir_database.json")
        self.initialize_database()
        
    def initialize_database(self):
        if not os.path.exists(self.data_directory):
            os.makedirs(self.data_directory)
        if not os.path.exists(self.database_file):
            self.create_empty_database()
            
    def create_empty_database(self):
        empty_db = {
            "compounds": {},
            "version": "1.0"
        }
        with open(self.database_file, 'w') as f:
            json.dump(empty_db, f, indent=4)

    def identify_functional_groups(self, mol):
        """Moleküldeki fonksiyonel grupları tespit et"""
        functional_groups = {
            'alcohol': False,
            'carboxylic_acid': False,
            'ester': False,
            'aldehyde': False,
            'ketone': False,
            'ether': False,
            'alkene': False,
            'alkyne': False,
            'aromatic': False,
            'amine': False,
            'amide': False
        }
        
        # SMARTS patterns for functional groups
        patterns = {
            'alcohol': '[OH1][CX4]',
            'carboxylic_acid': '[CX3](=O)[OH]',
            'ester': '[CX3](=O)[OX2H0]',
            'aldehyde': '[CX3H1](=O)[#6]',
            'ketone': '[CX3](=O)[#6]',
            'ether': '[OX2]([#6])[#6]',
            'alkene': '[CX3]=[CX3]',
            'alkyne': '[CX2]#[CX2]',
            'aromatic': 'c1ccccc1',
            'amine': '[NX3;H2,H1;!$(NC=O)]',
            'amide': '[NX3][CX3](=[OX1])[#6]'
        }
        
        # Check each pattern
        for group, pattern in patterns.items():
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                functional_groups[group] = True
                
        return functional_groups

    def generate_realistic_spectrum(self, mol):
        wavenumbers = np.linspace(400, 4000, 1000)
        spectrum = np.zeros_like(wavenumbers)
        
        # Fonksiyonel grupları tespit et
        groups = self.identify_functional_groups(mol)
        
        # Her fonksiyonel grup için karakteristik pikleri ekle
        if groups['alcohol']:
            spectrum += self.add_peak(wavenumbers, 3350, 100, 0.6)  # O-H stretch
            spectrum += self.add_peak(wavenumbers, 1050, 40, 0.4)   # C-O stretch
            
        if groups['carboxylic_acid']:
            spectrum += self.add_peak(wavenumbers, 3000, 300, 0.5)  # O-H stretch (broad)
            spectrum += self.add_peak(wavenumbers, 1720, 30, 0.8)   # C=O stretch
            spectrum += self.add_peak(wavenumbers, 1250, 50, 0.4)   # C-O stretch
            
        if groups['ester']:
            spectrum += self.add_peak(wavenumbers, 1735, 25, 0.8)   # C=O stretch
            spectrum += self.add_peak(wavenumbers, 1250, 40, 0.4)   # C-O stretch
            spectrum += self.add_peak(wavenumbers, 1050, 40, 0.3)   # C-O stretch
            
        if groups['aldehyde']:
            spectrum += self.add_peak(wavenumbers, 1730, 20, 0.8)   # C=O stretch
            spectrum += self.add_peak(wavenumbers, 2720, 20, 0.3)   # C-H aldehyde
            
        if groups['ketone']:
            spectrum += self.add_peak(wavenumbers, 1715, 25, 0.8)   # C=O stretch
            
        if groups['ether']:
            spectrum += self.add_peak(wavenumbers, 1100, 50, 0.4)   # C-O stretch
            
        if groups['alkene']:
            spectrum += self.add_peak(wavenumbers, 1640, 30, 0.4)   # C=C stretch
            spectrum += self.add_peak(wavenumbers, 3080, 30, 0.3)   # =C-H stretch
            
        if groups['alkyne']:
            spectrum += self.add_peak(wavenumbers, 2120, 30, 0.3)   # C≡C stretch
            spectrum += self.add_peak(wavenumbers, 3300, 20, 0.4)   # ≡C-H stretch
            
        if groups['aromatic']:
            spectrum += self.add_peak(wavenumbers, 3030, 30, 0.2)   # Aromatic C-H
            spectrum += self.add_peak(wavenumbers, 1600, 20, 0.3)   # C=C aromatic
            spectrum += self.add_peak(wavenumbers, 1500, 20, 0.3)   # C=C aromatic
            spectrum += self.add_peak(wavenumbers, 750, 20, 0.4)    # Out-of-plane bend
            
        if groups['amine']:
            spectrum += self.add_peak(wavenumbers, 3400, 80, 0.4)   # N-H stretch
            spectrum += self.add_peak(wavenumbers, 1600, 30, 0.3)   # N-H bend
            
        if groups['amide']:
            spectrum += self.add_peak(wavenumbers, 3300, 80, 0.5)   # N-H stretch
            spectrum += self.add_peak(wavenumbers, 1650, 30, 0.8)   # C=O stretch
            spectrum += self.add_peak(wavenumbers, 1550, 30, 0.4)   # N-H bend

        # Aliphatic C-H stretching (if any carbon present)
        if any(atom.GetSymbol() == 'C' for atom in mol.GetAtoms()):
            spectrum += self.add_peak(wavenumbers, 2960, 30, 0.3)  # CH3 asymmetric
            spectrum += self.add_peak(wavenumbers, 2930, 30, 0.3)  # CH2 asymmetric
            spectrum += self.add_peak(wavenumbers, 2870, 30, 0.2)  # CH3 symmetric
            
        # Add minimal noise
        spectrum += np.random.normal(0, 0.01, len(wavenumbers))
        
        return {
            'wavenumbers': wavenumbers.tolist(),
            'absorbance': spectrum.tolist()
        }
        
    def add_peak(self, x, center, width, height):
        """Gaussian peak oluştur"""
        return height * np.exp(-(x - center)**2 / (2 * width**2))
        
    def predict_ftir(self, mol_input, input_format='mol'):
        try:
            mol = self.read_molecule(mol_input, input_format)
            if mol is None:
                return None
                
            spectrum_data = self.generate_realistic_spectrum(mol)
            
            if spectrum_data is None:
                return None
                
            # DataFrame'i doğru column isimleriyle oluştur
            return pd.DataFrame({
                'Wavenumber': spectrum_data['wavenumbers'],
                'Absorbance': spectrum_data['absorbance']
            })
                
        except Exception as e:
            print(f"Prediction error: {str(e)}")
            return None
            
    def read_molecule(self, input_string, input_format='mol'):
        try:
            if input_format.lower() in ['mol', 'sdf']:
                if os.path.exists(input_string):
                    mol = Chem.MolFromMolFile(input_string)
                else:
                    return None
            else:
                raise ValueError(f"Unsupported format: {input_format}")
                
            if mol is None:
                raise ValueError("Failed to parse molecule")
                
            return mol
        except Exception as e:
            print(f"Molecule reading error: {str(e)}")
            return None

class ModernFTIRGUI:
    def __init__(self, root):
        self.root = root
        self.setup_window()
        self.setup_colors_and_fonts()
        
        # Initialize variables
        self.current_molecule_groups = None
        self.peak_details = {}
        self.db = FTIRDatabase()
        self.predictor = FTIRPredictor()
        self.current_mol = None
        self.highlighted_atoms = set()
        self.highlighted_bonds = set()
        self.atom_to_group = {}
        self.loading = False
        
        # Peak assignments dictionary
        self.peak_assignments = {
            (3200, 3600): "O-H stretch\nAlcohols and carboxylic acids",
            (1670, 1780): "C=O stretch\nAldehydes, ketones, carboxylic acids",
            (2850, 3000): "C-H stretch\nMost organic compounds",
            (3300, 3500): "N-H stretch\nAmines and amides",
            (1030, 1230): "C-N stretch\nAmines and amides",
            (1350, 1480): "C-H bending\nMethyl and methylene groups",
            (500, 1500): "Fingerprint region\nUnique molecular pattern"
        }
        
        self.create_modern_layout()

    def setup_window(self):
        """Configure main window with modern styling"""
        self.root.title("FTIR Spectrum Predictor - Advanced Analysis Suite")
        self.root.geometry("1600x900")
        self.root.minsize(1200, 700)
        
        # Configure window icon and styling
        self.root.configure(bg='#1a1a2e')
        
        # Center window on screen
        self.center_window()

    def center_window(self):
        """Center the window on screen"""
        self.root.update_idletasks()
        width = self.root.winfo_width()
        height = self.root.winfo_height()
        x = (self.root.winfo_screenwidth() // 2) - (width // 2)
        y = (self.root.winfo_screenheight() // 2) - (height // 2)
        self.root.geometry(f'{width}x{height}+{x}+{y}')

    def setup_colors_and_fonts(self):
        """Setup modern color scheme and fonts"""
        self.colors = {
            'primary': '#4f46e5',      # Indigo
            'primary_dark': '#3730a3',
            'secondary': '#06b6d4',     # Cyan
            'accent': '#f59e0b',        # Amber
            'success': '#10b981',       # Emerald
            'danger': '#ef4444',        # Red
            'warning': '#f59e0b',       # Amber
            'bg_primary': '#1a1a2e',    # Dark blue
            'bg_secondary': '#16213e',  # Darker blue
            'bg_tertiary': '#eee2dc',   # Light cream
            'text_primary': '#ffffff',
            'text_secondary': '#94a3b8',
            'border': '#334155',
            'card_bg': '#0f172a',
            'glass_bg': 'rgba(255, 255, 255, 0.1)'
        }
        
        # Configure modern fonts
        self.fonts = {
            'title': ('Inter', 24, 'bold'),
            'subtitle': ('Inter', 14, 'normal'),
            'heading': ('Inter', 16, 'bold'),
            'body': ('Inter', 11, 'normal'),
            'small': ('Inter', 9, 'normal'),
            'button': ('Inter', 10, 'bold')
        }

    def create_modern_layout(self):
        """Create modern layout with cards and animations"""
        # Configure ttk styles
        self.configure_modern_styles()
        
        # Main container with gradient-like effect
        main_frame = tk.Frame(self.root, bg=self.colors['bg_primary'])
        main_frame.pack(fill='both', expand=True, padx=20, pady=20)
        
        # Header section
        self.create_header(main_frame)
        
        # Content area
        content_frame = tk.Frame(main_frame, bg=self.colors['bg_primary'])
        content_frame.pack(fill='both', expand=True, pady=(20, 0))
        
        # Create three-column layout
        self.create_three_column_layout(content_frame)

    def configure_modern_styles(self):
        """Configure modern ttk styles"""
        style = ttk.Style()
        
        # Configure modern button style
        style.configure('Modern.TButton',
                       background=self.colors['primary'],
                       foreground='white',
                       borderwidth=0,
                       focuscolor='none',
                       padding=(20, 10),
                       font=self.fonts['button'])
        
        style.map('Modern.TButton',
                 background=[('active', self.colors['primary_dark']),
                           ('pressed', self.colors['primary_dark'])])
        
        # Modern entry style
        style.configure('Modern.TEntry',
                       fieldbackground=self.colors['bg_secondary'],
                       foreground=self.colors['text_primary'],
                       borderwidth=1,
                       insertcolor=self.colors['text_primary'])
        
        # Modern label style
        style.configure('Modern.TLabel',
                       background=self.colors['bg_primary'],
                       foreground=self.colors['text_primary'],
                       font=self.fonts['body'])
        
        # Modern frame style
        style.configure('Card.TFrame',
                       background=self.colors['bg_secondary'],
                       relief='flat',
                       borderwidth=1)

    def create_header(self, parent):
        """Create modern header with title and navigation"""
        header_frame = tk.Frame(parent, bg=self.colors['bg_primary'], height=100)
        header_frame.pack(fill='x', pady=(0, 20))
        header_frame.pack_propagate(False)
        
        # Title section
        title_frame = tk.Frame(header_frame, bg=self.colors['bg_primary'])
        title_frame.pack(side='left', fill='y')
        
        title_label = tk.Label(title_frame, 
                              text="FTIR Spectrum Predictor",
                              font=self.fonts['title'],
                              fg=self.colors['text_primary'],
                              bg=self.colors['bg_primary'])
        title_label.pack(anchor='w', pady=(10, 0))
        
        subtitle_label = tk.Label(title_frame,
                                 text="Advanced Molecular Analysis & Spectral Prediction",
                                 font=self.fonts['subtitle'],
                                 fg=self.colors['text_secondary'],
                                 bg=self.colors['bg_primary'])
        subtitle_label.pack(anchor='w')
        
        # Status indicator
        status_frame = tk.Frame(header_frame, bg=self.colors['bg_primary'])
        status_frame.pack(side='right', fill='y', padx=20)
        
        self.status_indicator = tk.Label(status_frame,
                                        text="● Ready",
                                        font=self.fonts['body'],
                                        fg=self.colors['success'],
                                        bg=self.colors['bg_primary'])
        self.status_indicator.pack(anchor='e', pady=20)

    def create_three_column_layout(self, parent):
        """Create modern three-column layout"""
        # Left panel - File input and molecule info
        left_panel = self.create_card_frame(parent)
        left_panel.pack(side='left', fill='both', expand=False, padx=(0, 10))
        left_panel.configure(width=400)
        
        # Middle panel - Molecule viewer
        middle_panel = self.create_card_frame(parent)
        middle_panel.pack(side='left', fill='both', expand=True, padx=(0, 10))
        
        # Right panel - Spectrum viewer
        right_panel = self.create_card_frame(parent)
        right_panel.pack(side='right', fill='both', expand=True)
        
        self.setup_left_panel(left_panel)
        self.setup_middle_panel(middle_panel)
        self.setup_right_panel(right_panel)

    def create_card_frame(self, parent):
        """Create a modern card-style frame"""
        card = tk.Frame(parent, 
                       bg=self.colors['bg_secondary'],
                       relief='solid',
                       bd=1,
                       highlightbackground=self.colors['border'],
                       highlightthickness=1)
        return card

    def setup_left_panel(self, panel):
        """Setup left panel with file input and controls"""
        # Panel header
        header = tk.Frame(panel, bg=self.colors['bg_secondary'], height=60)
        header.pack(fill='x', padx=20, pady=(20, 0))
        header.pack_propagate(False)
        
        header_label = tk.Label(header,
                               text="Molecule Input",
                               font=self.fonts['heading'],
                               fg=self.colors['text_primary'],
                               bg=self.colors['bg_secondary'])
        header_label.pack(anchor='w', pady=10)
        
        # File selection section
        file_section = tk.Frame(panel, bg=self.colors['bg_secondary'])
        file_section.pack(fill='x', padx=20, pady=10)
        
        file_label = tk.Label(file_section,
                             text="MOL File:",
                             font=self.fonts['body'],
                             fg=self.colors['text_secondary'],
                             bg=self.colors['bg_secondary'])
        file_label.pack(anchor='w', pady=(0, 5))
        
        # File input frame
        file_input_frame = tk.Frame(file_section, bg=self.colors['bg_secondary'])
        file_input_frame.pack(fill='x', pady=(0, 10))
        
        self.file_path = tk.StringVar()
        file_entry = tk.Entry(file_input_frame,
                             textvariable=self.file_path,
                             font=self.fonts['body'],
                             bg=self.colors['card_bg'],
                             fg=self.colors['text_primary'],
                             insertbackground=self.colors['text_primary'],
                             relief='flat',
                             bd=10)
        file_entry.pack(side='left', fill='x', expand=True, ipady=8)
        
        browse_btn = self.create_modern_button(file_input_frame, "Browse", self.browse_file, self.colors['secondary'])
        browse_btn.pack(side='right', padx=(10, 0))
        
        # Predict button
        predict_btn = self.create_modern_button(file_section, "🔬 Predict Spectrum", self.predict_spectrum, self.colors['primary'])
        predict_btn.pack(fill='x', pady=10)
        
        # Loading indicator
        self.loading_frame = tk.Frame(file_section, bg=self.colors['bg_secondary'])
        
        self.loading_label = tk.Label(self.loading_frame,
                                     text="⟳ Analyzing molecule...",
                                     font=self.fonts['body'],
                                     fg=self.colors['accent'],
                                     bg=self.colors['bg_secondary'])
        self.loading_label.pack()
        
        # Molecule info section
        info_section = tk.Frame(panel, bg=self.colors['bg_secondary'])
        info_section.pack(fill='x', padx=20, pady=20)
        
        info_header = tk.Label(info_section,
                              text="Molecule Information",
                              font=self.fonts['heading'],
                              fg=self.colors['text_primary'],
                              bg=self.colors['bg_secondary'])
        info_header.pack(anchor='w', pady=(0, 10))
        
        # Info display
        self.mol_info_frame = tk.Frame(info_section, 
                                      bg=self.colors['card_bg'],
                                      relief='solid',
                                      bd=1)
        self.mol_info_frame.pack(fill='x', pady=5)
        
        self.mol_info = tk.Label(self.mol_info_frame,
                                text="Load a MOL file to view molecular properties",
                                font=self.fonts['small'],
                                fg=self.colors['text_secondary'],
                                bg=self.colors['card_bg'],
                                wraplength=350,
                                justify='left')
        self.mol_info.pack(padx=15, pady=15)
        
        # Functional groups section
        groups_section = tk.Frame(panel, bg=self.colors['bg_secondary'])
        groups_section.pack(fill='both', expand=True, padx=20, pady=(0, 20))
        
        groups_header = tk.Label(groups_section,
                                text="Detected Functional Groups",
                                font=self.fonts['heading'],
                                fg=self.colors['text_primary'],
                                bg=self.colors['bg_secondary'])
        groups_header.pack(anchor='w', pady=(0, 10))
        
        # Groups display
        self.groups_frame = tk.Frame(groups_section, bg=self.colors['bg_secondary'])
        self.groups_frame.pack(fill='both', expand=True)

    def setup_middle_panel(self, panel):
        """Setup middle panel with molecule viewer"""
        # Panel header
        header = tk.Frame(panel, bg=self.colors['bg_secondary'], height=60)
        header.pack(fill='x', padx=20, pady=(20, 0))
        header.pack_propagate(False)
        
        header_label = tk.Label(header,
                               text="Molecule Structure",
                               font=self.fonts['heading'],
                               fg=self.colors['text_primary'],
                               bg=self.colors['bg_secondary'])
        header_label.pack(anchor='w', pady=10)
        
        # Molecule viewer
        viewer_frame = tk.Frame(panel, bg=self.colors['card_bg'])
        viewer_frame.pack(fill='both', expand=True, padx=20, pady=(10, 20))
        
        # Create matplotlib figure for molecule
        plt.style.use('dark_background')
        self.mol_fig, self.mol_ax = plt.subplots(figsize=(6, 6), 
                                                facecolor=self.colors['card_bg'])
        self.mol_ax.set_facecolor(self.colors['card_bg'])
        self.mol_canvas = FigureCanvasTkAgg(self.mol_fig, master=viewer_frame)
        self.mol_canvas.draw()
        self.mol_canvas.get_tk_widget().pack(fill='both', expand=True)
        self.mol_ax.axis('off')
        
        # Add placeholder text
        self.mol_ax.text(0.5, 0.5, 'Load a MOL file to view\nmolecular structure',
                        ha='center', va='center',
                        transform=self.mol_ax.transAxes,
                        fontsize=14,
                        color=self.colors['text_secondary'],
                        alpha=0.7)

    def setup_right_panel(self, panel):
        """Setup right panel with spectrum viewer"""
        # Panel header
        header = tk.Frame(panel, bg=self.colors['bg_secondary'], height=60)
        header.pack(fill='x', padx=20, pady=(20, 0))
        header.pack_propagate(False)
        
        header_label = tk.Label(header,
                               text="FTIR Spectrum",
                               font=self.fonts['heading'],
                               fg=self.colors['text_primary'],
                               bg=self.colors['bg_secondary'])
        header_label.pack(side='left', pady=10)
        
        # Spectrum controls
        controls_frame = tk.Frame(header, bg=self.colors['bg_secondary'])
        controls_frame.pack(side='right', pady=10)
        
        export_btn = self.create_small_button(controls_frame, "Export", self.export_spectrum)
        export_btn.pack(side='right', padx=(5, 0))
        
        # Spectrum viewer
        spectrum_frame = tk.Frame(panel, bg=self.colors['card_bg'])
        spectrum_frame.pack(fill='both', expand=True, padx=20, pady=(10, 20))
        
        # Create matplotlib figure for spectrum
        self.fig, self.ax = plt.subplots(figsize=(10, 6), 
                                        facecolor=self.colors['card_bg'])
        self.ax.set_facecolor(self.colors['card_bg'])
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=spectrum_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill='both', expand=True)
        
        # Initial plot settings with modern styling
        self.ax.set_xlabel('Wavenumber (cm⁻¹)', fontsize=12, color=self.colors['text_primary'])
        self.ax.set_ylabel('Transmittance (%)', fontsize=12, color=self.colors['text_primary'])
        self.ax.grid(True, linestyle='--', alpha=0.2, color=self.colors['text_secondary'])
        self.ax.set_xlim(4000, 400)
        self.ax.set_ylim(0, 100)
        
        # Style the plot
        self.ax.spines['bottom'].set_color(self.colors['text_secondary'])
        self.ax.spines['left'].set_color(self.colors['text_secondary'])
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['right'].set_visible(False)
        self.ax.tick_params(colors=self.colors['text_secondary'])
        
        # Add placeholder text
        self.ax.text(0.5, 0.5, 'Load a molecule to view\npredicted FTIR spectrum',
                    ha='center', va='center',
                    transform=self.ax.transAxes,
                    fontsize=14,
                    color=self.colors['text_secondary'],
                    alpha=0.7)
        
        # Add cursor and annotations
        self.cursor = Cursor(self.ax, useblit=True, color=self.colors['accent'], linewidth=1)
        self.setup_spectrum_annotations()

    def create_small_button(self, parent, text, command):
        """Create a small modern button"""
        btn = tk.Button(parent,
                       text=text,
                       command=command,
                       font=self.fonts['small'],
                       bg=self.colors['secondary'],
                       fg='white',
                       relief='flat',
                       bd=0,
                       padx=10,
                       pady=5,
                       cursor='hand2')
        return btn

    def setup_spectrum_annotations(self):
        """Setup spectrum annotations"""
        self.annot = self.ax.annotate("", 
                                     xy=(0, 0), 
                                     xytext=(15, 15),
                                     textcoords="offset points",
                                     bbox=dict(boxstyle='round,pad=0.4',
                                             fc=self.colors['accent'],
                                             ec=self.colors['primary'],
                                             lw=2,
                                             alpha=0.95),
                                     fontsize=10,
                                     color='white',
                                     fontweight='bold')
        self.annot.set_visible(False)
        
        # Connect events
        self.canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.canvas.mpl_connect('axes_leave_event', self.on_leave_axes)

    def browse_file(self):
        """Modern file browser with better UX"""
        filetypes = [("MOL files", "*.mol"), ("All files", "*.*")]
        filename = filedialog.askopenfilename(
            title="Select MOL file",
            filetypes=filetypes,
            initialdir=os.getcwd()
        )
        if filename:
            self.file_path.set(filename)
            self.update_status("File loaded", "success")
            self.display_mol_content(filename)

    def display_mol_content(self, filename):
        """Display molecule with modern styling"""
        try:
            mol = Chem.MolFromMolFile(filename)
            if mol is None:
                self.show_error("Failed to read MOL file")
                return

            self.current_mol = mol
            
            # Generate 2D coordinates if needed
            if not mol.GetNumConformers():
                AllChem.Compute2DCoords(mol)

            # Clear and update molecule viewer
            self.mol_ax.clear()
            self.mol_ax.set_facecolor(self.colors['card_bg'])
            
            # Create molecule image with modern styling
            img = Draw.MolToImage(mol, size=(400, 400), 
                                 kekulize=True, 
                                 wedgeBonds=True,
                                 highlightAtoms=[],
                                 highlightBonds=[])
            
            self.mol_ax.imshow(img)
            self.mol_ax.axis('off')
            self.mol_canvas.draw()

            # Update molecule info with modern card design
            formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
            weight = Chem.rdMolDescriptors.CalcExactMolWt(mol)
            
            info_text = f"Formula: {formula}\nMolecular Weight: {weight:.2f} g/mol\nAtoms: {mol.GetNumAtoms()}\nBonds: {mol.GetNumBonds()}"
            self.mol_info.config(text=info_text, fg=self.colors['text_primary'])
            
            # Identify and display functional groups
            self.identify_and_display_groups(mol)
            
        except Exception as e:
            self.show_error(f"Failed to read MOL file: {str(e)}")

    def identify_and_display_groups(self, mol):
        """Identify and display functional groups with modern UI"""
        # Clear previous groups
        for widget in self.groups_frame.winfo_children():
            widget.destroy()
        
        # Identify functional groups
        groups = self.predictor.identify_functional_groups(mol)
        
        # Display groups as modern cards
        for group_name, present in groups.items():
            if present:
                group_card = tk.Frame(self.groups_frame,
                                    bg=self.colors['success'],
                                    relief='flat',
                                    bd=1)
                group_card.pack(fill='x', pady=2)
                
                group_label = tk.Label(group_card,
                                     text=f"✓ {group_name.replace('_', ' ').title()}",
                                     font=self.fonts['small'],
                                     fg='white',
                                     bg=self.colors['success'])
                group_label.pack(padx=10, pady=5)

    def predict_spectrum(self):
        """Predict spectrum with modern loading animation"""
        if not self.file_path.get():
            self.show_error("Please select a MOL file")
            return
        
        # Show loading state
        self.show_loading(True)
        self.update_status("Predicting spectrum...", "warning")
        
        # Run prediction in thread to avoid UI freezing
        threading.Thread(target=self._predict_spectrum_thread, daemon=True).start()

    def _predict_spectrum_thread(self):
        """Thread function for spectrum prediction"""
        try:
            mol_input = self.file_path.get()
            mol = self.predictor.read_molecule(mol_input, 'mol')
            
            if mol is None:
                self.root.after(0, lambda: self.show_error("Failed to read molecule"))
                return

            # Get functional groups
            self.current_molecule_groups = self.predictor.identify_functional_groups(mol)
            
            # Predict spectrum
            spectrum = self.predictor.predict_ftir(mol_input, input_format='mol')
            
            if spectrum is not None:
                # Update UI in main thread
                self.root.after(0, lambda: self._update_spectrum_display(spectrum, mol))
            else:
                self.root.after(0, lambda: self.show_error("Failed to predict spectrum"))
                
        except Exception as e:
            self.root.after(0, lambda: self.show_error(str(e)))
        finally:
            self.root.after(0, lambda: self.show_loading(False))

    def _update_spectrum_display(self, spectrum, mol):
        try:
            # Prepare peak details
            self.prepare_peak_details(mol)
            
            # Store spectrum data
            self.spectrum = spectrum  # Spectrum verisini sınıf değişkeni olarak sakla
            
            # Convert to transmittance
            transmittance = 100 * np.power(10, -np.array(spectrum['Absorbance']))
            wavenumbers = np.array(spectrum['Wavenumber'])
            
            # Update plot with modern styling
            self.ax.clear()
            self.ax.set_facecolor(self.colors['card_bg'])
            
            # Plot spectrum with gradient effect
            self.ax.plot(wavenumbers, transmittance,
                        color=self.colors['primary'], linewidth=2.5, alpha=0.9)
            
            # Fill under curve for modern look
            self.ax.fill_between(wavenumbers, transmittance, 0,
                               color=self.colors['primary'], alpha=0.1)
            
            # Modern plot styling
            self.ax.set_xlabel('Wavenumber (cm⁻¹)', fontsize=12, 
                             color=self.colors['text_primary'], fontweight='bold')
            self.ax.set_ylabel('Transmittance (%)', fontsize=12, 
                             color=self.colors['text_primary'], fontweight='bold')
            self.ax.grid(True, linestyle='--', alpha=0.3, color=self.colors['text_secondary'])
            self.ax.set_xlim(4000, 400)
            self.ax.set_ylim(0, 100)
            
            # Style spines and ticks
            for spine in self.ax.spines.values():
                spine.set_color(self.colors['text_secondary'])
                spine.set_linewidth(1)
            self.ax.spines['top'].set_visible(False)
            self.ax.spines['right'].set_visible(False)
            self.ax.tick_params(colors=self.colors['text_secondary'], labelsize=10)
            
            # Store peaks for interaction
            self.peaks = wavenumbers
            
            # Redraw
            self.canvas.draw()
            
            # Update status
            self.update_status("Spectrum predicted successfully", "success")
            
        except Exception as e:
            self.show_error(f"Error updating spectrum: {str(e)}")

    def show_loading(self, show):
        """Show/hide loading indicator"""
        if show:
            self.loading_frame.pack(fill='x', pady=5)
            self.animate_loading()
        else:
            self.loading_frame.pack_forget()
            
    def animate_loading(self):
        """Simple loading animation"""
        if self.loading_frame.winfo_viewable():
            current_text = self.loading_label.cget('text')
            dots = current_text.count('.')
            if dots >= 3:
                new_text = "⟳ Analyzing molecule"
            else:
                new_text = current_text + "."
            self.loading_label.config(text=new_text)
            self.root.after(500, self.animate_loading)

    def update_status(self, message, status_type="info"):
        """Update status indicator with color coding"""
        colors = {
            "success": self.colors['success'],
            "warning": self.colors['warning'],
            "error": self.colors['danger'],
            "info": self.colors['secondary']
        }
        
        symbols = {
            "success": "●",
            "warning": "⚠",
            "error": "✕",
            "info": "ℹ"
        }
        
        color = colors.get(status_type, self.colors['secondary'])
        symbol = symbols.get(status_type, "●")
        
        self.status_indicator.config(text=f"{symbol} {message}", fg=color)

    def show_error(self, message):
        """Show error with modern styling"""
        self.update_status("Error occurred", "error")
        messagebox.showerror("Error", message)

    def create_modern_button(self, parent, text, command, color=None):
        """Create a modern styled button"""
        if color is None:
            color = self.colors['primary']
        
        btn = tk.Button(parent,
                       text=text,
                       command=command,
                       font=self.fonts['button'],
                       bg=color,
                       fg='white',
                       relief='flat',
                       bd=0,
                       padx=20,
                       pady=10,
                       cursor='hand2',
                       activebackground=self.colors['primary_dark'],
                       activeforeground='white')
        
        # Add hover effects
        btn.bind('<Enter>', lambda e: btn.configure(bg=self.colors['primary_dark']))
        btn.bind('<Leave>', lambda e: btn.configure(bg=color))
        
        return btn

    def prepare_peak_details(self, mol):
        """Prepare detailed peak information"""
        self.peak_details = {}
        
        # Add general peak assignments
        self.peak_details.update(self.peak_assignments)
        
        # Add specific peak details based on functional groups
        if self.current_molecule_groups:
            if self.current_molecule_groups.get('ester'):
                self.peak_details.update({
                    (1750, 1730): {
                        'group': 'Ester C=O',
                        'description': 'Ester carbonyl stretch',
                        'intensity': 'Strong',
                        'structural_info': 'R-O-C(=O)-R'
                    },
                    (1300, 1000): {
                        'group': 'Ester C-O',
                        'description': 'Ester C-O stretch',
                        'intensity': 'Strong',
                        'structural_info': 'R-O-C(=O)-R'
                    }
                })
            
            if self.current_molecule_groups.get('alcohol'):
                self.peak_details.update({
                    (3200, 3600): {
                        'group': 'O-H stretch',
                        'description': 'Alcohol O-H stretch',
                        'intensity': 'Broad, strong',
                        'structural_info': 'R-OH'
                    },
                    (1050, 1200): {
                        'group': 'C-O stretch',
                        'description': 'Alcohol C-O stretch',
                        'intensity': 'Strong',
                        'structural_info': 'R-OH'
                    }
                })
            
            if self.current_molecule_groups.get('carbonyl'):
                self.peak_details.update({
                    (1670, 1780): {
                        'group': 'C=O stretch',
                        'description': 'Carbonyl stretch',
                        'intensity': 'Very strong',
                        'structural_info': 'Depends on carbonyl type'
                    }
                })
            
            if self.current_molecule_groups.get('aromatic'):
                self.peak_details.update({
                    (3080, 3030): {
                        'group': 'Aromatic C-H',
                        'description': 'Aromatic C-H stretch',
                        'intensity': 'Weak to medium',
                        'structural_info': 'Ar-H'
                    },
                    (1600, 1585): {
                        'group': 'Aromatic C=C',
                        'description': 'Aromatic ring stretch',
                        'intensity': 'Medium',
                        'structural_info': 'Conjugated ring'
                    },
                    (850, 700): {
                        'group': 'Aromatic C-H',
                        'description': 'Out-of-plane bending',
                        'intensity': 'Strong',
                        'structural_info': 'Adjacent H atoms pattern'
                    }
                })
            
            if self.current_molecule_groups.get('amine'):
                self.peak_details.update({
                    (3300, 3500): {
                        'group': 'N-H stretch',
                        'description': 'Amine N-H stretch',
                        'intensity': 'Medium',
                        'structural_info': 'Primary/secondary amine'
                    },
                    (1000, 1250): {
                        'group': 'C-N stretch',
                        'description': 'Amine C-N stretch',
                        'intensity': 'Medium',
                        'structural_info': 'R-NH2 or R2NH'
                    }
                })
                
    def on_mouse_move(self, event):
        """Handle mouse movement over spectrum"""
        if event.inaxes == self.ax and hasattr(self, 'peaks'):
            x = event.xdata
            y = event.ydata
            
            if x is None or y is None:
                return
            
            found_peak = False
            for (wave_min, wave_max), details in self.peak_details.items():
                if wave_min <= x <= wave_max:
                    description = (
                        f"🔬 {details['group']}\n"
                        f"📊 Range: {wave_min}-{wave_max} cm⁻¹\n"
                        f"⚡ Intensity: {details['intensity']}\n"
                        f"🧪 Structure: {details['structural_info']}"
                    )
                    
                    self.annot.xy = (x, y)
                    self.annot.set_text(description)
                    if not self.annot.get_visible():
                        self.annot.set_visible(True)
                        self.canvas.draw_idle()
                    found_peak = True
                    break
            
            if not found_peak and self.annot.get_visible():
                self.annot.set_visible(False)
                self.canvas.draw_idle()

    def export_spectrum(self):
        """Export spectrum data"""
        if not hasattr(self, 'peaks'):
            self.show_error("No spectrum to export")
            return
        
        filetypes = [
            ("PNG Image", "*.png"),
            ("PDF Document", "*.pdf"),
            ("SVG Vector", "*.svg"),
            ("CSV Data", "*.csv")
        ]
        
        filename = filedialog.asksaveasfilename(
            title="Export Spectrum",
            filetypes=filetypes,
            defaultextension=".png"
        )
        
        if filename:
            try:
                if filename.endswith('.csv'):
                    data = {
                        'Wavenumber': self.peaks,
                        'Transmittance': 100 * np.power(10, -np.array(self.spectrum['Absorbance']))
                    }
                    pd.DataFrame(data).to_csv(filename, index=False)
                else:
                    self.fig.savefig(filename, dpi=300, bbox_inches='tight',
                                   facecolor=self.colors['card_bg'])
                
                self.update_status(f"Exported to {os.path.basename(filename)}", "success")
            except Exception as e:
                self.show_error(f"Export failed: {str(e)}")

    def initialize_complete_gui(self):
        """Complete GUI initialization"""
        self.add_menu_bar()
        self.add_keyboard_shortcuts()
        self.create_context_menu()
        self.create_status_bar()
        self.apply_final_styling()

    def add_keyboard_shortcuts(self):
        """Add keyboard shortcuts"""
        self.root.bind('<Control-o>', lambda e: self.browse_file())
        self.root.bind('<Control-p>', lambda e: self.predict_spectrum())
        self.root.bind('<Control-s>', lambda e: self.export_spectrum())
        self.root.bind('<F5>', lambda e: self.refresh_display())
        self.root.bind('<Escape>', lambda e: self.root.quit())

    def create_context_menu(self):
        """Create right-click context menu"""
        context_menu = tk.Menu(self.root, tearoff=0)
        context_menu.add_command(label="Open File", command=self.browse_file)
        context_menu.add_command(label="Predict", command=self.predict_spectrum)
        context_menu.add_separator()
        context_menu.add_command(label="Export", command=self.export_spectrum)
        context_menu.add_command(label="Refresh", command=self.refresh_display)
        context_menu.add_separator()
        context_menu.add_command(label="Exit", command=self.root.quit)
        
        def show_context_menu(event):
            try:
                context_menu.tk_popup(event.x_root, event.y_root)
            finally:
                context_menu.grab_release()
        
        self.root.bind("<Button-3>", show_context_menu)

    def on_leave_axes(self, event):
        """Handle mouse leaving the plot area"""
        # Hide the annotation when mouse leaves the plot
        if hasattr(self, 'annot') and self.annot.get_visible():
            self.annot.set_visible(False)
            self.canvas.draw_idle()        

def main():
    try:
        root = tk.Tk()
        app = ModernFTIRGUI(root)
        
        # Set window properties
        try:
            root.state('zoomed')  # Windows
        except:
            root.attributes('-zoomed', True)  # Linux/Mac
        
        root.mainloop()
    except Exception as e:
        print(f"Error starting application: {e}")

if __name__ == "__main__":
    main()