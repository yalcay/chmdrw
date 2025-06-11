import os
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
from io import BytesIO
import xml.etree.ElementTree as ET

# RDKit uyarılarını gizle
RDLogger.DisableLog('rdApp.*')

# RDKit'te molekül oluşturup kaydetme (alternatif yöntem)
mol = Chem.MolFromSmiles('C=O')  # Örnek formaldehit
AllChem.Compute2DCoords(mol)
Chem.MolToMolFile(mol, 'output.mol', forceV3000=False)  # V2000 formatında kaydet

class FTIRPredictor:
    def __init__(self):
        self.wavenumbers = np.linspace(400, 4000, 1000)
        self.ir_data = self.load_complete_ir_data()
        
    def load_complete_ir_data(self):
        """Excel'deki tüm FTIR verilerini yükle"""
        ir_data = {
            'alcohol': [
                {'range': (3700, 3584), 'center': 3642, 'width': 58, 'height': self.determine_height('medium, sharp'), 
                 'appearance': 'medium, sharp', 'group': 'O-H stretching', 'compound_class': 'alcohol', 'comments': 'free'},
                {'range': (3550, 3200), 'center': 3375, 'width': 175, 'height': self.determine_height('strong, broad'), 
                 'appearance': 'strong, broad', 'group': 'O-H stretching', 'compound_class': 'alcohol', 'comments': 'intermolecular bonded'},
                {'range': (3200, 2700), 'center': 2950, 'width': 250, 'height': self.determine_height('weak, broad'), 
                 'appearance': 'weak, broad', 'group': 'O-H stretching', 'compound_class': 'alcohol', 'comments': 'intramolecular bonded'}
            ],
            'amine': [
                {'range': (3500, 3500), 'center': 3500, 'width': 20, 'height': self.determine_height('medium'), 
                 'appearance': 'medium', 'group': 'N-H stretching', 'compound_class': 'primary amine', 'comments': ''},
                {'range': (3400, 3300), 'center': 3350, 'width': 50, 'height': self.determine_height('medium'), 
                 'appearance': 'medium', 'group': 'N-H stretching', 'compound_class': 'aliphatic primary amine', 'comments': ''},
                {'range': (3350, 3310), 'center': 3330, 'width': 20, 'height': self.determine_height('medium'), 
                 'appearance': 'medium', 'group': 'N-H stretching', 'compound_class': 'secondary amine', 'comments': ''},
                {'range': (3000, 2800), 'center': 2900, 'width': 100, 'height': self.determine_height('strong, broad'), 
                 'appearance': 'strong, broad', 'group': 'N-H stretching', 'compound_class': 'amine salt', 'comments': ''}
            ],
            'carboxylic_acid': [
                {'range': (3300, 2500), 'center': 3000, 'width': 400, 'height': self.determine_height('strong, broad'), 
                 'appearance': 'strong, broad', 'group': 'O-H stretching', 'compound_class': 'carboxylic acid', 'comments': 'usually centered on 3000 cm-1'},
                {'range': (1720, 1706), 'center': 1713, 'width': 7, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=O stretching', 'compound_class': 'carboxylic acid', 'comments': 'dimer'},
                {'range': (1760, 1760), 'center': 1760, 'width': 10, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=O stretching', 'compound_class': 'carboxylic acid', 'comments': 'monomer'},
                {'range': (1440, 1395), 'center': 1417.5, 'width': 22.5, 'height': self.determine_height('medium'), 
                 'appearance': 'medium', 'group': 'O-H bending', 'compound_class': 'carboxylic acid', 'comments': ''}
            ],
            'alkyne': [
                {'range': (3333, 3267), 'center': 3300, 'width': 33, 'height': self.determine_height('strong, sharp'), 
                 'appearance': 'strong, sharp', 'group': 'C-H stretching', 'compound_class': 'alkyne', 'comments': ''},
                {'range': (2260, 2190), 'center': 2225, 'width': 35, 'height': self.determine_height('weak'), 
                 'appearance': 'weak', 'group': 'C≡C stretching', 'compound_class': 'alkyne', 'comments': 'disubstituted'},
                {'range': (2140, 2100), 'center': 2120, 'width': 20, 'height': self.determine_height('weak'), 
                 'appearance': 'weak', 'group': 'C≡C stretching', 'compound_class': 'alkyne', 'comments': 'monosubstituted'}
            ],
            'alkene': [
                {'range': (3100, 3000), 'center': 3050, 'width': 50, 'height': self.determine_height('medium'), 
                 'appearance': 'medium', 'group': 'C-H stretching', 'compound_class': 'alkene', 'comments': ''},
                {'range': (1678, 1668), 'center': 1673, 'width': 5, 'height': self.determine_height('weak'), 
                 'appearance': 'weak', 'group': 'C=C stretching', 'compound_class': 'alkene', 'comments': 'disubstituted (trans)'},
                {'range': (1675, 1665), 'center': 1670, 'width': 5, 'height': self.determine_height('weak'), 
                 'appearance': 'weak', 'group': 'C=C stretching', 'compound_class': 'alkene', 'comments': 'trisubstituted'},
                {'range': (1675, 1665), 'center': 1670, 'width': 5, 'height': self.determine_height('weak'), 
                 'appearance': 'weak', 'group': 'C=C stretching', 'compound_class': 'alkene', 'comments': 'tetrasubstituted'},
                {'range': (1662, 1626), 'center': 1644, 'width': 18, 'height': self.determine_height('medium'), 
                 'appearance': 'medium', 'group': 'C=C stretching', 'compound_class': 'alkene', 'comments': 'disubstituted (cis)'},
                {'range': (1658, 1648), 'center': 1653, 'width': 5, 'height': self.determine_height('medium'), 
                 'appearance': 'medium', 'group': 'C=C stretching', 'compound_class': 'alkene', 'comments': 'vinylidene'},
                {'range': (1650, 1600), 'center': 1625, 'width': 25, 'height': self.determine_height('medium'), 
                 'appearance': 'medium', 'group': 'C=C stretching', 'compound_class': 'conjugated alkene', 'comments': ''},
                {'range': (1648, 1638), 'center': 1643, 'width': 5, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=C stretching', 'compound_class': 'alkene', 'comments': 'monosubstituted'},
                {'range': (995, 985), 'center': 990, 'width': 5, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=C bending', 'compound_class': 'alkene', 'comments': 'monosubstituted'},
                {'range': (980, 960), 'center': 970, 'width': 10, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=C bending', 'compound_class': 'alkene', 'comments': 'disubstituted (trans)'},
                {'range': (895, 885), 'center': 890, 'width': 5, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=C bending', 'compound_class': 'alkene', 'comments': 'vinylidene'},
                {'range': (840, 790), 'center': 815, 'width': 25, 'height': self.determine_height('medium'), 
                 'appearance': 'medium', 'group': 'C=C bending', 'compound_class': 'alkene', 'comments': 'trisubstituted'},
                {'range': (730, 665), 'center': 697.5, 'width': 32.5, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=C bending', 'compound_class': 'alkene', 'comments': 'disubstituted (cis)'}
            ],
            'alkane': [
                {'range': (3000, 2840), 'center': 2920, 'width': 80, 'height': self.determine_height('medium'), 
                 'appearance': 'medium', 'group': 'C-H stretching', 'compound_class': 'alkane', 'comments': ''},
                {'range': (1465, 1465), 'center': 1465, 'width': 10, 'height': self.determine_height('medium'), 
                 'appearance': 'medium', 'group': 'C-H bending', 'compound_class': 'alkane', 'comments': 'methylene group'},
                {'range': (1450, 1450), 'center': 1450, 'width': 10, 'height': self.determine_height('medium'), 
                 'appearance': 'medium', 'group': 'C-H bending', 'compound_class': 'alkane', 'comments': 'methyl group'},
                {'range': (1385, 1380), 'center': 1382.5, 'width': 2.5, 'height': self.determine_height('medium'), 
                 'appearance': 'medium', 'group': 'C-H bending', 'compound_class': 'alkane', 'comments': 'gem dimethyl'}
            ],
            'aldehyde': [
                {'range': (2830, 2695), 'center': 2762.5, 'width': 67.5, 'height': self.determine_height('medium'), 
                 'appearance': 'medium', 'group': 'C-H stretching', 'compound_class': 'aldehyde', 'comments': 'doublet'},
                {'range': (1740, 1720), 'center': 1730, 'width': 10, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=O stretching', 'compound_class': 'aldehyde', 'comments': 'normal aldehyde'},
                {'range': (1710, 1685), 'center': 1697.5, 'width': 12.5, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=O stretching', 'compound_class': 'aldehyde', 'comments': 'conjugated aldehyde'},
                {'range': (1390, 1380), 'center': 1385, 'width': 5, 'height': self.determine_height('medium'), 
                 'appearance': 'medium', 'group': 'C-H bending', 'compound_class': 'aldehyde', 'comments': ''}
            ],
            'ketone': [
                {'range': (1725, 1705), 'center': 1715, 'width': 10, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=O stretching', 'compound_class': 'aliphatic ketone', 'comments': 'or cyclohexanone or cyclopentenone'},
                {'range': (1685, 1666), 'center': 1675.5, 'width': 9.5, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=O stretching', 'compound_class': 'conjugated ketone', 'comments': ''},
                {'range': (1620, 1610), 'center': 1615, 'width': 5, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=C stretching', 'compound_class': 'α,β-unsaturated ketone', 'comments': ''}
            ],
            'ester': [
                {'range': (1750, 1735), 'center': 1742.5, 'width': 7.5, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=O stretching', 'compound_class': 'esters', 'comments': '6-membered lactone'},
                {'range': (1750, 1735), 'center': 1742.5, 'width': 7.5, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=O stretching', 'compound_class': 'δ-lactone', 'comments': 'γ: 1770'},
                {'range': (1770, 1780), 'center': 1775, 'width': 5, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=O stretching', 'compound_class': 'vinyl / phenyl ester', 'comments': ''},
                {'range': (1730, 1715), 'center': 1722.5, 'width': 7.5, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=O stretching', 'compound_class': 'α,β-unsaturated ester', 'comments': 'or formates'},
                {'range': (1210, 1163), 'center': 1186.5, 'width': 23.5, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C-O stretching', 'compound_class': 'ester', 'comments': ''}
            ],
            'amide': [
                {'range': (1690, 1690), 'center': 1690, 'width': 20, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=O stretching', 'compound_class': 'primary amide', 'comments': 'free (associated: 1650)'},
                {'range': (1680, 1680), 'center': 1680, 'width': 20, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=O stretching', 'compound_class': 'secondary amide', 'comments': 'free (associated: 1640)'},
                {'range': (1680, 1680), 'center': 1680, 'width': 20, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=O stretching', 'compound_class': 'tertiary amide', 'comments': 'free (associated: 1630)'},
                {'range': (1650, 1650), 'center': 1650, 'width': 20, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=O stretching', 'compound_class': 'δ-lactam', 'comments': ''},
                {'range': (1750, 1700), 'center': 1725, 'width': 25, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=O stretching', 'compound_class': 'γ-lactam', 'comments': ''},
                {'range': (1760, 1730), 'center': 1745, 'width': 15, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=O stretching', 'compound_class': 'β-lactam', 'comments': ''},
                {'range': (1650, 1580), 'center': 1615, 'width': 35, 'height': self.determine_height('medium'), 
                 'appearance': 'medium', 'group': 'N-H bending', 'compound_class': 'amine', 'comments': ''}
            ],
            'nitrile': [
                {'range': (2260, 2222), 'center': 2241, 'width': 19, 'height': self.determine_height('weak'), 
                 'appearance': 'weak', 'group': 'C≡N stretching', 'compound_class': 'nitrile', 'comments': ''}
            ],
            'thiol': [
                {'range': (2600, 2550), 'center': 2575, 'width': 25, 'height': self.determine_height('weak'), 
                 'appearance': 'weak', 'group': 'S-H stretching', 'compound_class': 'thiol', 'comments': ''}
            ],
            'ether': [
                {'range': (1275, 1200), 'center': 1237.5, 'width': 37.5, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C-O stretching', 'compound_class': 'alkyl aryl ether', 'comments': ''},
                {'range': (1225, 1200), 'center': 1212.5, 'width': 12.5, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C-O stretching', 'compound_class': 'vinyl ether', 'comments': ''},
                {'range': (1150, 1085), 'center': 1117.5, 'width': 32.5, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C-O stretching', 'compound_class': 'aliphatic ether', 'comments': ''}
            ],
            'alcohol_oh_bending': [
                {'range': (1420, 1330), 'center': 1375, 'width': 45, 'height': self.determine_height('medium'), 
                 'appearance': 'medium', 'group': 'O-H bending', 'compound_class': 'alcohol', 'comments': ''},
                {'range': (1205, 1124), 'center': 1164.5, 'width': 40.5, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C-O stretching', 'compound_class': 'tertiary alcohol', 'comments': ''},
                {'range': (1124, 1087), 'center': 1105.5, 'width': 18.5, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C-O stretching', 'compound_class': 'secondary alcohol', 'comments': ''},
                {'range': (1085, 1050), 'center': 1067.5, 'width': 17.5, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C-O stretching', 'compound_class': 'primary alcohol', 'comments': ''}
            ],
            'anhydride': [
                {'range': (1818, 1818), 'center': 1818, 'width': 20, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=O stretching', 'compound_class': 'anhydride', 'comments': ''},
                {'range': (1775, 1775), 'center': 1775, 'width': 20, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=O stretching', 'compound_class': 'conjugated anhydride', 'comments': ''},
                {'range': (1050, 1040), 'center': 1045, 'width': 5, 'height': self.determine_height('strong, broad'), 
                 'appearance': 'strong, broad', 'group': 'CO-O-CO stretching', 'compound_class': 'anhydride', 'comments': ''}
            ],
            'acid_halide': [
                {'range': (1815, 1785), 'center': 1800, 'width': 15, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=O stretching', 'compound_class': 'acid halide', 'comments': ''},
                {'range': (1800, 1770), 'center': 1785, 'width': 15, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=O stretching', 'compound_class': 'conjugated acid halide', 'comments': ''}
            ],
            'carbon_dioxide': [
                {'range': (2349, 2349), 'center': 2349, 'width': 20, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'O=C=O stretching', 'compound_class': 'carbon dioxide', 'comments': ''}
            ],
            'isocyanate': [
                {'range': (2275, 2250), 'center': 2262.5, 'width': 12.5, 'height': self.determine_height('strong, broad'), 
                 'appearance': 'strong, broad', 'group': 'N=C=O stretching', 'compound_class': 'isocyanate', 'comments': ''}
            ],
            'thiocyanate': [
                {'range': (2175, 2140), 'center': 2157.5, 'width': 17.5, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'S-C≡N stretching', 'compound_class': 'thiocyanate', 'comments': ''}
            ],
            'azide': [
                {'range': (2160, 2120), 'center': 2140, 'width': 20, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'N=N=N stretching', 'compound_class': 'azide', 'comments': ''}
            ],
            'ketene': [
                {'range': (2150, 2150), 'center': 2150, 'width': 20, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=C=O stretching', 'compound_class': 'ketene', 'comments': ''}
            ],
            'carbodiimide': [
                {'range': (2145, 2120), 'center': 2132.5, 'width': 12.5, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'N=C=N stretching', 'compound_class': 'carbodiimide', 'comments': ''}
            ],
            'isothiocyanate': [
                {'range': (2140, 1990), 'center': 2065, 'width': 75, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'N=C=S stretching', 'compound_class': 'isothiocyanate', 'comments': ''}
            ],
            'allene': [
                {'range': (2000, 1900), 'center': 1950, 'width': 50, 'height': self.determine_height('medium'), 
                 'appearance': 'medium', 'group': 'C=C=C stretching', 'compound_class': 'allene', 'comments': ''}
            ],
            'ketenimine': [
                {'range': (2000, 2000), 'center': 2000, 'width': 20, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C=C=N stretching', 'compound_class': 'ketenimine', 'comments': ''}
            ],
            'aromatic': [
                {'range': (2000, 1650), 'center': 1825, 'width': 175, 'height': self.determine_height('weak'), 
                 'appearance': 'weak', 'group': 'C-H bending', 'compound_class': 'aromatic compound', 'comments': 'overtone'},
                {'range': (1600, 1585), 'center': 1592.5, 'width': 7.5, 'height': self.determine_height('medium'), 
                 'appearance': 'medium', 'group': 'C=C stretching', 'compound_class': 'aromatic', 'comments': ''},
                {'range': (1500, 1500), 'center': 1500, 'width': 20, 'height': self.determine_height('medium'), 
                 'appearance': 'medium', 'group': 'C=C stretching', 'compound_class': 'aromatic', 'comments': ''},
                {'range': (880, 860), 'center': 870, 'width': 10, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C-H bending', 'compound_class': '1,2,4-trisubstituted', 'comments': ''},
                {'range': (880, 860), 'center': 870, 'width': 10, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C-H bending', 'compound_class': '1,3-disubstituted', 'comments': ''},
                {'range': (810, 790), 'center': 800, 'width': 10, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C-H bending', 'compound_class': '1,4-disubstituted', 'comments': ''},
                {'range': (780, 760), 'center': 770, 'width': 10, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C-H bending', 'compound_class': '1,2,3-trisubstituted', 'comments': ''},
                {'range': (755, 755), 'center': 755, 'width': 10, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C-H bending', 'compound_class': '1,2-disubstituted', 'comments': ''},
                {'range': (750, 750), 'center': 750, 'width': 10, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C-H bending', 'compound_class': 'monosubstituted', 'comments': ''}
            ],
            'nitro': [
                {'range': (1550, 1500), 'center': 1525, 'width': 25, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'N-O stretching', 'compound_class': 'nitro compound', 'comments': ''}
            ],
            'sulfur_compounds': [
                {'range': (1415, 1380), 'center': 1397.5, 'width': 17.5, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'S=O stretching', 'compound_class': 'sulfate', 'comments': ''},
                {'range': (1410, 1380), 'center': 1395, 'width': 15, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'S=O stretching', 'compound_class': 'sulfonyl chloride', 'comments': ''},
                {'range': (1372, 1335), 'center': 1353.5, 'width': 18.5, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'S=O stretching', 'compound_class': 'sulfonate', 'comments': ''},
                {'range': (1370, 1335), 'center': 1352.5, 'width': 17.5, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'S=O stretching', 'compound_class': 'sulfonamide', 'comments': ''},
                {'range': (1350, 1342), 'center': 1346, 'width': 4, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'S=O stretching', 'compound_class': 'sulfonic acid', 'comments': 'anhydrous'},
                {'range': (1350, 1300), 'center': 1325, 'width': 25, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'S=O stretching', 'compound_class': 'sulfone', 'comments': ''},
                {'range': (1070, 1030), 'center': 1050, 'width': 20, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'S=O stretching', 'compound_class': 'sulfoxide', 'comments': ''}
            ],
            'halo_compounds': [
                {'range': (1400, 1000), 'center': 1200, 'width': 200, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C-F stretching', 'compound_class': 'fluoro compound', 'comments': ''},
                {'range': (850, 550), 'center': 700, 'width': 150, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C-Cl stretching', 'compound_class': 'halo compound', 'comments': ''},
                {'range': (690, 515), 'center': 602.5, 'width': 87.5, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C-Br stretching', 'compound_class': 'halo compound', 'comments': ''},
                {'range': (600, 500), 'center': 550, 'width': 50, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C-I stretching', 'compound_class': 'halo compound', 'comments': ''}
            ],
            'amine_cn_stretching': [
                {'range': (1342, 1266), 'center': 1304, 'width': 38, 'height': self.determine_height('strong'), 
                 'appearance': 'strong', 'group': 'C-N stretching', 'compound_class': 'aromatic amine', 'comments': ''},
                {'range': (1250, 1020), 'center': 1135, 'width': 115, 'height': self.determine_height('medium'), 
                 'appearance': 'medium', 'group': 'C-N stretching', 'compound_class': 'amine', 'comments': ''}
            ],
            'phenol': [
                {'range': (3650, 3600), 'center': 3625, 'width': 25, 'height': self.determine_height('strong, sharp'),
                 'appearance': 'strong, sharp', 'group': 'O-H stretching', 'compound_class': 'phenol', 'comments': 'free'},
                {'range': (3400, 3200), 'center': 3300, 'width': 100, 'height': self.determine_height('strong, broad'),
                 'appearance': 'strong, broad', 'group': 'O-H stretching', 'compound_class': 'phenol', 'comments': 'H-bonded'},
                {'range': (1260, 1180), 'center': 1220, 'width': 40, 'height': self.determine_height('strong'),
                 'appearance': 'strong', 'group': 'C-O stretching', 'compound_class': 'phenol', 'comments': ''}
            ],
            'imine': [
                {'range': (1690, 1640), 'center': 1665, 'width': 25, 'height': self.determine_height('medium'), 
                 'appearance': 'medium', 'group': 'C=N stretching', 'compound_class': 'imine / oxime', 'comments': ''}
            ]
        }
        
        return ir_data
    
    def determine_height(self, appearance):
        """Görünüm bilgisine göre pik yüksekliğini belirle"""
        appearance = str(appearance).lower()
        
        if 'strong' in appearance and 'broad' in appearance:
            return 0.9
        elif 'strong' in appearance and 'sharp' in appearance:
            return 0.85
        elif 'strong' in appearance:
            return 0.8
        elif 'medium' in appearance and 'broad' in appearance:
            return 0.6
        elif 'medium' in appearance and 'sharp' in appearance:
            return 0.55
        elif 'medium' in appearance:
            return 0.5
        elif 'weak' in appearance and 'broad' in appearance:
            return 0.3
        elif 'weak' in appearance and 'sharp' in appearance:
            return 0.25
        elif 'weak' in appearance:
            return 0.2
        else:
            return 0.5  # varsayılan
    
    def identify_functional_groups(self, mol):
        """Moleküldeki tüm fonksiyonel grupları tespit et"""
        functional_groups = {
            'alcohol': False,
            'amine': False,
            'carboxylic_acid': False,
            'ester': False,
            'aldehyde': False,
            'ketone': False,
            'ether': False,
            'alkene': False,
            'alkyne': False,
            'aromatic': False,
            'amide': False,
            'nitrile': False,
            'thiol': False,
            'halide': False,
            'nitro': False,
            'sulfoxide': False,
            'sulfone': False,
            'anhydride': False,
            'imide': False,
            'azide': False,
            'isocyanate': False,
            'thiocyanate': False,
            'ketene': False,
            'carbodiimide': False,
            'isothiocyanate': False,
            'phenol': False,
            'allene': False,
            'ketenimine': False
        }
        
        # SMARTS desenleri
        patterns = {
            'alcohol': '[OH1][CX4]',
            'amine': '[NX3;H2,H1;!$(NC=O)]',
            'carboxylic_acid': '[CX3](=O)[OH]',
            'ester': '[CX3](=O)[OX2H0][#6]',
            'aldehyde': '[CX3H1](=O)[#6]',
            'ketone': '[CX3](=O)[#6]',
            'ether': '[OD2]([#6])[#6]',
            'alkene': '[CX3]=[CX3]',
            'alkyne': '[CX2]#[CX2]',
            'aromatic': 'c1ccccc1',
            'amide': '[NX3][CX3](=[OX1])[#6]',
            'nitrile': '[NX1]#[CX2]',
            'thiol': '[SH]',
            'halide': '[F,Cl,Br,I]',
            'nitro': '[NX3](=O)=O',
            'sulfoxide': '[SX3](=O)([#6])[#6]',
            'sulfone': '[SX4](=O)(=O)([#6])[#6]',
            'anhydride': '[CX3](=O)[OX2][CX3](=O)',
            'imide': '[CX3](=O)[NX3][CX3](=O)',
            'azide': '[NX3]=[NX2+]=[NX3-]',
            'isocyanate': '[NX2]=[CX2]=[OX1]',
            'thiocyanate': '[SX2][CX2]#[NX1]',
            'ketene': '[CX2]=[CX2]=[OX1]',
            'carbodiimide': '[NX2]=[CX2]=[NX2]',
            'isothiocyanate': '[NX2]=[CX2]=[SX1]',
            'allene': '[CX2]=[CX2]=[CX2]',
            'phenol': '[OH1][c]',  
            'ketenimine': '[CX2]=[CX2]=[NX2]'
        }
        
        # Her deseni kontrol et
        for group, pattern in patterns.items():
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                functional_groups[group] = True
                
        return functional_groups
    
    def generate_spectrum(self, mol):
        """Transmittance (%) tabanlı FTIR spektrumu oluşturur"""
        try:
            # 1. Temel spektrum ve atamaları hazırla
            absorbance = np.zeros_like(self.wavenumbers)
            peak_assignments = []
            
            # 2. Molekül geçerlilik kontrolü
            if mol is None or mol.GetNumAtoms() == 0:
                raise ValueError("Invalid molecule - no atoms present")
            
            # 3. Fonksiyonel grupları tespit et
            groups = self.identify_functional_groups(mol)
            if not groups:
                raise ValueError("No functional groups identified")
            
            # 4. Her grup için ilgili pikleri ekle
            for group, present in groups.items():
                if present and group in self.ir_data:
                    for peak_info in self.ir_data[group]:
                        try:
                            center = float(peak_info['center'])
                            width = float(peak_info['width'])
                            height = float(peak_info['height'])
                            
                            if not (400 <= center <= 4000):
                                continue  # Geçersiz dalga sayısı aralığı
                            
                            # Gauss pik ekle
                            peak = self.add_peak(
                                self.wavenumbers,
                                center,
                                width,
                                height
                            )
                            absorbance += peak
                            
                            # Pik atamalarını kaydet
                            peak_assignments.append({
                                'group': group,
                                'center': center,
                                'range': peak_info['range'],
                                'appearance': peak_info['appearance'],
                                'type': peak_info['compound_class'],
                                'height': height,
                                'comments': peak_info['comments']
                            })
                        except (ValueError, KeyError) as e:
                            print(f"Skipping invalid peak for {group}: {str(e)}")
                            continue
            
            # 5. Alifatik C-H gerilmeleri (her karbonlu molekül için)
            if any(atom.GetSymbol() == 'C' for atom in mol.GetAtoms()):
                # CH3 asimetrik
                absorbance += self.add_peak(self.wavenumbers, 2960, 30, 0.3)
                peak_assignments.append({
                    'group': 'alkane',
                    'center': 2960,
                    'range': (2940, 2980),
                    'appearance': 'medium',
                    'type': 'CH3 asymmetric stretch',
                    'height': 0.3,
                    'comments': ''
                })
                
                # CH2 asimetrik
                absorbance += self.add_peak(self.wavenumbers, 2930, 30, 0.3)
                peak_assignments.append({
                    'group': 'alkane',
                    'center': 2930,
                    'range': (2910, 2950),
                    'appearance': 'medium',
                    'type': 'CH2 asymmetric stretch',
                    'height': 0.3,
                    'comments': ''
                })
                
                # CH3 simetrik
                absorbance += self.add_peak(self.wavenumbers, 2870, 30, 0.2)
                peak_assignments.append({
                    'group': 'alkane',
                    'center': 2870,
                    'range': (2850, 2890),
                    'appearance': 'medium',
                    'type': 'CH3 symmetric stretch',
                    'height': 0.2,
                    'comments': ''
                })
            
            # 6. Normalizasyon ve gürültü ekleme
            if np.max(absorbance) > 0:
                absorbance = absorbance / np.max(absorbance)
            
            absorbance += np.random.normal(0, 0.01, len(self.wavenumbers))
            
            # 7. Absorbans'tan Transmittance'a dönüşüm
            transmittance = 100 * np.power(10, -np.array(absorbance))
            
            # 8. Son kontrol
            if np.all(absorbance == 0):
                raise ValueError("Generated spectrum is empty - no peaks added")
            
            return {
                'wavenumbers': self.wavenumbers.tolist(),
                'transmittance': transmittance.tolist(),
                'absorbance': absorbance.tolist(),
                'functional_groups': groups,
                'peak_assignments': peak_assignments
            }
            
        except Exception as e:
            print(f"Spectrum generation failed: {str(e)}")
            # Fallback mekanizması (temel C-H pikleri)
            basic_absorbance = np.zeros_like(self.wavenumbers)
            basic_absorbance += self.add_peak(self.wavenumbers, 2960, 30, 0.3)
            basic_absorbance += self.add_peak(self.wavenumbers, 2930, 30, 0.3)
            basic_absorbance += self.add_peak(self.wavenumbers, 2870, 30, 0.2)
            basic_transmittance = 100 * np.power(10, -basic_absorbance)
            
            return {
                'wavenumbers': self.wavenumbers.tolist(),
                'transmittance': basic_transmittance.tolist(),
                'absorbance': basic_absorbance.tolist(),
                'functional_groups': groups if 'groups' in locals() else {},
                'peak_assignments': [{
                    'group': 'alkane',
                    'center': 2960,
                    'range': (2940, 2980),
                    'appearance': 'medium',
                    'type': 'Basic C-H peaks',
                    'height': 0.3,
                    'comments': 'Fallback spectrum'
                }]
            }

    def predict_ftir(self, mol_input, input_format='mol'):
        """Geliştirilmiş FTIR tahmin fonksiyonu"""
        try:
            # Molekülü oku
            mol = self.read_molecule(mol_input, input_format)
            if mol is None:
                raise ValueError("Molecule could not be parsed from the input file")
            
            # Molekülün geçerli olup olmadığını kontrol et
            if mol.GetNumAtoms() == 0:
                raise ValueError("Molecule contains no atoms")
            
            # 2D koordinatları hesapla
            AllChem.Compute2DCoords(mol)
            
            # Spektrum oluştur
            spectrum_data = self.generate_spectrum(mol)
            
            if spectrum_data is None:
                raise ValueError("Spectrum generation failed")
            
            # DataFrame oluştur
            df = pd.DataFrame({
                'Wavenumber': spectrum_data['wavenumbers'],
                'Absorbance': spectrum_data['absorbance']
            })
            
            # Meta verileri ekle
            df.attrs['functional_groups'] = spectrum_data['functional_groups']
            df.attrs['peak_assignments'] = spectrum_data['peak_assignments']
            df.attrs['molecule'] = mol
            
            return df
            
        except Exception as e:
            print(f"Prediction error: {str(e)}")
            
            # Hata durumunda fallback mekanizması
            try:
                if mol is not None:
                    # Temel bir spektrum oluştur
                    basic_spectrum = np.zeros_like(self.wavenumbers)
                    if any(atom.GetSymbol() == 'C' for atom in mol.GetAtoms()):
                        basic_spectrum += self.add_peak(self.wavenumbers, 2960, 30, 0.3)
                        basic_spectrum += self.add_peak(self.wavenumbers, 2930, 30, 0.3)
                        basic_spectrum += self.add_peak(self.wavenumbers, 2870, 30, 0.2)
                    
                    df = pd.DataFrame({
                        'Wavenumber': self.wavenumbers.tolist(),
                        'Absorbance': basic_spectrum.tolist()
                    })
                    
                    groups = self.identify_functional_groups(mol) if mol else {}
                    df.attrs['functional_groups'] = groups
                    df.attrs['peak_assignments'] = [{
                        'group': 'alkane',
                        'center': 2960,
                        'range': (2940, 2980),
                        'appearance': 'medium',
                        'type': 'Basic C-H peaks',
                        'height': 0.3,
                        'comments': 'Fallback spectrum'
                    }]
                    df.attrs['molecule'] = mol
                    
                    return df
            except:
                pass
            
            return None
            
    def add_peak(self, x, center, width, height):
        """
        Add a Gaussian peak to the spectrum
        x: wavenumber array
        center: peak center
        width: peak width (sigma)
        height: peak height
        """
        # Convert width to sigma (standard deviation)
        sigma = width / 2.355  # FWHM to sigma conversion
        
        # Calculate Gaussian peak
        peak = height * np.exp(-((x - center) ** 2) / (2 * sigma ** 2))
        
        return peak
    
    def read_molecule(self, filepath):
        """Geliştirilmiş MOL dosyası okuma fonksiyonu"""
        try:
            # 1. Dosya formatını kontrol et
            if not os.path.exists(filepath):
                raise FileNotFoundError(f"File not found: {filepath}")
                
            file_ext = os.path.splitext(filepath)[1].lower()
            
            # 2. Dosya içeriğini oku
            with open(filepath, 'r', encoding='utf-8') as f:
                content = f.read()
                
            # 3. MOL formatı versiyonunu belirle
            mol = None
            if 'V2000' in content:
                # V2000 formatı
                mol = Chem.MolFromMolBlock(content, sanitize=False, removeHs=False)
            elif 'V3000' in content:
                # V3000 formatı
                mol = Chem.MolFromMolBlock(content, sanitize=False, removeHs=False)
            else:
                # Versiyon belirtilmemiş, V2000 varsayalım
                mol = Chem.MolFromMolFile(filepath, sanitize=False, removeHs=False)
                
            if mol is None:
                # SDF denemesi
                suppl = Chem.SDMolSupplier(filepath, sanitize=False, removeHs=False)
                mol = next(suppl, None)
                
            # 4. Molekül doğrulama ve temizleme
            if mol is not None:
                try:
                    # Hidrojenler korunarak sanitize
                    Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL^Chem.SANITIZE_ADJUSTHS)
                    
                    # Temel özellik kontrolü
                    if mol.GetNumAtoms() == 0:
                        raise ValueError("Molecule contains no atoms")
                        
                    # 2D koordinatları hesapla
                    AllChem.Compute2DCoords(mol)
                    
                    # Kerato kontrolü (özel durumlar için)
                    mol.UpdatePropertyCache(strict=False)
                    
                    return mol
                    
                except Exception as e:
                    print(f"Sanitization warning: {str(e)}")
                    # Sanitize başarısız olursa, minimum temizleme ile devam et
                    mol.UpdatePropertyCache(strict=False)
                    return mol
            
            raise ValueError("Could not parse molecule from file")
            
        except Exception as e:
            print(f"Error reading molecule: {str(e)}")
            return None

    def validate_mol_file(self, filepath):
        """MOL dosyası validasyonu"""
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                lines = f.readlines()
                
            # 1. Minimum satır sayısı kontrolü
            if len(lines) < 4:
                return False, "Invalid MOL file: Too few lines"
                
            # 2. Başlık satırları kontrolü
            header = lines[0].strip()
            if not header:
                return False, "Invalid MOL file: Missing header"
                
            # 3. Counts Line kontrolü (4. satır)
            counts_line = lines[3].strip()
            if len(counts_line) < 6:
                return False, "Invalid MOL file: Invalid counts line"
                
            try:
                # Atom ve bağ sayılarını çıkar
                atom_count = int(counts_line[0:3])
                bond_count = int(counts_line[3:6])
                
                # 4. Atom ve bağ sayıları kontrolü
                if atom_count <= 0:
                    return False, "Invalid MOL file: No atoms defined"
                    
                # 5. Dosya uzunluğu kontrolü
                expected_length = 4 + atom_count + bond_count
                if len(lines) < expected_length:
                    return False, "Invalid MOL file: File too short for atom/bond count"
                    
                # 6. Atom bloku formatı kontrolü
                for i in range(4, 4 + atom_count):
                    if len(lines[i]) < 31:  # Minimum atom satırı uzunluğu
                        return False, f"Invalid MOL file: Invalid atom line {i-3}"
                        
                # 7. Bağ bloku formatı kontrolü
                for i in range(4 + atom_count, 4 + atom_count + bond_count):
                    if len(lines[i]) < 9:  # Minimum bağ satırı uzunluğu
                        return False, f"Invalid MOL file: Invalid bond line {i-3-atom_count}"
                        
            except ValueError:
                return False, "Invalid MOL file: Invalid number format in counts line"
                
            return True, "MOL file is valid"
            
        except Exception as e:
            return False, f"Error validating MOL file: {str(e)}"

    def load_molecule(self, filepath):
        """Molekül yükleme ana fonksiyonu"""
        try:
            # 1. Dosya validasyonu
            is_valid, message = self.validate_mol_file(filepath)
            if not is_valid:
                raise ValueError(message)
                
            # 2. Molekül okuma
            mol = self.read_molecule(filepath)
            if mol is None:
                raise ValueError("Failed to read molecule")
                
            # 3. Molekül özellikleri kontrolü
            if mol.GetNumAtoms() == 0:
                raise ValueError("Molecule contains no atoms")
                
            # 4. Atom tipleri kontrolü
            valid_atoms = set(['H', 'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I'])
            for atom in mol.GetAtoms():
                if atom.GetSymbol() not in valid_atoms:
                    print(f"Warning: Unusual atom type found: {atom.GetSymbol()}")
                    
            # 5. 2D koordinat kontrolü
            if not mol.GetConformer().Is3D():
                AllChem.Compute2DCoords(mol)
                
            return mol
            
        except Exception as e:
            print(f"Error loading molecule: {str(e)}")
            return None

    def _read_cdxml(self, filepath):
        """CDXML dosyasını oku ve RDKit Molekülüne çevir"""
        try:
            tree = ET.parse(filepath)
            root = tree.getroot()
            
            # CDXML'den Mol blok çıkarma (basitleştirilmiş)
            ns = {'cdx': 'http://www.cambridgesoft.com/xml/cdxml.dtd'}
            molblock = ""
            
            for fragment in root.findall('.//cdx:fragment', ns):
                for mol in fragment.findall('.//cdx:mol', ns):
                    molblock = ET.tostring(mol, encoding='unicode')
                    break
                if molblock: break
                
            if molblock:
                # CDXML'den gelen veriyi temizle
                molblock = molblock.replace('<?xml version="1.0" ?>', '')
                mol = Chem.MolFromMolBlock(molblock, sanitize=False)
                if mol:
                    try:
                        Chem.SanitizeMol(mol)
                        AllChem.Compute2DCoords(mol)
                        return mol
                    except:
                        mol.UpdatePropertyCache(strict=False)
                        return mol
        except Exception as e:
            print(f"CDXML parsing error: {str(e)}")
        return None

    def find_significant_peaks(self, spectrum_data, threshold=0.3):
        """Spektrumdaki önemli pikleri tespit et"""
        peaks = []
        wavenumbers = np.array(spectrum_data['wavenumbers'])
        absorbance = np.array(spectrum_data['absorbance'])
        
        # Scipy'nin peak_finding algoritmasını kullan
        from scipy.signal import find_peaks
        peak_indices, _ = find_peaks(absorbance, height=threshold, distance=20)
        
        for idx in peak_indices:
            peaks.append(wavenumbers[idx])
        
        return peaks

    def validate_spectrum_data(self, spectrum_data):
        """Validate spectrum data before plotting"""
        if not isinstance(spectrum_data, dict):
            return False
        
        required_keys = ['wavenumbers', 'transmittance', 'absorbance', 
                        'functional_groups', 'peak_assignments']
        
        if not all(key in spectrum_data for key in required_keys):
            return False
        
        # Check data lengths
        if len(spectrum_data['wavenumbers']) != len(spectrum_data['transmittance']):
            return False
        
        # Check value ranges
        if not (np.all(np.array(spectrum_data['transmittance']) >= 0) and 
                np.all(np.array(spectrum_data['transmittance']) <= 100)):
            return False
        
        return True
    
class FTIRApp:
    def __init__(self, root):
        self.root = root
        self.predictor = FTIRPredictor()
        self.current_mol = None
        self.current_spectrum = None
        self.setup_gui()
        self.peak_text = None
        self.highlighted_atoms = []
        self.original_mol_img = None
    
    def setup_gui(self):
        """GUI'yi kur"""
        self.root.title("FTIR Spectrum Predictor")
        self.root.geometry("1200x800")
        
        # Ana çerçeve
        main_frame = tk.Frame(self.root)
        main_frame.pack(fill='both', expand=True, padx=10, pady=10)
        
        # Sol panel - Giriş ve bilgiler
        left_panel = tk.Frame(main_frame, width=300, bg='#f0f0f0')
        left_panel.pack(side='left', fill='y', padx=(0, 10))
        left_panel.pack_propagate(False)
        
        # Sağ panel - Spektrum ve molekül görüntüleyici
        right_panel = tk.Frame(main_frame)
        right_panel.pack(side='right', fill='both', expand=True)
        
        # Molekül görüntüleyici
        mol_frame = tk.Frame(right_panel, height=400)
        mol_frame.pack(fill='x', pady=(0, 10))
        
        self.mol_fig, self.mol_ax = plt.subplots(figsize=(8, 4))
        self.mol_canvas = FigureCanvasTkAgg(self.mol_fig, master=mol_frame)
        self.mol_canvas.get_tk_widget().pack(fill='both', expand=True)
        self.mol_ax.axis('off')
        self.mol_ax.text(0.5, 0.5, 'Molecule will appear here',
                        ha='center', va='center', fontsize=12, color='gray')
        
        # Spektrum görüntüleyici
        spec_frame = tk.Frame(right_panel)
        spec_frame.pack(fill='both', expand=True)
        
        self.spec_fig, self.spec_ax = plt.subplots(figsize=(8, 4))
        self.spec_canvas = FigureCanvasTkAgg(self.spec_fig, master=spec_frame)
        self.spec_canvas.get_tk_widget().pack(fill='both', expand=True)
        self.setup_spectrum_plot()
        
        # Giriş kontrolleri
        input_frame = tk.Frame(left_panel, bg='#f0f0f0')
        input_frame.pack(fill='x', pady=(0, 10))
        
        tk.Label(input_frame, text="MOL File:", bg='#f0f0f0').pack(anchor='w')
        
        self.file_entry = tk.Entry(input_frame)
        self.file_entry.pack(fill='x', pady=(0, 5))
        
        browse_btn = tk.Button(input_frame, text="Browse", command=self.browse_file)
        browse_btn.pack(fill='x')
        
        predict_btn = tk.Button(input_frame, text="Predict Spectrum", 
                              command=self.predict_spectrum, bg='#4CAF50', fg='white')
        predict_btn.pack(fill='x', pady=(10, 0))
        
        # Molekül bilgileri
        info_frame = tk.Frame(left_panel, bg='#f0f0f0')
        info_frame.pack(fill='x', pady=(10, 0))
        
        tk.Label(info_frame, text="Molecule Info:", bg='#f0f0f0', 
                font=('Arial', 10, 'bold')).pack(anchor='w')
        
        self.mol_info = tk.Label(info_frame, text="No molecule loaded", 
                               bg='#f0f0f0', wraplength=280, justify='left')
        self.mol_info.pack(fill='x', pady=(5, 0))
        
        # Fonksiyonel gruplar
        groups_frame = tk.Frame(left_panel, bg='#f0f0f0')
        groups_frame.pack(fill='both', expand=True, pady=(10, 0))
        
        tk.Label(groups_frame, text="Functional Groups:", bg='#f0f0f0',
               font=('Arial', 10, 'bold')).pack(anchor='w')
        
        self.groups_text = tk.Text(groups_frame, height=10, wrap='word', 
                                 bg='white', font=('Arial', 9))
        self.groups_text.pack(fill='both', expand=True, pady=(5, 0))
        
        # İmleç etkileşimi
        self.cursor = Cursor(self.spec_ax, useblit=True, color='red', linewidth=1)
        self.annotation = self.spec_ax.annotate("", xy=(0, 0), xytext=(20, 20),
                                              textcoords="offset points",
                                              bbox=dict(boxstyle="round", fc="w"),
                                              arrowprops=dict(arrowstyle="->"))
        self.annotation.set_visible(False)
        
        self.spec_canvas.mpl_connect("motion_notify_event", self.on_mouse_move)
        self.spec_canvas.mpl_connect("axes_leave_event", self.on_leave_axes)
        
    def browse_file(self):
        """Modern file browser with better UX from ftir.py"""
        filetypes = [("MOL files", "*.mol"), ("SDF files", "*.sdf"), ("All files", "*.*")]
        filename = filedialog.askopenfilename(
            title="Select MOL file",
            filetypes=filetypes,
            initialdir=os.getcwd()
        )
        if filename:
            self.file_entry.delete(0, tk.END)
            self.file_entry.insert(0, filename)
            self.display_molecule(filename)

    def display_molecule(self, filepath):
        """Molekülü göster"""
        try:
            # Molekülü yükle
            mol = self.predictor.read_molecule(filepath)
            
            if mol is None:
                # Alternatif okuma yöntemleri deneyelim
                with open(filepath, 'r') as f:
                    content = f.read()
                    mol = Chem.MolFromMolBlock(content)
                
                if mol is None:
                    raise ValueError("RDKit could not parse the molecule file")
            
            self.current_mol = mol
            
            # 2D koordinatları hesapla
            AllChem.Compute2DCoords(mol)
            
            # Molekül bilgilerini güncelle
            formula = rdMolDescriptors.CalcMolFormula(mol)
            weight = rdMolDescriptors.CalcExactMolWt(mol)
            info_text = f"Formula: {formula}\nMolecular Weight: {weight:.2f}"
            self.mol_info.config(text=info_text)
            
            # Molekülü çiz ve orijinal görüntüyü sakla
            self.mol_ax.clear()
            img = Draw.MolToImage(mol, size=(400, 300))
            self.original_mol_img = img  # Store the original image
            self.mol_ax.imshow(img)
            self.mol_ax.axis('off')
            self.mol_canvas.draw()
            
            # Fonksiyonel grupları tespit et
            groups = self.predictor.identify_functional_groups(mol)
            groups_text = "Detected groups:\n"
            for group, present in groups.items():
                if present:
                    groups_text += f"• {group.replace('_', ' ').title()}\n"
            
            self.groups_text.delete(1.0, tk.END)
            self.groups_text.insert(tk.END, groups_text)
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to display molecule:\n{str(e)}")
            print(f"Error displaying molecule: {str(e)}")
            # Reset display
            self.mol_ax.clear()
            self.mol_ax.text(0.5, 0.5, 'Failed to load molecule',
                           ha='center', va='center', fontsize=12, color='red')
            self.mol_ax.axis('off')
            self.mol_canvas.draw()
            self.mol_info.config(text="No molecule loaded")
            self.groups_text.delete(1.0, tk.END)
    
    def on_mouse_move(self, event):
        """Geliştirilmiş fare hareketi işleyicisi"""
        if event.inaxes != self.spec_ax or not hasattr(self, 'peak_assignments'):
            return
        
        x = event.xdata
        y = event.ydata
        
        # En yakın pik atamasını bul
        closest_peak = None
        min_dist = float('inf')
        
        for peak in self.peak_assignments:
            dist = abs(x - peak['center'])
            if dist < min_dist and dist < 50:  # 50 cm⁻¹ içinde
                min_dist = dist
                closest_peak = peak
        
        if closest_peak:
            # Açıklamayı güncelle
            text = (f"{closest_peak['type'].title()}\n"
                   f"{closest_peak['group'].replace('_', ' ').title()}\n"
                   f"{closest_peak['range'][0]:.0f}-{closest_peak['range'][1]:.0f} cm⁻¹\n"
                   f"Appearance: {closest_peak['appearance']}\n"
                   f"Comments: {closest_peak['comments']}")
            
            self.annotation.xy = (x, y)
            self.annotation.set_text(text)
            self.annotation.set_visible(True)
            
            # İlgili bağları vurgula
            self.highlight_bonds(closest_peak)
        else:
            self.annotation.set_visible(False)
            # Vurgulamaları temizle ve orijinal molekülü göster
            if self.original_mol_img is not None:
                self.mol_ax.clear()
                self.mol_ax.imshow(self.original_mol_img)
                self.mol_ax.axis('off')
                self.mol_canvas.draw()
        
        self.spec_canvas.draw_idle()
    
    def highlight_bonds(self, peak_info):
        # Renk tanımları
        COLOR_SCHEME = {
            'O-H': (0.2, 0.4, 0.8),    # Mavi
            'N-H': (0.4, 0.8, 0.4),     # Yeşil
            'C-H': {
                'alkane': (0.8, 0.8, 0.2),   # Sarı
                'aromatic': (0.8, 0.6, 0.2), # Turuncu
                'alkene': (0.8, 0.4, 0.6),   # Pembe
                'alkyne': (0.6, 0.2, 0.8),   # Mor
                'aldehyde': (0.8, 0.5, 0.5)  # Açık kırmızı
            },
            'C=O': (0.8, 0.2, 0.2),    # Kırmızı
            'C=C': (0.8, 0.4, 0.4),     # Açık kırmızı
            'C≡C': (0.6, 0.2, 0.6),     # Mor
            'C-N': (0.4, 0.6, 0.2),     # Zeytin yeşili
            'N-H bending': (0.2, 0.8, 0.4)  # Parlak yeşil
        }

        """Molekül üzerinde ilgili bağları ve fonksiyonel grupları vurgula"""
        if self.current_mol is None or peak_info is None:
            return
            
        mol = self.current_mol
        
        # Vurgulanacak atom ve bağları belirle
        highlight_atoms = set()
        highlight_bonds = set()
        
        # Pik bilgisine göre vurgulanacak kısımları belirle
        group = peak_info.get('group', '')
        bond_type = peak_info.get('type', '')
        
        # Fonksiyonel grup ve bağ tipine göre SMARTS desenleri
        smarts_patterns = []
        
        
        # N-H bending için özel vurgulama (1650-1580 cm⁻¹ aralığı)
        if 'N-H bending' in bond_type or ('amine' in group and 1580 <= peak_info.get('center', 0) <= 1650):
            # Primary amine N-H (yeşil)
            smarts_patterns.append('[NX3H2]')
            # Secondary amine N-H (açık yeşil)
            smarts_patterns.append('[NX3H1]')
            
        # Genel fonksiyonel grup desenleri
        if 'alcohol' in group or 'O-H' in bond_type:
            smarts_patterns.append('[OH]')
        if 'amine' in group or 'N-H' in bond_type:
            # Primary amine
            smarts_patterns.append('[NX3H2]')
            # Secondary amine
            smarts_patterns.append('[NX3H1]([#6])[#6]')
            # Tertiary amine
            smarts_patterns.append('[NX3H0]([#6])([#6])[#6]')
        # N-H bending için özel vurgulama
        if 'N-H bending' in peak_info.get('type', ''):
            # Primary amine N-H bending (yeşil)
            smarts_patterns.append('[NX3H2]')
            # Secondary amine N-H bending (açık yeşil)
            smarts_patterns.append('[NX3H1]')
            
            atom_colors.update({i: (0.4, 0.8, 0.4) for i in highlight_atoms})  # Yeşil
            bond_colors.update({i: (0.4, 0.8, 0.4) for i in highlight_bonds})  # Yeşil
            
        if 'carboxylic_acid' in group or 'C=O' in bond_type:
            smarts_patterns.append('C(=O)[OH]')
        if 'ester' in group:
            smarts_patterns.append('C(=O)[OX2]')
        if 'aldehyde' in group:
            smarts_patterns.append('[CX3H1](=O)')
        if 'ketone' in group:
            smarts_patterns.append('[CX3](=O)[#6]')
        if 'ether' in group or 'C-O-C' in bond_type:
            smarts_patterns.append('[OD2]([#6])[#6]')
        if 'alkene' in group or 'C=C' in bond_type:
            smarts_patterns.append('[CX3]=[CX3]')
        if 'alkyne' in group or 'C#C' in bond_type:
            smarts_patterns.append('[CX2]#[CX2]')
        if 'aromatic' in group:
            smarts_patterns.append('c1ccccc1')
        if 'amide' in group:
            smarts_patterns.append('C(=O)[NX3]')
        if 'nitrile' in group or 'C#N' in bond_type:
            smarts_patterns.append('[NX1]#[CX2]')
        if 'thiol' in group or 'S-H' in bond_type:
            smarts_patterns.append('[SH]')
        if 'halide' in group:
            smarts_patterns.append('[F,Cl,Br,I]')
        if 'nitro' in group or 'N=O' in bond_type:
            smarts_patterns.append('[NX3](=O)=O')
        if 'sulfoxide' in group or 'S=O' in bond_type:
            smarts_patterns.append('[SX3](=O)')
        if 'sulfone' in group:
            smarts_patterns.append('[SX4](=O)(=O)')
        if 'anhydride' in group:
            smarts_patterns.append('C(=O)OC(=O)')
        if 'imide' in group:
            smarts_patterns.append('C(=O)[NX3]C(=O)')
        
        # Özel bağ tipleri
        if 'O=C=O' in bond_type:
            smarts_patterns.append('O=C=O')
        if 'N=C=O' in bond_type:
            smarts_patterns.append('N=C=O')
        if 'N=C=S' in bond_type:
            smarts_patterns.append('N=C=S')
        if 'N=N=N' in bond_type:
            smarts_patterns.append('N=N=N')
        if 'C=C=C' in bond_type:
            smarts_patterns.append('C=C=C')
        
        # Tüm desenleri kontrol et
        for pattern in smarts_patterns:
            patt = Chem.MolFromSmarts(pattern)
            if patt:
                matches = mol.GetSubstructMatches(patt)
                for match in matches:
                    highlight_atoms.update(match)
        
        # Vurgulanan atomlar arasındaki bağları bul
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if begin_idx in highlight_atoms and end_idx in highlight_atoms:
                highlight_bonds.add(bond.GetIdx())
        
        # Molekülü vurgulanmış şekilde çiz
        highlight_atoms = list(highlight_atoms)
        highlight_bonds = list(highlight_bonds)
        
        # Renkleri belirle (N-H bending için özel renk)
        atom_colors = {}
        bond_colors = {}
        
        # N-H bending için yeşil renk
        if 'N-H bending' in bond_type or ('amine' in group and 1580 <= peak_info.get('center', 0) <= 1650):
            for i in highlight_atoms:
                atom_colors[i] = COLOR_SCHEME['N-H bending']
            for i in highlight_bonds:
                bond_colors[i] = COLOR_SCHEME['N-H bending']
        else:
            # Diğer gruplar için varsayılan renkler
            for i in highlight_atoms:
                atom_colors[i] = (0.8, 0.4, 0.4)  # Kırmızımsı
            for i in highlight_bonds:
                bond_colors[i] = (0.8, 0.4, 0.4)  # Kırmızımsı
        
        # Molekülü çiz
        img = Draw.MolToImage(mol, size=(400, 300), 
                             highlightAtoms=highlight_atoms,
                             highlightBonds=highlight_bonds,
                             highlightAtomColors=atom_colors,
                             highlightBondColors=bond_colors)
        
        # Görüntüyü güncelle
        self.mol_ax.clear()
        self.mol_ax.imshow(img)
        self.mol_ax.axis('off')
        self.mol_canvas.draw()
      
    def on_leave_axes(self, event):
        """Eksenlerden çıkınca açıklamayı gizle"""
        if self.annotation.get_visible():
            self.annotation.set_visible(False)
            self.spec_canvas.draw_idle()

    def setup_spectrum_plot(self):
        """Standart IR grafiği (0% altta, 100% üstte)"""
        self.spec_ax.clear()
        self.spec_ax.set_xlabel('Wavenumber (cm⁻¹)')
        self.spec_ax.set_ylabel('Transmittance (%)')
        self.spec_ax.set_xlim(4000, 400)
        self.spec_ax.set_ylim(0, 100)
        self.spec_ax.invert_yaxis()  # Standart IR görünümü
        self.spec_ax.grid(True, linestyle='--', alpha=0.5)
        self.spec_canvas.draw()
    
    def predict_spectrum(self):
        """Improved spectrum prediction method"""
        filepath = self.file_entry.get()
        try:
            # Load molecule
            mol = self.predictor.read_molecule(filepath)
            if mol is None:
                raise ValueError("Molecule could not be parsed")
            
            # Generate spectrum
            result = self.predictor.generate_spectrum(mol)
            if not result:
                raise ValueError("Spectrum generation failed")
            
            # Store peak assignments for hover functionality
            self.peak_assignments = result['peak_assignments']
            
            # Clear and plot new spectrum
            self.spec_ax.clear()
            wavenumbers = np.array(result['wavenumbers'])
            transmittance = np.array(result['transmittance'])
            
            # Plot with proper axis orientation
            self.spec_ax.plot(wavenumbers, transmittance, 'b-', linewidth=1)
            
            # Set up axes
            self.spec_ax.set_xlabel('Wavenumber (cm⁻¹)')
            self.spec_ax.set_ylabel('Transmittance (%)')
            self.spec_ax.set_xlim(4000, 400)  # Reverse x-axis for IR convention
            self.spec_ax.set_ylim(0, 100)
            
            # Major and minor ticks
            self.spec_ax.xaxis.set_major_locator(plt.MultipleLocator(500))
            self.spec_ax.xaxis.set_minor_locator(plt.MultipleLocator(100))
            self.spec_ax.xaxis.set_tick_params(which='minor', labelbottom=False)
            
            self.spec_ax.grid(True, alpha=0.3, which='both')
            
            # Update the plot
            self.spec_canvas.draw()
            
            # Update molecule info
            formula = rdMolDescriptors.CalcMolFormula(mol)
            weight = rdMolDescriptors.CalcExactMolWt(mol)
            
            # Update functional groups text
            groups_text = "Detected functional groups:\n"
            for group, present in result['functional_groups'].items():
                if present:
                    groups_text += f"• {group.replace('_', ' ').title()}\n"
            
            self.groups_text.delete(1.0, tk.END)
            self.groups_text.insert(tk.END, groups_text)
            
            # Update molecule info
            self.mol_info.config(text=f"Formula: {formula}\nMolecular Weight: {weight:.2f} g/mol")
            
            # Önemli pikleri tespit et (intensity threshold'a göre)
            significant_peaks = []
            for peak in self.peak_assignments:
                if peak['height'] > 0.3:  # Önemli pik eşiği
                    significant_peaks.append(peak['center'])
        
        except Exception as e:
            messagebox.showerror("Error", f"Failed to generate spectrum:\n{str(e)}")
            print(f"Error generating spectrum: {str(e)}")
            self.setup_spectrum_plot()  # Reset to default view         

def main():
    root = tk.Tk()
    app = FTIRApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()