class Element:
    """单个化学元素类"""
    
    def __init__(self, symbol, name, atomic_number, atomic_mass, group, period, block):
        self.symbol = symbol
        self.name = name
        self.atomic_number = atomic_number
        self.atomic_mass = atomic_mass
        self.group = group
        self.period = period
        self.block = block
    
    def __str__(self):
        return f"{self.name} ({self.symbol}), atomic_number: {self.atomic_number}, atomic_mass: {self.atomic_mass}"
    
    def __repr__(self):
        return f"Element('{self.symbol}', '{self.name}', {self.atomic_number}, {self.atomic_mass}, {self.group}, {self.period}, '{self.block}')"

class PeriodicTable:
    """元素周期表类"""
    
    def __init__(self):
        # 直接在初始化时定义所有元素
        self.elements = {
            "H": Element("H", "Hydrogen", 1, 1.008, 1, 1, "s"),
            "He": Element("He", "Helium", 2, 4.003, 18, 1, "s"),
            "Li": Element("Li", "Lithium", 3, 6.941, 1, 2, "s"),
            "Be": Element("Be", "Beryllium", 4, 9.012, 2, 2, "s"),
            "B": Element("B", "Boron", 5, 10.811, 13, 2, "p"),
            "C": Element("C", "Carbon", 6, 12.011, 14, 2, "p"),
            "N": Element("N", "Nitrogen", 7, 14.007, 15, 2, "p"),
            "O": Element("O", "Oxygen", 8, 15.999, 16, 2, "p"),
            "F": Element("F", "Fluorine", 9, 18.998, 17, 2, "p"),
            "Ne": Element("Ne", "Neon", 10, 20.180, 18, 2, "p"),
            "Na": Element("Na", "Sodium", 11, 22.990, 1, 3, "s"),
            "Mg": Element("Mg", "Magnesium", 12, 24.305, 2, 3, "s"),
            "Al": Element("Al", "Aluminum", 13, 26.982, 13, 3, "p"),
            "Si": Element("Si", "Silicon", 14, 28.085, 14, 3, "p"),
            "P": Element("P", "Phosphorus", 15, 30.974, 15, 3, "p"),
            "S": Element("S", "Sulfur", 16, 32.065, 16, 3, "p"),
            "Cl": Element("Cl", "Chlorine", 17, 35.453, 17, 3, "p"),
            "Ar": Element("Ar", "Argon", 18, 39.948, 18, 3, "p"),
            "K": Element("K", "Potassium", 19, 39.098, 1, 4, "s"),
            "Ca": Element("Ca", "Calcium", 20, 40.078, 2, 4, "s"),
            "Sc": Element("Sc", "Scandium", 21, 44.956, 3, 4, "d"),
            "Ti": Element("Ti", "Titanium", 22, 47.867, 4, 4, "d"),
            "V": Element("V", "Vanadium", 23, 50.942, 5, 4, "d"),
            "Cr": Element("Cr", "Chromium", 24, 51.996, 6, 4, "d"),
            "Mn": Element("Mn", "Manganese", 25, 54.938, 7, 4, "d"),
            "Fe": Element("Fe", "Iron", 26, 55.845, 8, 4, "d"),
            "Co": Element("Co", "Cobalt", 27, 58.933, 9, 4, "d"),
            "Ni": Element("Ni", "Nickel", 28, 58.693, 10, 4, "d"),
            "Cu": Element("Cu", "Copper", 29, 63.546, 11, 4, "d"),
            "Zn": Element("Zn", "Zinc", 30, 65.38, 12, 4, "d"),
            "Ga": Element("Ga", "Gallium", 31, 69.723, 13, 4, "p"),
            "Ge": Element("Ge", "Germanium", 32, 72.630, 14, 4, "p"),
            "As": Element("As", "Arsenic", 33, 74.922, 15, 4, "p"),
            "Se": Element("Se", "Selenium", 34, 78.971, 16, 4, "p"),
            "Br": Element("Br", "Bromine", 35, 79.904, 17, 4, "p"),
            "Kr": Element("Kr", "Krypton", 36, 83.798, 18, 4, "p"),
            "Rb": Element("Rb", "Rubidium", 37, 85.468, 1, 5, "s"),
            "Sr": Element("Sr", "Strontium", 38, 87.62, 2, 5, "s"),
            "Y": Element("Y", "Yttrium", 39, 88.906, 3, 5, "d"),
            "Zr": Element("Zr", "Zirconium", 40, 91.224, 4, 5, "d"),
            "Nb": Element("Nb", "Niobium", 41, 92.906, 5, 5, "d"),
            "Mo": Element("Mo", "Molybdenum", 42, 95.95, 6, 5, "d"),
            "Tc": Element("Tc", "Technetium", 43, 98, 7, 5, "d"),
            "Ru": Element("Ru", "Ruthenium", 44, 101.07, 8, 5, "d"),
            "Rh": Element("Rh", "Rhodium", 45, 102.906, 9, 5, "d"),
            "Pd": Element("Pd", "Palladium", 46, 106.42, 10, 5, "d"),
            "Ag": Element("Ag", "Silver", 47, 107.868, 11, 5, "d"),
            "Cd": Element("Cd", "Cadmium", 48, 112.414, 12, 5, "d"),
            "In": Element("In", "Indium", 49, 114.818, 13, 5, "p"),
            "Sn": Element("Sn", "Tin", 50, 118.710, 14, 5, "p"),
            "Sb": Element("Sb", "Antimony", 51, 121.760, 15, 5, "p"),
            "Te": Element("Te", "Tellurium", 52, 127.60, 16, 5, "p"),
            "I": Element("I", "Iodine", 53, 126.904, 17, 5, "p"),
            "Xe": Element("Xe", "Xenon", 54, 131.293, 18, 5, "p"),
            "Cs": Element("Cs", "Cesium", 55, 132.905, 1, 6, "s"),
            "Ba": Element("Ba", "Barium", 56, 137.327, 2, 6, "s"),
            "La": Element("La", "Lanthanum", 57, 138.905, 3, 6, "d"),
            "Ce": Element("Ce", "Cerium", 58, 140.116, 3, 6, "f"),
            "Pr": Element("Pr", "Praseodymium", 59, 140.908, 3, 6, "f"),
            "Nd": Element("Nd", "Neodymium", 60, 144.242, 3, 6, "f"),
            "Pm": Element("Pm", "Promethium", 61, 145, 3, 6, "f"),
            "Sm": Element("Sm", "Samarium", 62, 150.36, 3, 6, "f"),
            "Eu": Element("Eu", "Europium", 63, 151.964, 3, 6, "f"),
            "Gd": Element("Gd", "Gadolinium", 64, 157.25, 3, 6, "f"),
            "Tb": Element("Tb", "Terbium", 65, 158.925, 3, 6, "f"),
            "Dy": Element("Dy", "Dysprosium", 66, 162.500, 3, 6, "f"),
            "Ho": Element("Ho", "Holmium", 67, 164.930, 3, 6, "f"),
            "Er": Element("Er", "Erbium", 68, 167.259, 3, 6, "f"),
            "Tm": Element("Tm", "Thulium", 69, 168.934, 3, 6, "f"),
            "Yb": Element("Yb", "Ytterbium", 70, 173.054, 3, 6, "f"),
            "Lu": Element("Lu", "Lutetium", 71, 174.967, 3, 6, "d"),
            "Hf": Element("Hf", "Hafnium", 72, 178.49, 4, 6, "d"),
            "Ta": Element("Ta", "Tantalum", 73, 180.948, 5, 6, "d"),
            "W": Element("W", "Tungsten", 74, 183.84, 6, 6, "d"),
            "Re": Element("Re", "Rhenium", 75, 186.207, 7, 6, "d"),
            "Os": Element("Os", "Osmium", 76, 190.23, 8, 6, "d"),
            "Ir": Element("Ir", "Iridium", 77, 192.217, 9, 6, "d"),
            "Pt": Element("Pt", "Platinum", 78, 195.084, 10, 6, "d"),
            "Au": Element("Au", "Gold", 79, 196.967, 11, 6, "d"),
            "Hg": Element("Hg", "Mercury", 80, 200.592, 12, 6, "d"),
            "Tl": Element("Tl", "Thallium", 81, 204.383, 13, 6, "p"),
            "Pb": Element("Pb", "Lead", 82, 207.2, 14, 6, "p"),
            "Bi": Element("Bi", "Bismuth", 83, 208.980, 15, 6, "p"),
            "Po": Element("Po", "Polonium", 84, 209, 16, 6, "p"),
            "At": Element("At", "Astatine", 85, 210, 17, 6, "p"),
            "Rn": Element("Rn", "Radon", 86, 222, 18, 6, "p"),
            "Fr": Element("Fr", "Francium", 87, 223, 1, 7, "s"),
            "Ra": Element("Ra", "Radium", 88, 226, 2, 7, "s"),
            "Ac": Element("Ac", "Actinium", 89, 227, 3, 7, "d"),
            "Th": Element("Th", "Thorium", 90, 232.038, 3, 7, "f"),
            "Pa": Element("Pa", "Protactinium", 91, 231.036, 3, 7, "f"),
            "U": Element("U", "Uranium", 92, 238.029, 3, 7, "f"),
            "Np": Element("Np", "Neptunium", 93, 237, 3, 7, "f"),
            "Pu": Element("Pu", "Plutonium", 94, 244, 3, 7, "f"),
            "Am": Element("Am", "Americium", 95, 243, 3, 7, "f"),
            "Cm": Element("Cm", "Curium", 96, 247, 3, 7, "f"),
            "Bk": Element("Bk", "Berkelium", 97, 247, 3, 7, "f"),
            "Cf": Element("Cf", "Californium", 98, 251, 3, 7, "f"),
            "Es": Element("Es", "Einsteinium", 99, 252, 3, 7, "f"),
            "Fm": Element("Fm", "Fermium", 100, 257, 3, 7, "f"),
            "Md": Element("Md", "Mendelevium", 101, 258, 3, 7, "f"),
            "No": Element("No", "Nobelium", 102, 259, 3, 7, "f"),
            "Lr": Element("Lr", "Lawrencium", 103, 262, 3, 7, "d"),
            "Rf": Element("Rf", "Rutherfordium", 104, 267, 4, 7, "d"),
            "Db": Element("Db", "Dubnium", 105, 268, 5, 7, "d"),
            "Sg": Element("Sg", "Seaborgium", 106, 271, 6, 7, "d"),
            "Bh": Element("Bh", "Bohrium", 107, 272, 7, 7, "d"),
            "Hs": Element("Hs", "Hassium", 108, 270, 8, 7, "d"),
            "Mt": Element("Mt", "Meitnerium", 109, 276, 9, 7, "d"),
            "Ds": Element("Ds", "Darmstadtium", 110, 281, 10, 7, "d"),
            "Rg": Element("Rg", "Roentgenium", 111, 280, 11, 7, "d"),
            "Cn": Element("Cn", "Copernicium", 112, 285, 12, 7, "d"),
            "Nh": Element("Nh", "Nihonium", 113, 284, 13, 7, "p"),
            "Fl": Element("Fl", "Flerovium", 114, 289, 14, 7, "p"),
            "Mc": Element("Mc", "Moscovium", 115, 288, 15, 7, "p"),
            "Lv": Element("Lv", "Livermorium", 116, 293, 16, 7, "p"),
            "Ts": Element("Ts", "Tennessine", 117, 294, 17, 7, "p"),
            "Og": Element("Og", "Oganesson", 118, 294, 18, 7, "p")
        }
        
        # 初始化索引
        self._initialize_indices()
    
    def _initialize_indices(self):
        """初始化各种索引"""
        # 名称到元素的映射
        self.by_name = {element.name.lower(): element for element in self.elements.values()}
        
        # 原子序数到元素的映射
        self.by_number = {element.atomic_number: element for element in self.elements.values()}
    
    def get_element(self, identifier):
        """
        根据标识符获取元素
        
        参数:
        - identifier: 可以是元素符号、名称或原子序数
        
        返回:
        - 元素对象，如果找不到则返回None
        """
        # 根据不同类型的标识符查找元素
        if isinstance(identifier, int):
            # 原子序数
            return self.by_number.get(identifier)
        elif identifier in self.elements:
            # 元素符号
            return self.elements[identifier]
        elif identifier.lower() in self.by_name:
            # 元素名称
            return self.by_name[identifier.lower()]
        else:
            print(f"警告: 找不到元素 '{identifier}'")
            return None
    
    def get_atomic_mass(self, identifier):
        """获取指定元素的相对原子质量"""
        element = self.get_element(identifier)
        if element:
            return element.atomic_mass
        return None
    
    def get_atomic_number(self, identifier):
        """获取指定元素的原子序数"""
        element = self.get_element(identifier)
        if element:
            return element.atomic_number
        return None
    
    def get_element_name(self, identifier):
        """获取指定元素的全名"""
        element = self.get_element(identifier)
        if element:
            return element.name
        return None
    
    def get_all_elements(self):
        """获取所有元素的符号列表"""
        return list(self.elements.keys())
    
    def get_elements_by_group(self, group):
        """获取指定族的所有元素"""
        return [e for e in self.elements.values() if e.group == group]
    
    def get_elements_by_period(self, period):
        """获取指定周期的所有元素"""
        return [e for e in self.elements.values() if e.period == period]
    
    def get_elements_by_block(self, block):
        """获取指定区块的所有元素"""
        return [e for e in self.elements.values() if e.block == block]

# 创建周期表单例对象
periodic_table = PeriodicTable()

# 便捷函数
def get_element(identifier):
    """根据标识符获取元素"""
    return periodic_table.get_element(identifier)

def get_atomic_mass(element):
    """获取元素的原子质量"""
    return periodic_table.get_atomic_mass(element)

def get_atomic_number(element):
    """获取元素的原子序数"""
    return periodic_table.get_atomic_number(element)

def get_element_name(element):
    """获取元素的名称"""
    return periodic_table.get_element_name(element)

# 示例使用
if __name__ == "__main__":
    # 获取铜的信息
    cu = get_element("Cu")
    print(f"铜的信息: {cu}")
    
    # 获取铜的原子质量
    cu_mass = get_atomic_mass("Cu")
    print(f"铜的原子质量: {cu_mass}")
    
    # 使用原子序数获取元素
    fe = get_element(26)
    print(f"原子序数26是: {fe.name}")
    
    # 获取所有元素
    elements = periodic_table.get_all_elements()
    print(f"元素总数: {len(elements)}")
    print(f"前5个元素: {elements[:5]}")
    
    # 获取第1族元素
    group1 = periodic_table.get_elements_by_group(1)
    print(f"第1族元素: {[e.symbol for e in group1]}")
