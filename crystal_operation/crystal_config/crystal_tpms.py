import numpy as np

##### Schwarz-D 晶体 ######相关函数
def DSchwarz(Alloy, Range, Center = [1/2, 1/2, 1/2], Index = 1, eq = 2, angle = 0):
    Other = []
    X, Y, Z = Range[0], Range[1], Range[2]
    Structure = np.copy(Alloy).tolist()
    biasX, biasY, biasZ =  Center[0] - 1/2, Center[1] - 1/2, Center[2] - 1/2
    Center = [Center[0] * Range[0], Center[1] * Range[1], Center[2] * Range[2]]
    # 生成旋转矩阵
    rotation_mat = rotation_matrix(0, 0, angle)
    for s in Structure:
        sc =  rotation_mat @ np.array(s)
        if eq == 1:
            size = 2
            a, b, c = (sc[0] - biasX) * size, (sc[1] - biasY) * size, (sc[2]- biasZ) * size
            DSC = np.cos(a * np.pi / X) * np.cos(b * np.pi / Y) * np.cos(c * np.pi / Z) -\
              np.sin(a * np.pi / X) * np.sin(b * np.pi / Y) * np.sin(c * np.pi / Z)
        else:
            size = 1
            a, b, c = (sc[0] - biasX) * size, (sc[1] - biasY) * size, (sc[2]- biasZ) * size
            DSC = (np.sin(a * np.pi / X) * np.sin(b * np.pi / Y) * np.sin(c * np.pi / Z) +
                   np.sin(a * np.pi / X) * np.cos(b * np.pi / Y) * np.cos(c * np.pi / Z) +
                   np.cos(a * np.pi / X) * np.sin(b * np.pi / Y) * np.cos(c * np.pi / Z) +
                   np.cos(a * np.pi / X) * np.cos(b * np.pi / Y) * np.sin(c * np.pi / Z))
        if Index == 1:
            if DSC < 0:
                Other.append(s)
        elif Index ==2:
            if (DSC <= 0.1) and (DSC >= -0.1):
                Other.append(s)
        else:
            if DSC > 0:
                Other.append(s)
    return Other

##### Schwarz-G 晶体 ######相关函数
def GSchwarz(Alloy, Range, Center = [1/2, 1/2, 1/2], Index = 1):
    Other = []
    X, Y, Z = Range[0], Range[1], Range[2]
    Structure = np.copy(Alloy).tolist()
    biasX, biasY, biasZ =  Center[0] - 1/2, Center[1] - 1/2, Center[2] - 1/2
    for s in Structure:
        a, b, c = s[0] - biasX, s[1] - biasY, s[2]- biasZ
        size = 2
        GSC = (np.sin(size * a*np.pi/X)*np.cos(size * b*np.pi/Y) +
               np.sin(size * c*np.pi/Z)*np.cos(size * a*np.pi/X) +
               np.sin(size * b*np.pi/Y)*np.cos(size * c*np.pi/Z))
        if Index == 1:
            if GSC <= 0:
                Other.append(s)
        elif Index ==2:
            if (GSC <= 0.2) and (GSC >= -0.2):
                Other.append(s)
        else:
            if GSC > 0:
                Other.append(s)
    return Other

##### Schwarz-p 晶体 ######相关函数
def PSchwarz(Alloy, Range, Center = [1/2, 1/2, 1/2], Index = 1):
    Other = []
    X, Y, Z = Range[0], Range[1], Range[2]
    Structure = np.copy(Alloy).tolist()
    biasX, biasY, biasZ = Center[0] - 1/2, Center[1] - 1/2, Center[2] - 1/2
    for s in Structure:
        size = 2
        a, b, c = s[0] - biasX, s[1] - biasY, s[2]- biasZ
        PSC = np.cos(size*a*np.pi/X) + np.cos(size*b*np.pi/Y) + np.cos(size*c*np.pi/Z)
        if Index == 1:
            if PSC <= 0:
                Other.append(s)
        elif Index == 2:
            if (PSC <= 0.2) and (PSC >= -0.2):
                Other.append(s)
        else:
            if PSC > 0:
                Other.append(s)
    return Other

##### Schwarz-IWP 晶体 ######相关函数
def IWPSchwarz(Alloy, Range, Center = [1/2, 1/2, 1/2], Index = 1):
    Other = []
    X, Y, Z = Range[0], Range[1], Range[2]
    Structure = np.copy(Alloy).tolist()
    biasX, biasY, biasZ = Center[0] - 1/2, Center[1] - 1/2, Center[2] - 1/2
    for s in Structure:
        size = 2
        a, b, c = (s[0] - biasX) * size, (s[1] - biasY) * size, (s[2] - biasZ) * size
        ISC = 2*(np.cos(a*np.pi/X) * np.cos(b*np.pi/Y) +
                np.cos(b*np.pi/Y) * np.cos(c*np.pi/Z) +
                np.cos(c*np.pi/Z) * np.cos(a*np.pi/X)) - \
              (np.cos(2*a*np.pi/X) + np.cos(2*b*np.pi/Y) + np.cos(2*c*np.pi/Z))
        if Index == 1:
            if ISC <= 0:
                Other.append(s)
        elif Index == 2:
            if (ISC <= 0.2) and (ISC >= -0.2):
                Other.append(s)
        else:
            if ISC > 0:
                Other.append(s)
    return Other

##### Schwarz D 空间分布晶体 ######相关函数
def Math_Schwarz_M1(Alloy, Range, Center = [0, 0, 0], Index = 1, angle = 0):
    delt = 0.0000
    Other = []
    Structure = np.copy(Alloy).tolist()
    Center = [Center[0] * Range[0], Center[1] * Range[1], Center[2] * Range[2]]
    # 生成旋转矩阵
    rotation_mat = rotation_matrix(0, 0, angle)
    for s in Structure:
        a, b, c = s[0], s[1], s[2]
        # [11-1]
        Top = [Center[0], Center[1], Center[2] + Range[2] / 8]
        Top = rotation_mat @ np.array(Top)
        Vector = [1, 1, -1]
        Vector = rotation_mat @ np.array(Vector)
        B1 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) >= delt
        Top = [Center[0] + Range[2], Center[1], Center[2]]
        Top = rotation_mat @ np.array(Top)
        B2 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= -delt
        # [-111]
        Top = [Center[0] + Range[2] / 8, Center[1], Center[2]]
        Top = rotation_mat @ np.array(Top)
        Vector = [-1, 1, 1]
        Vector = rotation_mat @ np.array(Vector)
        B3 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) >= delt
        Top = [Center[0], Center[1] + Range[2], Center[2]]
        Top = rotation_mat @ np.array(Top)
        B4 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= -delt
        # [1-11]
        Top = [Center[0], Center[1] + Range[2] / 8, Center[2]]
        Top = rotation_mat @ np.array(Top)
        Vector = [1, -1, 1]
        Vector = rotation_mat @ np.array(Vector)
        B5 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) >= delt
        Top = [Center[0] + Range[2], Center[1], Center[2]]
        Top = rotation_mat @ np.array(Top)
        B6 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= -delt
        # [111]
        Top = [Center[0], Center[1] + Range[2], Center[2] + Range[2] / 8]
        Top = rotation_mat @ np.array(Top)
        Vector = [1, 1, 1]
        Vector = rotation_mat @ np.array(Vector)
        B7 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= delt
        Top = [Center[0], Center[1], Center[2]]
        Top = rotation_mat @ np.array(Top)
        B8 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) >= -delt
        if (B1 and B2 and B3 and B4 and B5 and B6 and B7 and B8):
            Other.append(s)
    return Other

def Math_Schwarz_M2(Alloy, Range, Center = [0, 0, 0], Index = 1, angle = 0):
    delt = 0.0000
    Other = []
    Structure = np.copy(Alloy).tolist()
    Center = [Center[0] * Range[0], Center[1] * Range[1], Center[2] * Range[2]]
    # 生成旋转矩阵
    rotation_mat = rotation_matrix(0, 0, angle)
    for s in Structure:
        a, b, c = s[0], s[1], s[2]
        # [11-1]
        Top = [Center[0] + Range[0], Center[1] + Range[1], Center[2] + Range[2] / 8 * 7]
        Top = rotation_mat @ np.array(Top)
        Vector = [1, 1, -1]
        Vector = rotation_mat @ np.array(Vector)
        B1 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= delt
        Top = [Center[0] + Range[0], Center[1], Center[2] + Range[2]]
        Top = rotation_mat @ np.array(Top)
        B2 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) >= -delt
        # [-111]
        Top = [Center[0] + Range[0] / 8 * 7, Center[1] + Range[1], Center[2] + Range[2]]
        Top = rotation_mat @ np.array(Top)
        Vector = [-1, 1, 1]
        Vector = rotation_mat @ np.array(Vector)
        B3 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= delt
        Top = [Center[0] + Range[0], Center[1] + Range[1], Center[2]]
        Top = rotation_mat @ np.array(Top)
        B4 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) >= -delt
        # [1-11]
        Top = [Center[0] + Range[0], Center[1] + Range[1] / 8 * 7, Center[2] + Range[2]]
        Top = rotation_mat @ np.array(Top)
        Vector = [1, -1, 1]
        Vector = rotation_mat @ np.array(Vector)
        B5 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= delt
        Top = [Center[0] + Range[0], Center[1] + Range[1], Center[2]]
        Top = rotation_mat @ np.array(Top)
        B6 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) >= -delt
        # [111]
        Top = [Center[0] + Range[0], Center[1] + Range[1] / 8 * 7, Center[2]]
        Top = rotation_mat @ np.array(Top)
        Vector = [1, 1, 1]
        Vector = rotation_mat @ np.array(Vector)
        B7 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) >= delt
        Top = [Center[0] + Range[0], Center[1] + Range[1], Center[2] + Range[2]]
        Top = rotation_mat @ np.array(Top)
        B8 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= -delt
        if (B1 and B2 and B3 and B4 and B5 and B6 and B7 and B8):
            Other.append(s)
    return Other

def Math_Schwarz_T1(Alloy, Range, Center = [0, 0, 0], Index = 1, angle = 0):
    delt = 0
    lamb = 0
    Other = []
    Structure = np.copy(Alloy).tolist()
    Center = [Center[0] * Range[0], Center[1] * Range[1], Center[2] * Range[2]]
    # 生成旋转矩阵
    rotation_mat = rotation_matrix(0, 0, angle)
    for s in Structure:
        a, b, c = s[0], s[1], s[2]
        # [111]
        Top = [Center[0] + Range[0], Center[1] + Range[1] / 8 * 7, Center[2]]
        Top = rotation_mat @ np.array(Top)
        Vector = [1, 1, 1]
        Vector = rotation_mat @ np.array(Vector)
        B1 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) < delt
        Top = [Center[0], Center[1] + Range[2], Center[2] + Range[2] / 8]
        Top = rotation_mat @ np.array(Top)
        B2 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) > lamb
        if (B1 and B2):
            Other.append(s)
    Other = Box_Slice(Other, Range, Sign=1)
    return DSchwarz(Other, Range, Index=3, eq=2, Center=Center)

def Math_Schwarz_T2(Alloy, Range, Center = [0, 0, 0], angle = 0):
    delt = 0.0
    lamb = 0
    Other = []
    Structure = np.copy(Alloy).tolist()
    Center = [Center[0] * Range[0], Center[1] * Range[1], Center[2] * Range[2]]
    # 生成旋转矩阵
    rotation_mat = rotation_matrix(0, 0, angle)
    for s in Structure:
        a, b, c = s[0], s[1], s[2]
        # [11-1]
        Top = [Center[0] + Range[0], Center[1] + Range[1] / 8 * 7, Center[2] + Range[2]]
        Top = rotation_mat @ np.array(Top)
        Vector = [1, 1, -1]
        Vector = rotation_mat @ np.array(Vector)
        B1 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) < delt
        Top = [Center[0] + Range[2] / 8, Center[1], Center[2]]
        Top = rotation_mat @ np.array(Top)
        B2 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) > lamb
        if (B1 and B2):
            Other.append(s)
    Other = Box_Slice(Other, Range, Sign = 1)
    return DSchwarz(Other, Range, Index=3, eq=2, Center=Center, angle = 180)

def Math_Schwarz_T3(Alloy, Range, Center = [0, 0, 0], angle = 0):
    delt = 0
    lamb = 0.0
    Other = []
    Structure = np.copy(Alloy).tolist()
    Center = [Center[0] * Range[0], Center[1] * Range[1], Center[2] * Range[2]]
    # 生成旋转矩阵
    rotation_mat = rotation_matrix(0, 0, angle)
    for s in Structure:
        a, b, c = s[0], s[1], s[2]
        # [-111]
        Top = [Center[0] + Range[0] / 8 * 7, Center[1], Center[2] + Range[2]]
        Top = rotation_mat @ np.array(Top)
        Vector = [-1, 1, 1]
        Vector = rotation_mat @ np.array(Vector)
        B1 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) > lamb
        Top = [Center[0] + Range[1] / 8, Center[1], Center[2] + Range[2]]
        Top = rotation_mat @ np.array(Top)
        B2 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) < delt
        if (B1 and B2):
            Other.append(s)
    Other = Box_Slice(Other, Range, Sign = 1)
    return DSchwarz(Other, Range, Index=1, eq=2, Center=Center, angle = -90)

def Math_Schwarz_T4(Alloy, Range, Center = [0, 0, 0], angle = 0):
    delt = 0
    lamb = 0.0
    Other = []
    Structure = np.copy(Alloy).tolist()
    Center = [Center[0] * Range[0], Center[1] * Range[1], Center[2] * Range[2]]
    # 生成旋转矩阵
    rotation_mat = rotation_matrix(0, 0, angle)
    for s in Structure:
        a, b, c = s[0], s[1], s[2]
        # [-111]
        Top = [Center[0] + Range[0], Center[1] + Range[1] / 8 * 1, Center[2]]
        Top = rotation_mat @ np.array(Top)
        Vector = [1, -1, 1]
        Vector = rotation_mat @ np.array(Vector)
        B1 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) < delt
        Top = [Center[0] + Range[0], Center[1] + Range[1] / 8 * 7, Center[2]]
        Top = rotation_mat @ np.array(Top)
        B2 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) > lamb
        if (B1 and B2):
            Other.append(s)
    Other = Box_Slice(Other, Range, Sign = 1)
    return DSchwarz(Other, Range, Index=1, eq=2, Center=Center, angle = 90)

def Math_D_SCI(Range=None, a=3.61, M = 63.546, SeedModel=0):
    if Range is None:
        value = random.random()  # 生成0到1之间的随机浮点数
        print("Range Not Input, Random value is:", value)
        Range = [value, value, value]  
    Schwarz = {}
    Fcc = FCC()
    center = [1 / 2, 1 / 2, 1 / 2]
    # 用中心生成坐标点，设置晶相来建立晶体
    axis = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    Cube = Grow_Crystal(Range, Fcc, axis, Center = center)
    Schwarz[0] = Math_Schwarz_M1(Cube, Range)
    if SeedModel == 1: Cfg(Schwarz[0], Range, a, M, 'D-SCI-1')
        
    axis = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    Cube = Grow_Crystal(Range, Fcc, axis, Center = center)
    Schwarz[1] = Math_Schwarz_M2(Cube, Range)
    if SeedModel == 1: Cfg(Schwarz[1], Range, a, M, 'D_SCI-2')

    # 建立111Kelvin晶体
    axis = [[-1, 2, 2], [2, -1, 2], [2, 2, -1]]
    Cube = Grow_Crystal(Range, Fcc, axis, Center = center)
    Schwarz[2] = Math_Schwarz_T1(Cube, Range)
    if SeedModel == 1: Cfg(Schwarz[2], Range, a, M, 'D_SCI-3')

    # 建立11-1Kelvin晶体
    axis = [[-1, 2, -2], [2, -1, -2], [-2, -2, -1]]
    Cube = Grow_Crystal(Range, Fcc, axis, Center = center)
    Schwarz[3] = Math_Schwarz_T2(Cube, Range)
    Schwarz[3] = Move(Schwarz[3], [-1 / 2, -1 / 2, 0], Range)
    if SeedModel == 1: Cfg(Schwarz[3], Range, a, M, 'D_SCI-4')

    # 建立-111Kelvin晶体
    axis = [[-1, -2, -2], [-2, -1, 2], [-2, 2, -1]]
    Cube = Grow_Crystal(Range, Fcc, axis, Center = center)
    Schwarz[4] = Math_Schwarz_T3(Cube, Range)
    Schwarz[4] = Move(Schwarz[4], [0, -1 / 2, -1 / 2], Range)
    if SeedModel == 1: Cfg(Schwarz[4], Range, a, M, 'D_SCI-5')

    # 建立1-11Kelvin晶体
    axis = [[-1, -2, 2], [-2, -1, -2], [2, -2, -1]]
    Cube = Grow_Crystal(Range, Fcc, axis, Center = center)
    Schwarz[5] = Math_Schwarz_T4(Cube, Range)
    Schwarz[5] = Move(Schwarz[5], [-1 / 2, 0, -1 / 2], Range)
    if SeedModel == 1: Cfg(Schwarz[5], Range, a, M, 'D_SCI-6')

    All = Merge_Group(Schwarz)
    # Range = np.dot(Range)
    # All = Lt.Box_Slice(All, Range, Sign = 1)
    N = len(All)
    print("Num of Atoms:", N)
    Size = round(Range[0] * a / 2 * math.sqrt(2) / 10, 2)
    print("D-SC Ds:", Size)
    return All

def Math_D_SCI_NT(Range=None, a=3.61, M = 63.546, pitch=None, head=None, roll=None, SeedModel=0):
    if Range is None:
        value = random.random()  # 生成0到1之间的随机浮点数
        print("Range Not Input, Random value is:", value)
        Range = [value, value, value]  
    Fcc = FCC()
    center = [1 / 2, 1 / 2, 1 / 2]
    DSC = {}
    axis = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    Cube = Grow_Crystal(Range, Fcc, axis, Center = [1/2, 1/2, 1/2])
    DSC[0] = DSchwarz(Cube, Range, Index = 3, eq = 1, Center = center)
    axis = rotation_matrix(pitch, head, roll)
    # axis = rotation_matrix(30, 60, 27)
    # axis = [[6, -5, 0], [5, 6, 0], [0, 0, 1]]
    Cube = Grow_Crystal(Range, Fcc, axis, Center = [1/2, 1/2, 1/2])
    DSC[1] = DSchwarz(Cube, Range, Index = 1, eq = 1, Center = center)
    All = Merge_Group(DSC)
    All = Box_Slice(All, Range, Sign = 1)
    N = len(All)
    print("All:", N)
    print("Num of Atoms:", N)
    Size = round(Range[0] * a / 2 * math.sqrt(2) / 10, 2)
    print("D-SC Ds:", Size)
    return All