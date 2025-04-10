import numpy as np
import math
import Crystal_Operation as Lt
# 定义FCC坐标
FCC = Lt.FCC()
# 定义晶格常数
a = 3.61
# 定义摩尔质量
M = 63.546
# 定义XYZ轴单位长度
# X, Y, Z = 5, 5, 5
# X, Y, Z = X * 3, Y * 3, Z * 3
num = 10.9
X, Y, Z = num, num, num

Range = [X, Y, Z]
Center = [Range[0] / 2, Range[1] / 2, Range[2] / 2]

GSC = {}
# 用中心生成坐标点，设置晶相来建立晶体
# 建立中心Kelvin晶体
axis = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
# axis = [[-1, 2, -2], [2, -1, -2], [-2, -2, -1]]
# Cube = Lt.Grow_Crystal(Range, FCC, axis, Center = [1/2, 1/2, 1/2])
center =  [0, 0, 0]
Cube = Lt.Grow_Crystal(Range, FCC, axis, Center = [1/2, 1/2, 1/2])
Cube = Lt.Box_Slice(Cube, Range, Sign = 1)
Lt.Cfg(Cube, Range, a, M, 'CB')

GSC = Lt.GSchwarz(Cube, Range, Index = 2, Center = center)
Lt.Cfg(GSC, Range, a, M, 'GSC')

DSC = Lt.DSchwarz(Cube, Range, Index = 2, eq = 2, Center = center)
Lt.Cfg(DSC, Range, a, M, 'DSC-1')

DSC = Lt.DSchwarz(Cube, Range, Index = 3, eq = 2, Center = center)
Lt.Cfg(DSC, Range, a, M, 'DSC-2')

DSC = Lt.DSchwarz(Cube, Range, Index = 1, eq = 2, Center = center)
Lt.Cfg(DSC, Range, a, M, 'DSC-3')

DSC = Lt.DSchwarz(Cube, Range, Index = 1, eq = 1, Center = center)
Lt.Cfg(DSC, Range, a, M, 'DSC-4')

DSC = Lt.DSchwarz(Cube, Range, Index = 3, eq = 1, Center = center)
Lt.Cfg(DSC, Range, a, M, 'DSC-5')

PSC = Lt.PSchwarz(Cube, Range, Index = 2, Center = center)
Lt.Cfg(PSC, Range, a, M, 'PSC')

ISC = Lt.IWPSchwarz(Cube, Range, Index = 2, Center = center)
Lt.Cfg(ISC, Range, a, M, 'ISC')

Octahedron = Lt.Octahedron(Cube, Range, Center = [1/2, 1/2, 1/2])
Lt.Cfg(Octahedron, Range, a, M, 'Octahedron')