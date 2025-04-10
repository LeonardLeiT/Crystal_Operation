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
num = 27
X, Y, Z = num, num, num

Range = [X, Y, Z]
Center = [Range[0] / 2, Range[1] / 2, Range[2] / 2]

PSC = {}
# 用中心生成坐标点，设置晶相来建立晶体
# 建立中心Kelvin晶体
axis = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
Cube = Lt.Grow_Crystal(Range, FCC, axis, Center = [1/2, 1/2, 1/2])
Cube = Lt.Box_Slice(Cube, Range, Sign = 1)
PSC[0] = Lt.PSchwarz(Cube, Range, Index = 0)
Lt.Cfg(PSC[0], Range, a, M, 'T-O')

# 建立010Kelvin晶体
axis = Lt.rotation_matrix(60, 90, 11)
print(axis)
Cube = Lt.Grow_Crystal(Range, FCC, axis, Center = [1/2, 1/2, 1/2])
Cube = Lt.Box_Slice(Cube, Range, Sign = 1)
PSC[1] = Lt.PSchwarz(Cube, Range, Index = 1)
Lt.Cfg(PSC[1], Range, a, M, 'T-E')


All = Lt.Merge_Group(PSC)
# origin = Lt.coordinate_min(All)
# print(origin)
# All = Lt.Move_to_origin(All, origin)
N = len(All)
print("Ave:", N / 2)
N = round(N / 1e4, 1)
print("All:", N)
Size = round((X + 0.5) * a / 10, 1)
print("Schwarz-p-Size:", Size)
name = 'Schwarz-P-' + str(N) + 'w-' + str(Size) + 'nm'
Lt.Cfg(All, Range, a, M, name)



