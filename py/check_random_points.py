from scipy import spatial as ss
import numpy as np
import datum_correction as dc
import csv


datum_points_file = '../data/points/usa_111k_datum_correction.csv'
tree_array, datum_array = dc.read_csv(datum_points_file)

tree = ss.KDTree(tree_array, 1)

pts = np.array([
  [-109.7811575956, 44.7235527953],
  [-110.0788108002, 44.1846681017],
  [-109.698944505, 44.6700888063],
  [-110.7611811343, 44.1718016999],
  [-110.5644522283, 44.4484130681],
  [-110.3774350439, 44.3090776079],
  [-110.345736631, 44.5508419213],
  [-110.2208009239, 44.6703142431],
  ]).astype(np.float32)


results = []
for pt in pts:
    out_dict = dc.idw(pt, tree, datum_array)
    results.append(out_dict)

for result in results:
    print result["estimate"]
