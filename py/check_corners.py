from scipy import spatial as ss
import numpy as np
import datum_correction as dc
import csv


datum_points_file = '../data/usa_111k_datum_correction.csv'
tree_array, datum_array = dc.read_csv(datum_points_file)

tree = ss.KDTree(tree_array, 1)

pts = np.array([
    [-111.00, 43.00],
    [-111.00, 44.00],
    [-111.00, 45.00],
    [-108.00, 45.00],
    [-108.00, 44.00],
    [-108.00, 43.00],
    ]).astype(np.float)

results = []
for pt in pts:
    out_dict = dc.idw(pt, tree, datum_array)
    results.append(out_dict)

for result in results:
    print result["estimate"]
