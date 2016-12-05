from scipy import spatial as ss
import numpy as np
import datum_correction as dc
import csv


def write_gdem(data, out_file):
    with open(out_file, 'wb') as f:
        writer = csv.writer(f)
        for row in data:
            writer.writerow(row)


datum_points_file = '../data/usa_111k_datum_correction.csv'
gdem_points_file = '../data/gdem_points.csv'

# Bring in datum info
tree_array, datum_array = dc.read_csv(datum_points_file)
# Bring in GDEM points
fields_dict = {
    "lon_index": 0,
    "lat_index": 1,
    "diff_index": 2
    }
gdem_tree_array, gdem_datum_array = dc.read_csv(gdem_points_file, fields_dict)

# Create Tree
datum_tree = ss.KDTree(tree_array, 1)

# Get estimate for each point, add it to the array
for pt in gdem_datum_array:
    coords = np.array([pt[0], pt[1]])
    out_dict = dc.idw(pt, tree, datum_array)
    pt[3] = out_dict["estimate"]

# Write to CSV
write_gdem(gdem_datum_array, '../data/gdem_correction.csv')
