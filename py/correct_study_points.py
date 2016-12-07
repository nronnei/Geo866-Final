from scipy import spatial as ss
import numpy as np
import datum_correction as dc
import csv
import os


def read_points(in_file):
  with open(in_file, 'rb') as f:
    reader = csv.reader(f, delimiter=",")
    data = []
    # Skip header
    reader.next()
    for row in reader:
      data.append(row)
    return data


def write_points(data, out_file):
  with open(out_file, 'wb') as f:
    writer = csv.writer(f)
    for row in data:
      writer.writerow(row)


os.chdir("/home/nronnei/gis/class/spatial_analysis/final_project/")
datum_points_file = './data/points/usa_111k_datum_correction.csv'
tree_array, datum_array = dc.read_csv(datum_points_file)
study_points_file = './data/points/study_points_raw.csv'
study_points = read_points(study_points_file)


tree = ss.KDTree(tree_array, 1)


for pt in study_points:
  loc = np.array([pt[0], pt[1]]).astype(np.float32)
  out_dict = dc.idw(loc, tree, datum_array)
  pt[2] = np.float32(out_dict["estimate"])
  corrected_elevation = np.float32(pt[2] + np.float32(pt[3]))
  pt.append(corrected_elevation)

head = ['x', 'y', 'dc', 'elev', 'c_elev']
study_points.insert(0, head)
corrected_points_file = './data/points/study_points_Python.csv'
write_points(study_points, corrected_points_file)
