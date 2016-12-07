from scipy import spatial as ss
import numpy as np
import csv
import sys
import time

def find_coord_names(header):
    """
    This function attempts to find the the indices of the pertinent fields based
    on whether or not common/logical names appear in the header row of a CSV.

    :param header: A list containing the values of the header row of a CSV.
    :return: A dictionary containing the field indices for lon, lat, and diff.
    """
    formatted_header = []
    field_dict = {
        "lon_index": None,
        "lat_index": None,
        "diff_index": None,
    }
    for h in header:
        formatted_header.append(h.lower())

    try:
        # Get lon index
        if "lon" in formatted_header:
            field_dict["lon_index"] = formatted_header.index("lon")
        elif "x" in formatted_header:
            field_dict["lon_index"] = formatted_header.index("x")
        elif "xloc" in formatted_header:
            field_dict["lon_index"] = formatted_header.index("xloc")
        elif "locx" in formatted_header:
            field_dict["lon_index"] = formatted_header.index("locx")

        # Get lat index
        if "lat" in formatted_header:
            field_dict["lat_index"] = formatted_header.index("lat")
        elif "y" in formatted_header:
            field_dict["lat_index"] = formatted_header.index("y")
        elif "yloc" in formatted_header:
            field_dict["lat_index"] = formatted_header.index("yloc")
        elif "locy" in formatted_header:
            field_dict["lat_index"] = formatted_header.index("locy")

        # Get diff index
        if "diff" in formatted_header:
            field_dict["diff_index"] = formatted_header.index("diff")
        elif "datumcorrection" in formatted_header:
            field_dict["diff_index"] = formatted_header.index("datumcorrection")
        elif "dat_cor" in formatted_header:
            field_dict["diff_index"] = formatted_header.index("dat_cor")
        elif "dat_corr" in formatted_header:
            field_dict["diff_index"] = formatted_header.index("dat_corr")
        elif "correction" in formatted_header:
            field_dict["diff_index"] = formatted_header.index("correction")
        elif "datum" in formatted_header:
            field_dict["diff_index"] = formatted_header.index("datum")

    except KeyError:
        print ("Key Error: failed to build fields dictionary.")
        sys.exit(1)

    for key in field_dict.iterkeys():
        if field_dict[key] is None:
            print ("Failed to build fields dictionary.")
            print ("No value found for key '{0}'.".format(key))
            print ("Please create a fields dictionary with the appropriate field indices "
                   "to pass into the read_csv function or rename your fields.")
            print ('Example of fields dictionary setup:\n'
                   '"fields" = {\n'
                   '    "lon_index": 0,\n'
                   '    "lat_index": 1,\n'
                   '    "diff_index": 2\n'
                   '}\n'
                   'NOTE: the keys MUST be the same as those listed here.)')
            sys.exit(1)

    return field_dict


def read_csv(in_file_name, fields=None):
    """
    Read in a CSV file and prepare the data for use with SciPy's KD Tree.

    :param in_file_name: Path to input CSV (string).
    :param fields: A dictionary containing the indices of the relevant fields.

     Example of fields dictionary setup:
     fields = {
         "lon_index": 0,
         "lat_index": 1,
         "diff_index": 2
         }
    NOTE: the keys MUST be the same as those listed here.

    :return: Three numpy arrays zipped for use in SciPy's KD Tree.
    """

    # Open the file, create CSV Reader
    reader = csv.reader(open(in_file_name, 'rb'), delimiter=',')

    # Set up empty lists to hold values for each field
    lon = []
    lat = []
    diff = []

    # If no fields are supplied, attempt to find them in the header
    head = next(reader)
    if fields is None:
        fields = find_coord_names(head)
    else:
        pass

    # Read in spreadsheet and place values in their respective lists
    for row in reader:
        lon.append(row[fields['lon_index']])
        lat.append(row[fields['lat_index']])
        diff.append(row[fields['diff_index']])

    # Make sure arrays are in proper type
    np_lon = np.array(lon).astype(np.float)
    np_lat = np.array(lat).astype(np.float)
    np_diff = np.array(diff).astype(np.float)
    tree_array = zip(np_lon, np_lat)
    datum_array = zip(np_lon, np_lat, np_diff)

    # Save some memory
    del lon, lat, diff, np_lon, np_lat, np_diff

    return tree_array, datum_array


def write_csv(out_file, data_dict, fields=None):
    """
    Writes a list of dictionaries to an output file. fields can be passed, or it will
    simply be extracted from the dictionary keys. The advantage of passing it is that
    the field order can be set, or a subset of fields can be saved. The disadvantage is
    that if fieldnames are passed that aren't keys, this won't work!

    :param out_file:
    :param data_dict:
    :param fields:
    :return:
    """
    import csv

    if fields is None:
        fields = [key for key in data_dict[0].iterkeys()]
    of = open(out_file, 'w')
    of.write(','.join(fields)+'\n')  # This writes the header row to the file.
    try:
        writer = csv.DictWriter(of, fieldnames=fields)
    except ValueError:  # if the field names don't match what's in the dict
        print "Field names don't match dictionary keys!"
        return
    for rec in data_dict:
        writer.writerow(rec)    # This writes all the rows in dataDict to file

    of.close()
    return


def calc_spherical_dist(lon1, lat1, lon2, lat2):
    """
    Calculates approximate ground distance, in kilometers, between the two locations.
    This calculates the spherical great circle distance between the two coordinates.
    For the algorithm see: http://en.wikipedia.org/wiki/Great-circle_distance

    :param lon1:
    :param lat1:
    :param lon2:
    :param lat2:
    :return:
    """
    import math
    earth_radius = 6371.01  # average radius, in km
    # First calculate coordinates in radians
    lon1 *= (math.pi / 180.0)
    lat1 *= (math.pi / 180.0)
    lon2 *= (math.pi / 180.0)
    lat2 *= (math.pi / 180.0)

    # Now calculate the angular distance
    cos_delta_lon = math.cos(lon2-lon1)
    sin_delta_lon = math.sin(lon2-lon1)
    right_numerator = ((math.cos(lat1) * math.sin(lat2)) - (math.sin(lat1) * math.cos(lat2) * cos_delta_lon))**2
    numerator = ((math.cos(lat2) * sin_delta_lon)**2 + right_numerator)**0.5
    denominator = (math.sin(lat1) * math.sin(lat2)) + (math.cos(lat1) * math.cos(lat2) * cos_delta_lon)
    ad = math.atan2(numerator, denominator)
    distance = earth_radius * ad
    return distance


def idw(input_point_array, tree, datum_array, k=2):
    """
    Performs IDW interpolation to a coordinate (Lon, Lat dictionary) using
    coordinates stored in a 2D point tree structure. Distances are
    spherical, as this function assumes coordinates are lat-lon.

    :param input_point_array: A NumPy array structured as follows: np.array([lat, lon])
    :param tree: A KD Tree created from scipy.spatial.KDTree()
    :param datum_array: The NumPy Array generated by read_csv() which contains ALL the csv data.
    :param k: The IDW power function. This could range continuously from 0-3.
    :return:
    """

    idw_dict = {'dist': [], 'dcv': [], 'w': [], 'wdcv': []}

    # Get the 4 nearest neighbors of input_point_query
    query_result = tree.query(input_point_array, 6)
    # For each of the nearest neighbors...
    for result in query_result[1]:
        datum_point_array = datum_array[result]
        # Calculate distance between the input point (lon1/lat1) and the nearest neighbor (lon2/lat2)
        idw_dict['dist'].append(calc_spherical_dist(input_point_array[0], input_point_array[1],
                                                    datum_point_array[0], datum_point_array[1]))
        # Calculate the Inverse Distance Weight based on the distance between the points and the power function
        idw_dict['w'].append(1 / (idw_dict['dist'][-1] + 0.0)**k)
        # Get the Datum Correction Value from the nearest neighbor
        idw_dict['dcv'].append(datum_point_array[2])
        # Calculate the Weighted Datum Correction Value
        idw_dict['wdcv'].append(idw_dict['w'][-1] * idw_dict['dcv'][-1])

    if sum(idw_dict['w']) > 0:
        idw_dict['estimate'] = sum(idw_dict['wdcv'])/sum(idw_dict['w'])
    else:
        idw_dict['estimate'] = None
    return idw_dict


# #########################
# # Main Code Starts Here #
# #########################
#
# #####
# # Setup input variables
# #####
# # Check for input arguments and set variables. No error handling here!
# if len(sys.argv) == 4:
#     datum_points_file = sys.argv[1]
#     sample_points_file = sys.argv[2]
#     out_file = sys.argv[3]
# else:
#     datum_points_file = '../original/usa_111k_datum_correction.csv'  # The only one that matters right now (5/25/16)
#     sample_points_file = 'us_37ksample_pts.csv'  # Doesn't matter right now (5/25/16)
#     out_file = 'datum_conversion/usa_37k_datum_correction.csv'  # Doesn't matter right now (5/25/16)
# # Load datum correction and coordinates into NumPy Arrays
# print('\nLoading datum file....')
# t0 = time.time()
# tree_array, datum_array = read_csv(datum_points_file)
# time_elapsed = time.time() - t0
# print ("Done.")
# print ("Time elapsed: " + str(time_elapsed))
# print ("------")
# print ("Example Datum Array Row:")
# print (datum_array[5])
# print ("Example Tree Array Row:")
# print (tree_array[5])
#
# #####
# # Generate the KD Tree using SciPy Spatial
# #####
# print ('\nGenerating tree...')
# t0 = time.time()
# tree = ss.KDTree(tree_array, 1)
# time_elapsed = time.time() - t0
# print ("Done.")
# print ("Time elapsed: " + str(time_elapsed))
#
# #####
# # Calculate IDW
# #####
# #  Set up test point as a NumPy array  for input to IDW/tree.query
# pt = np.array([-93.6, 30.05]).astype(np.float)
# #  Perform/Time IDW calculation
# print ('\nCalculating IDW...')
# t0 = time.time()
# out_dict = idw(pt, tree, datum_array)
# time_elapsed = time.time() - t0
# print ("Done.")
# print ("Time elapsed: " + str(time_elapsed))
# print ("------")
# print ("Estimated Datum Correction:")
# print (str(out_dict['estimate']) + " for point (" + str(pt[0]) + ", " + str(pt[1]) + ")")

'''
Testing approach on a single point:
Given a Lon, Lat coordinate pair, get the 4 nearest neighbors

Pennsylvania: [-77.42958, 39.59458]
Oregon: [-120, 45]
Louisiana: [-93.6, 30.05]
'''
