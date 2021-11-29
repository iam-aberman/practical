import numpy as np

from math import sqrt


def prepare_data(rover_filename, base_filename, columns_to_delete):
    full_rover_data = np.matrix(np.loadtxt(rover_filename))
    full_base_data = np.matrix(np.loadtxt(base_filename))
    rows_to_delete_for_rover = list(range(0, 870)) + list(range(1031, full_rover_data.shape[0]))
    row_to_delete_for_base   = list(range(0, 870)) + list(range(1031, full_base_data.shape[0]))
    full_rover_data = np.delete(full_rover_data, rows_to_delete_for_rover, axis=0)
    full_rover_data = np.delete(full_rover_data, columns_to_delete, axis=1)
    full_base_data = np.delete(full_base_data, row_to_delete_for_base, axis=0)
    full_base_data = np.delete(full_base_data, columns_to_delete, axis=1)
    ret_val = full_rover_data - full_base_data
    for i in range(ret_val.shape[0]):
        ret_val[i, 0] = full_rover_data[i, 0]
    return ret_val

def prepare_estimation_data(raw_measurements):
    y = np.matrix(raw_measurements[:, 1])
    h = raw_measurements
    for row in h:
        row[0, 1] = 1
    return y, h


def are_collinear(lhs, rhs):
    return abs(lhs[0, 0] * rhs[0, 1] - lhs[0, 1] * rhs[0, 0]) < 1e-12


def calculate_l1_norm(vector):
    l1_norm = 0
    for i in range(vector.shape[0]):
        l1_norm += abs(vector[i, 0])
    return l1_norm


def compute_root_mean_square_deviation(estimate, measurements):
    rmsd = 0
    for i in range(estimate.shape[0]):
        rmsd += (estimate[i, 0] - measurements[i, 0]) ** 2
    return sqrt(rmsd / estimate.shape[0])


# y = Hq, where:
#   y = (y_1, ..., y_n)^T,
#   H = ([x_1, 1], ..., [x_n, 1])^T,
#   q = (k, b)^T {y = kx + b}
def compute_estimate():
    raw_measurements = prepare_data("ROVER1.GK", "BASE1.GK", [0, 3])
    y, h = prepare_estimation_data(raw_measurements)
    min_val = float('inf')
    best_estimate = np.matrix([ [float('inf')], [float('inf')] ])

    n_rows = h.shape[0]
    for i in range(n_rows):
        for j in range(n_rows):
            if are_collinear(h[i], h[j]):
                continue
            cur_h = np.matrix([
                [h[i, 0], h[i, 1]],
                [h[j, 0], h[j, 1]]
            ])
            cur_y = np.matrix([
                [y[i, 0]],
                [y[j, 0]]
            ])
            cur_est = cur_h.I * cur_y
            cur_res = y - h * cur_est
            if (calculate_l1_norm(cur_res) < min_val):
                min_val = calculate_l1_norm(cur_res)
                best_estimate = cur_est
    return best_estimate, compute_root_mean_square_deviation(h * best_estimate, y)


if __name__ == "__main__":
    estimate, rmsd = compute_estimate()
    print("Best estimate found using LAD is ", estimate.flatten())
    print("Root-mean-square deviation = ", rmsd)
