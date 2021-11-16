import numpy as np

from estimate_line import prepare_data
from estimate_line import compute_root_mean_square_deviation


def prepare_estimation_data(raw_measurements):
    z = np.matrix(raw_measurements[:, 2])
    h = raw_measurements
    for row in h:
        row[0, 2] = 1
    return z, h


# z = Hq, where:
#   z = (z_1, ..., z_n)^T
#   H = ([x_1, y_1, 1], ..., [x_n, y_n, 1])^T
#   q = (A, B, D) {z = Ax + By + D}
def compute_estimate_using_least_squares():
    raw_measurements = prepare_data("ROVER1.GK", "BASE1.GK", [0])
    z, h = prepare_estimation_data(raw_measurements)
    estimate = (h.getT() * h).I * h.getT() * z
    return estimate, compute_root_mean_square_deviation(h * estimate, z)


def compute_estimate_using_extended_least_squares():
    raw_measurements = prepare_data("ROVER1.GK", "BASE1.GK", [0])
    z, h = prepare_estimation_data(raw_measurements)

    u, sigma, v = np.linalg.svd(h)
    min_singular_value = sigma.min()
    identity = np.matrix(np.identity(h.shape[1]))
    estimate = -1 * (h.getT() * h - (min_singular_value ** 2) * identity).I * h.getT() * z
    return estimate, compute_root_mean_square_deviation(h * estimate, z)


if __name__ == "__main__":
    estimate_lsm, rmsd_lsm = compute_estimate_using_least_squares()
    print("Best estimate found using LSM is ", estimate_lsm.flatten())
    print("Root-mean-square deviation = ", rmsd_lsm)
    print("\n")

    estimate_elsm, rmsd_elsm = compute_estimate_using_extended_least_squares()
    print("Best estimate found using ELSM is ", estimate_lsm.flatten())
    print("Root-mean-square deviation = ", rmsd_lsm)
