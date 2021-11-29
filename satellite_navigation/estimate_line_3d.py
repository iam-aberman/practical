from estimate_line import prepare_data

import numpy as np


# k = H * q
#   k = (x_1 - x_0, ..., x_n - x_0)^T, where x -> y -> z
#   H = (1, ..., n)^T
#   q = x -> y -> z
def prepare_estimation_data(raw_measurements, axis_name):
    axis_no = -1
    if axis_name == "x":
        axis_no = 1
    elif axis_name == "y":
        axis_no = 2
    elif axis_name == "z":
        axis_name = 3
    else:
        raise IndexError
    k = raw_measurements[1:, axis_no]
    for row in k:
        row[0, 0] -= raw_measurements[0, axis_no]

    h = raw_measurements[1:, 0]
    return k, h


def solve_lsm_estimation_task(k, h):
    return (h.getT() * h).I * h.getT() * k


def solve_glsm_estimation_task(k, h):
    extended_h = np.hstack((h, k))
    u, sigma, v = np.linalg.svd(extended_h)
    min_singular_value = sigma.min()
    identity = np.matrix(np.identity(h.shape[1]))
    return (h.getT() * h - (min_singular_value ** 2) * identity).I * h.getT() * k


def print_results(point, direction, annotation):
    print(annotation)
    print("Init point = (", point[0, 1], point[0, 2], point[0, 3], ")")
    print("Direction  = (", direction[0], direction[1], direction[2], ")")


def estimate_direction(raw_measurements):
    lsm_direction_estimate = []
    glsm_direction_estimate = []
    for axis_name in ["x", "y", "z"]:
        k, h = prepare_estimation_data(raw_measurements, axis_name)
        lsm_direction_estimate.append(solve_lsm_estimation_task(k, h)[0, 0])
        glsm_direction_estimate.append(solve_glsm_estimation_task(k, h)[0, 0])

    print_results(raw_measurements[0], lsm_direction_estimate, "LSM")
    print_results(raw_measurements[0], glsm_direction_estimate, "GLSM")


if __name__ == "__main__":
    raw_measurements = prepare_data("ROVER1.GK", "BASE1.GK", [ ])
    estimate_direction(raw_measurements)
