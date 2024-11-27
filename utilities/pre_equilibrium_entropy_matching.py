import sys
import numpy as np
from scipy.interpolate import interp1d
import math

HBARC = 0.197326979


def read_file(filename):
    data = []
    comment_line = ""

    try:
        with open(filename, "r") as file:
            comment_line = next(file).strip()
            for line in file:
                # Split each line into a list of floating-point numbers
                row = [float(value) for value in line.split()]
                if len(row) > 0:
                    data.append(row)
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e} in file '{filename}'. Make sure all values are numbers.")
        sys.exit(1)

    return np.array(data), comment_line


def s_matching(e_kompost, nu_eff, file_path_eos):
    # Kompost entropy
    s_kompost = (
        (4.0 / 3.0)
        * (nu_eff * math.pi * math.pi / 30.0) ** (1.0 / 4.0)
        * (e_kompost ** (3.0 / 4.0))
    )  # GeV^3/4 / fm^9/4
    s_kompost = s_kompost * (HBARC ** (1 / 4.0))  # GeV/fm^2

    # EOS
    raw = np.fromfile(file_path_eos, dtype=(float, 4))
    # 4 columns
    # (energy density, local pressure, entropy density, local temperature)
    e = raw[:, 0]  # (GeV/fm^3)
    p = raw[:, 1]  # (GeV/fm^3)
    T = raw[:, 3]  # (GeV)

    entropy = [HBARC * (p + e) / T for p, e, T in zip(p, e, T)]  # GeV/fm^2
    # energy as a function of entropy [e(s)]
    energy_interp = interp1d(entropy, e, kind="linear", fill_value="extrapolate")
    # pressure as a function of energy [p(e)]
    p_interp = interp1d(e, p, kind="linear", fill_value="extrapolate")

    # Calculate the new energy and new bulk pressure
    e_Music = energy_interp(s_kompost)  # GeV/fm^3

    Bulk_pressure = e_kompost / 3.0 - p_interp(e_Music)  # GeV/fm^3

    return e_Music, Bulk_pressure


def main():
    # Check if all parameters are provided as command line arguments
    if len(sys.argv) != 6:
        print(
            "Usage: python script.py type_of_matching <eos_file> <filename_in> <filename_out> <nu_eff>"
        )
        sys.exit(1)

    # parse input parameters
    type_of_matching = int(sys.argv[1])
    file_path_eos = sys.argv[2]
    filename_in = sys.argv[3]
    filename_out = sys.argv[4]
    nu_eff = float(sys.argv[5])

    # e matching
    if type_of_matching == 0:
        # read the input file
        data, comment = read_file(filename_in)

    # s matching
    elif type_of_matching == 1:
        # read the input file
        data, comment = read_file(filename_in)

        # get the energy and calculate the new energy
        energy_density = data[:, 3]  # GeV/fm^3
        e_new, bulk_pressure = s_matching(
            energy_density, nu_eff, file_path_eos
        )  # GeV/fm^3

        # add the bulk pressure as a column in data
        data = np.column_stack((data, bulk_pressure))

        # Replace the column 4: change the energy kompost to the new energy
        data[:, 3] = e_new

    else:
        print("Problem with the variable type_of_matching, should be equal to 0 or 1.")

    # Write the final file: eta, x, y, new_energy, u^tau, u^x, u^y, u^eta,
    # pi^tautau, pi^taux, pi^tauy, pi^taueta, pi^xx, pi^xy, pi^xeta, pi^yy,
    # pi^yeta, pi^etaeta and Pi (in s matching).
    with open(filename_out, "w") as file:
        file.write(comment + "\n")
        for row in data:
            file.write(" ".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main()
