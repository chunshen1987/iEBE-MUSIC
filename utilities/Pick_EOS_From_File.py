import sys
import pickle

# main function should accept a command line argument for the index of the EoS to be picked
def main():

    # get a command line argument for the index of the EoS to be picked
    if len(sys.argv) != 3:
        print("Usage: Pick_EOS_From_File.py filename index")
        sys.exit(1)

    filename = str(sys.argv[1])
    index = f"{int(sys.argv[2]):04d}"
    print(f"Pick_EOS_From_File: filename = {filename}, index = {index}")

    # load the EoS file
    with open(filename, 'rb') as f:
        EoS_dict = pickle.load(f)

    # pick the EoS with the given index
    if index in EoS_dict:
        data = EoS_dict[index]
    else:
        print(f"Index {index} not found in file {filename}")
        sys.exit(1)

    # write the EoS to a new binary file with columns e, P, T
    with open(f"EoS_{index}.bin", 'wb') as f:
        f.write(data.tobytes())

if __name__ == "__main__":
    main()
