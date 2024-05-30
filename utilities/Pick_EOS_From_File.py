import sys
import pickle

def fetch_an_EOS(database: str, id: int) -> str:
    """
        This function fetch an EoS with `id` from the database file.
        It return the eosFileName
    """
    index = f"{id:04d}"
    #print(f"Pick_EOS_From_File: filename = {database}, index = {index}")

    # load the EoS file
    with open(database, 'rb') as f:
        EoS_dict = pickle.load(f)

    # pick the EoS with the given index
    if index in EoS_dict:
        data = EoS_dict[index]
    else:
        print(f"Index {index} not found in file {database}")
        sys.exit(1)

    # write the EoS to a new binary file with columns e, P, T
    eosFileName = f"EoS_{index}.bin"
    with open(eosFileName, 'wb') as f:
        f.write(data.tobytes())

    return eosFileName

if __name__ == "__main__":
    try:
        filename = str(sys.argv[1])
        index = int(sys.argv[2])
    except IndexError:
        print("Usage: Pick_EOS_From_File.py filename index")
        sys.exit(1)

    fetch_an_EOS(filename, index)
