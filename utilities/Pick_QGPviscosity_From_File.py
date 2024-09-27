import sys
import pickle


def fetchShearViscosity1D(database: str, id: int) -> str:
    """
        This function fetches a eta/s(T) with `id` from the database file.
        It returns the shearFileName
    """
    index = f"{id:04d}"

    # load the shear viscosity file
    with open(database, 'rb') as f:
        shear_dict = pickle.load(f)

    # pick the EoS with the given index
    if index in shear_dict:
        data = shear_dict[index]
    else:
        print(f"Index {index} not found in file {database}")
        sys.exit(1)

    # write the EoS to a new binary file with columns e, P, T
    shearFileName = f"shear_{index}.bin"
    with open(shearFileName, 'wb') as f:
        f.write(data.tobytes())

    return shearFileName


def fetchBulkViscosity1D(database: str, id: int) -> str:
    """
        This function fetches a zeta/s(T) with `id` from the database file.
        It returns the bulkFileName
    """
    index = f"{id:04d}"

    # load the bulk viscosity file
    with open(database, 'rb') as f:
        bulk_dict = pickle.load(f)

    # pick the EoS with the given index
    if index in bulk_dict:
        data = bulk_dict[index]
    else:
        print(f"Index {index} not found in file {database}")
        sys.exit(1)

    # write the EoS to a new binary file with columns e, P, T
    bulkFileName = f"bulk_{index}.bin"
    with open(bulkFileName, 'wb') as f:
        f.write(data.tobytes())

    return bulkFileName


if __name__ == "__main__":
    try:
        filename1 = str(sys.argv[1])
        filename2 = str(sys.argv[2])
        index = int(sys.argv[3])
    except IndexError:
        print(f"Usage: {sys.argv[0]} filename1 filename2 index")
        sys.exit(1)

    fetchShearViscosity1D(filename1, index)
    fetchBulkViscosity1D(filename2, index)
