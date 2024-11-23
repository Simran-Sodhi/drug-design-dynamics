import h5py
import sys

if len(sys.argv) < 4:
    print("Please provide three strings as input.")
else:
    file = sys.argv[1]
    datasetName = sys.argv[2]
    AttributeName = sys.argv[3]

    
def get_string(file, datasetName, AttributeName):
    # Open the HDF5 file in read mode (file ends with .hdf5)
    with h5py.File(file, 'r') as file:
        # Access the dataset containing the SMILES string
        path = datasetName + '/' + AttributeName
        print(path)
        smiles_dataset = file[path]
        
        # Read the SMILES string
        smiles_array = smiles_dataset[()]

        smiles_string = smiles_array[0].decode('utf-8')
        
        return smiles_string

if __name__ == "__main__":
    smiles_string = get_string(file, datasetName, AttributeName)
    print(smiles_string)  # Output the string for Go to capture
