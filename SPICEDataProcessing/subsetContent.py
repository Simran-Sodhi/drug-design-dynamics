import h5py

# Open the HDF5 file in read mode
with h5py.File('SPICE-2.0.1.hdf5', 'r') as file:
    # Access the specific dataset
    subset_data1 = file['404339841/smiles']
    subset_data2 = file['404339841/subset']
    
    # Print the content of the 'SER/subset' dataset
    print("Content of 404339841/smiles:")
    print(subset_data1[()])  # Reads the entire dataset
    print("Content of 404339841/subset:")
    print(subset_data2[()])  # Reads the entire dataset