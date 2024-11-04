import h5py
import sys
import json

def categorize_datasets(file_name="SPICE-2.0.1.hdf5", attribute_name="subset"):
    subsets_dict = {}

    with h5py.File(file_name, 'r') as file:
        for group_name in file:
            group = file[group_name]
            # Check if the group has the "subset" attribute
            if attribute_name in group:
                subset_value = group[attribute_name][()][0]
                if isinstance(subset_value, bytes):
                    subset_value = subset_value.decode()  # Decode bytes to string if necessary

                # Add group name to the list of groups under this subset value
                if subset_value not in subsets_dict:
                    subsets_dict[subset_value] = []
                subsets_dict[subset_value].append(group_name)

                # Print the path and content for monitoring
                # print(f"Accessed {group_name}/{attribute_name}: {subset_value}")

    return subsets_dict

if __name__ == "__main__":
    # Expect the file name and attribute name as arguments
    file_name = sys.argv[1]
    attribute_name = sys.argv[2]

    # Call the function and print results as JSON
    subsets = categorize_datasets(file_name, attribute_name)