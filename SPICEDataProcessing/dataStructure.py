import h5py

filename_hdf = 'SPICE-2.0.1.hdf5'
output_file = 'dataStructure.txt'

def h5_tree(val, pre='', file=None):
    items = len(val)
    for key, val in val.items():
        items -= 1
        if items == 0:
            # the last item
            if type(val) == h5py._hl.group.Group:
                file.write(pre + '└── ' + key + '\n')
                h5_tree(val, pre + '    ', file)
            else:
                try:
                    file.write(pre + '└── ' + key + ' (%d)\n' % len(val))
                except TypeError:
                    file.write(pre + '└── ' + key + ' (scalar)\n')
        else:
            if type(val) == h5py._hl.group.Group:
                file.write(pre + '├── ' + key + '\n')
                h5_tree(val, pre + '│   ', file)
            else:
                try:
                    file.write(pre + '├── ' + key + ' (%d)\n' % len(val))
                except TypeError:
                    file.write(pre + '├── ' + key + ' (scalar)\n')

with h5py.File(filename_hdf, 'r') as hf, open(output_file, 'w') as file:
    file.write(str(hf) + '\n')
    h5_tree(hf, file=file)
