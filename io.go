package main

import (
	"encoding/json"
	"log"
	"os/exec"
	"strings"

	// "gonum.org/v1/hdf5"
	"github.com/sbinet/go-hdf5"
)

/*
atomic_numbers, Shape: (38,), Dtype: int16
conformations, Shape: (4, 38, 3), Dtype: float32
dft_total_energy, Shape: (4,), Dtype: float64
dft_total_gradient, Shape: (4, 38, 3), Dtype: float32
formation_energy, Shape: (4,), Dtype: float64
mayer_indices, Shape: (4, 38, 38), Dtype: float32
mbis_charges, Shape: (4, 38, 1), Dtype: float32
mbis_dipoles, Shape: (4, 38, 3), Dtype: float32
mbis_octupoles, Shape: (4, 38, 3, 3, 3), Dtype: float32
mbis_quadrupoles, Shape: (4, 38, 3, 3), Dtype: float32
scf_dipole, Shape: (4, 3), Dtype: float32
scf_quadrupole, Shape: (4, 3, 3), Dtype: float32
smiles, Shape: (1,), Dtype: object
subset, Shape: (1,), Dtype: object
wiberg_lowdin_indices, Shape: (4, 38, 38), Dtype: float32
*/

// CategorizeDatasets() takes fileName and an attributeName as input, and returns a map whose key is name of subsets (stored in attributeName), and value is the collection of datasets names who fall into the subset.
func CategorizeDatasets(fileName, attributeName string) map[string][]string {
	cmd := exec.Command("python3", "categorize_datasets.py", fileName, attributeName)

	output, err := cmd.Output()
	if err != nil {
		log.Fatalf("Error executing Python script: %v", err)
	}

	// Parse the JSON output from Python into a Go map
	var subsetMap map[string][]string
	if err := json.Unmarshal(output, &subsetMap); err != nil {
		log.Fatalf("Error parsing JSON output: %v", err)
	}

	return subsetMap
}

// GetSmiles() takes into a dataset name, and returns smiles in the dataset.
// smiles, Shape: (1,), Dtype: string
func GetSmiles(file, datasetName, attributeName string) string {
	// Command to run the Python script
	cmd := exec.Command("python3", "get_string.py", file, datasetName, attributeName)

	// Run the command and capture the output
	out, err := cmd.Output()
	if err != nil {
		log.Fatalf("Error running Python script: %v", err)
	}

	// Convert the output to a string and trim any trailing whitespace
	smilesString := strings.TrimSpace(string(out))
	// fmt.Printf("SMILES string from Python: %s\n", smilesString)

	return smilesString
}

// GetScf_quadrupole() takes into a dataset name, and returns scf_quadrupole in the dataset.
// scf_quadrupole, Shape: (4, 3, 3), Dtype: float32
func GetScf_quadrupole(fileName, datasetName, attributeName string) [][][]float32 {
	file, err := hdf5.OpenFile(fileName, hdf5.F_ACC_RDONLY)
	if err != nil {
		log.Fatal("Error opening file:", err)
	}
	defer file.Close()

	// access the dataset
	path := datasetName + "/" + attributeName
	dataset, err := file.OpenDataset(path)
	if err != nil {
		log.Fatal("Error opening dataset:", err)
	}
	defer dataset.Close()

	// Get the dimensions of the 3-D matrix
	space := dataset.Space()
	dims, _, err := space.SimpleExtentDims()
	if err != nil {
		log.Fatalf("failed to get dataset dimensions: %v", err)
	}
	totalSize := int(dims[0] * dims[1] * dims[2])

	// Create a 1-D slice to hold the data
	flatMatrix := make([]float32, totalSize)

	// Read the data from the dataset into the flat slice
	if err := dataset.Read(&flatMatrix); err != nil {
		log.Fatalf("failed to read dataset: %v", err)
	}

	// Reshape the flat slice into a 3-D slice
	matrix := make([][][]float32, dims[0])
	for i := range matrix {
		matrix[i] = make([][]float32, dims[1])
		for j := range matrix[i] {
			matrix[i][j] = flatMatrix[(i*int(dims[1])*int(dims[2]))+(j*int(dims[2])) : (i*int(dims[1])*int(dims[2]))+(j+1)*int(dims[2])]
		}
	}

	// Print the 3-D matrix (or process it as needed)
	// fmt.Println("last element of 3-D matrix data:", matrix)

	return matrix
}

// GetScf_dipole() takes into a dataset name, and returns scf_dipole in the dataset.
// scf_dipole, Shape: (4, 3), Dtype: float32
func GetScf_dipole(fileName, datasetName, attributeName string) [][]float32 {
	file, err := hdf5.OpenFile(fileName, hdf5.F_ACC_RDONLY)
	if err != nil {
		log.Fatal("Error opening file:", err)
	}
	defer file.Close()

	// Access the dataset
	path := datasetName + "/" + attributeName
	dataset, err := file.OpenDataset(path)
	if err != nil {
		log.Fatal("Error opening dataset:", err)
	}
	defer dataset.Close()

	// Get the dimensions of the 2D matrix
	space := dataset.Space()
	dims, _, err := space.SimpleExtentDims()
	if err != nil {
		log.Fatalf("Failed to get dataset dimensions: %v", err)
	}
	rows, cols := int(dims[0]), int(dims[1])
	totalSize := rows * cols

	// Create a 1-D slice to hold the data
	flatMatrix := make([]float32, totalSize)

	// Read the data from the dataset into the flat slice
	if err := dataset.Read(&flatMatrix); err != nil {
		log.Fatalf("Failed to read dataset: %v", err)
	}

	// Reshape the flat slice into a 2D slice
	matrix := make([][]float32, rows)
	for i := 0; i < rows; i++ {
		matrix[i] = flatMatrix[i*cols : (i+1)*cols]
	}

	return matrix
}

// GetMbis_quadrupoles() takes into a dataset name, and returns mbis_octupoles in the dataset.
// mbis_quadrupoles, Shape: (4, 38, 3, 3), Dtype: float32
func GetMbis_quadrupoles(fileName, datasetName, attributeName string) [][][][]float32 {
	file, err := hdf5.OpenFile(fileName, hdf5.F_ACC_RDONLY)
	if err != nil {
		log.Fatal("Error opening file:", err)
	}
	defer file.Close()

	// Access the dataset
	path := datasetName + "/" + attributeName
	dataset, err := file.OpenDataset(path)
	if err != nil {
		log.Fatal("Error opening dataset:", err)
	}
	defer dataset.Close()

	// Get the dimensions of the 4D matrix
	space := dataset.Space()
	dims, _, err := space.SimpleExtentDims()
	if err != nil {
		log.Fatalf("Failed to get dataset dimensions: %v", err)
	}
	dim0, dim1, dim2, dim3 := int(dims[0]), int(dims[1]), int(dims[2]), int(dims[3])
	totalSize := dim0 * dim1 * dim2 * dim3

	// Create a 1-D slice to hold the data
	flatMatrix := make([]float32, totalSize)

	// Read the data from the dataset into the flat slice
	if err := dataset.Read(&flatMatrix); err != nil {
		log.Fatalf("Failed to read dataset: %v", err)
	}

	// Reshape the flat slice into a 4D slice
	matrix := make([][][][]float32, dim0)
	for i := 0; i < dim0; i++ {
		matrix[i] = make([][][]float32, dim1)
		for j := 0; j < dim1; j++ {
			matrix[i][j] = make([][]float32, dim2)
			for k := 0; k < dim2; k++ {
				start := (i*dim1*dim2*dim3 + j*dim2*dim3 + k*dim3)
				end := start + dim3
				matrix[i][j][k] = flatMatrix[start:end]
			}
		}
	}

	return matrix
}

// GetMbis_octupoles() takes into a dataset name, and returns mbis_octupoles in the dataset.
// mbis_octupoles, Shape: (4, 38, 3, 3, 3), Dtype: float32
func GetMbis_octupoles(fileName, datasetName, attributeName string) [][][][][]float32 {
	file, err := hdf5.OpenFile(fileName, hdf5.F_ACC_RDONLY)
	if err != nil {
		log.Fatal("Error opening file:", err)
	}
	defer file.Close()

	// Access the dataset
	path := datasetName + "/" + attributeName
	dataset, err := file.OpenDataset(path)
	if err != nil {
		log.Fatal("Error opening dataset:", err)
	}
	defer dataset.Close()

	// Get the dimensions of the 5D matrix
	space := dataset.Space()
	dims, _, err := space.SimpleExtentDims()
	if err != nil {
		log.Fatalf("Failed to get dataset dimensions: %v", err)
	}

	// Ensure that we have exactly 5 dimensions
	if len(dims) != 5 {
		log.Fatalf("Expected a 5-dimensional dataset, but got %d dimensions", len(dims))
	}
	d1, d2, d3, d4, d5 := int(dims[0]), int(dims[1]), int(dims[2]), int(dims[3]), int(dims[4])
	totalSize := d1 * d2 * d3 * d4 * d5

	// Create a 1-D slice to hold the data
	flatMatrix := make([]float32, totalSize)

	// Read the data from the dataset into the flat slice
	if err := dataset.Read(&flatMatrix); err != nil {
		log.Fatalf("Failed to read dataset: %v", err)
	}

	// Reshape the flat slice into a 5D slice
	matrix := make([][][][][]float32, d1)
	for i := 0; i < d1; i++ {
		matrix[i] = make([][][][]float32, d2)
		for j := 0; j < d2; j++ {
			matrix[i][j] = make([][][]float32, d3)
			for k := 0; k < d3; k++ {
				matrix[i][j][k] = make([][]float32, d4)
				for l := 0; l < d4; l++ {
					start := ((i * d2 * d3 * d4 * d5) + (j * d3 * d4 * d5) + (k * d4 * d5) + (l * d5))
					end := start + d5
					matrix[i][j][k][l] = flatMatrix[start:end]
				}
			}
		}
	}

	return matrix
}

// GetWiberg_lowdin_indices() takes into a dataset name, and returns wiberg_lowdin_indices in the dataset.
// wiberg_lowdin_indices, Shape: (4, 38, 38), Dtype: float32
func GetWiberg_lowdin_indices(fileName, datasetName, attributeName string) [][][]float32 {
	file, err := hdf5.OpenFile(fileName, hdf5.F_ACC_RDONLY)
	if err != nil {
		log.Fatal("Error opening file:", err)
	}
	defer file.Close()

	// access the dataset
	path := datasetName + "/" + attributeName
	dataset, err := file.OpenDataset(path)
	if err != nil {
		log.Fatal("Error opening dataset:", err)
	}
	defer dataset.Close()

	// Get the dimensions of the 3-D matrix
	space := dataset.Space()
	dims, _, err := space.SimpleExtentDims()
	if err != nil {
		log.Fatalf("failed to get dataset dimensions: %v", err)
	}
	totalSize := int(dims[0] * dims[1] * dims[2])

	// Create a 1-D slice to hold the data
	flatMatrix := make([]float32, totalSize)

	// Read the data from the dataset into the flat slice
	if err := dataset.Read(&flatMatrix); err != nil {
		log.Fatalf("failed to read dataset: %v", err)
	}

	// Reshape the flat slice into a 3-D slice
	matrix := make([][][]float32, dims[0])
	for i := range matrix {
		matrix[i] = make([][]float32, dims[1])
		for j := range matrix[i] {
			matrix[i][j] = flatMatrix[(i*int(dims[1])*int(dims[2]))+(j*int(dims[2])) : (i*int(dims[1])*int(dims[2]))+(j+1)*int(dims[2])]
		}
	}

	// Print the 3-D matrix (or process it as needed)
	// fmt.Println("last element of 3-D matrix data:", matrix)

	return matrix
}

// GetMbis_dipoles() takes into a dataset name, and returns mbis_dipoles in the dataset.
// mbis_dipoles, Shape: (4, 38, 3), Dtype: float32
func GetMbis_dipoles(fileName, datasetName, attributeName string) [][][]float32 {
	file, err := hdf5.OpenFile(fileName, hdf5.F_ACC_RDONLY)
	if err != nil {
		log.Fatal("Error opening file:", err)
	}
	defer file.Close()

	// access the dataset
	path := datasetName + "/" + attributeName
	dataset, err := file.OpenDataset(path)
	if err != nil {
		log.Fatal("Error opening dataset:", err)
	}
	defer dataset.Close()

	// Get the dimensions of the 3-D matrix
	space := dataset.Space()
	dims, _, err := space.SimpleExtentDims()
	if err != nil {
		log.Fatalf("failed to get dataset dimensions: %v", err)
	}
	totalSize := int(dims[0] * dims[1] * dims[2])

	// Create a 1-D slice to hold the data
	flatMatrix := make([]float32, totalSize)

	// Read the data from the dataset into the flat slice
	if err := dataset.Read(&flatMatrix); err != nil {
		log.Fatalf("failed to read dataset: %v", err)
	}

	// Reshape the flat slice into a 3-D slice
	matrix := make([][][]float32, dims[0])
	for i := range matrix {
		matrix[i] = make([][]float32, dims[1])
		for j := range matrix[i] {
			matrix[i][j] = flatMatrix[(i*int(dims[1])*int(dims[2]))+(j*int(dims[2])) : (i*int(dims[1])*int(dims[2]))+(j+1)*int(dims[2])]
		}
	}

	// Print the 3-D matrix (or process it as needed)
	// fmt.Println("last element of 3-D matrix data:", matrix)

	return matrix
}

// GetMbis_charges() takes into a dataset name, and returns mbis_charges in the dataset.
// mbis_charges, Shape: (4, 38, 1), Dtype: float32
func GetMbis_charges(fileName, datasetName, attributeName string) [][][]float32 {
	file, err := hdf5.OpenFile(fileName, hdf5.F_ACC_RDONLY)
	if err != nil {
		log.Fatal("Error opening file:", err)
	}
	defer file.Close()

	// access the dataset
	path := datasetName + "/" + attributeName
	dataset, err := file.OpenDataset(path)
	if err != nil {
		log.Fatal("Error opening dataset:", err)
	}
	defer dataset.Close()

	// Get the dimensions of the 3-D matrix
	space := dataset.Space()
	dims, _, err := space.SimpleExtentDims()
	if err != nil {
		log.Fatalf("failed to get dataset dimensions: %v", err)
	}
	totalSize := int(dims[0] * dims[1] * dims[2])

	// Create a 1-D slice to hold the data
	flatMatrix := make([]float32, totalSize)

	// Read the data from the dataset into the flat slice
	if err := dataset.Read(&flatMatrix); err != nil {
		log.Fatalf("failed to read dataset: %v", err)
	}

	// Reshape the flat slice into a 3-D slice
	matrix := make([][][]float32, dims[0])
	for i := range matrix {
		matrix[i] = make([][]float32, dims[1])
		for j := range matrix[i] {
			matrix[i][j] = flatMatrix[(i*int(dims[1])*int(dims[2]))+(j*int(dims[2])) : (i*int(dims[1])*int(dims[2]))+(j+1)*int(dims[2])]
		}
	}

	// Print the 3-D matrix (or process it as needed)
	// fmt.Println("last element of 3-D matrix data:", matrix)

	return matrix
}

// GetMayer_indices() takes into a dataset name, and returns mayer_indices in the dataset.
// mayer_indices, Shape: (4, 38, 38), Dtype: float32
func GetMayer_indices(fileName, datasetName, attributeName string) [][][]float32 {
	file, err := hdf5.OpenFile(fileName, hdf5.F_ACC_RDONLY)
	if err != nil {
		log.Fatal("Error opening file:", err)
	}
	defer file.Close()

	// access the dataset
	path := datasetName + "/" + attributeName
	dataset, err := file.OpenDataset(path)
	if err != nil {
		log.Fatal("Error opening dataset:", err)
	}
	defer dataset.Close()

	// Get the dimensions of the 3-D matrix
	space := dataset.Space()
	dims, _, err := space.SimpleExtentDims()
	if err != nil {
		log.Fatalf("failed to get dataset dimensions: %v", err)
	}
	totalSize := int(dims[0] * dims[1] * dims[2])

	// Create a 1-D slice to hold the data
	flatMatrix := make([]float32, totalSize)

	// Read the data from the dataset into the flat slice
	if err := dataset.Read(&flatMatrix); err != nil {
		log.Fatalf("failed to read dataset: %v", err)
	}

	// Reshape the flat slice into a 3-D slice
	matrix := make([][][]float32, dims[0])
	for i := range matrix {
		matrix[i] = make([][]float32, dims[1])
		for j := range matrix[i] {
			matrix[i][j] = flatMatrix[(i*int(dims[1])*int(dims[2]))+(j*int(dims[2])) : (i*int(dims[1])*int(dims[2]))+(j+1)*int(dims[2])]
		}
	}

	// Print the 3-D matrix (or process it as needed)
	// fmt.Println("last element of 3-D matrix data:", matrix)

	return matrix
}

// GetFormation_energy() takes into a dataset name, and returns formation_energy in the dataset.
// formation_energy, Shape: (4,), Dtype: float64
func GetFormation_energy(fileName, datasetName, attributeName string) []float64 {
	file, err := hdf5.OpenFile(fileName, hdf5.F_ACC_RDONLY)
	if err != nil {
		log.Fatal("Error opening file:", err)
	}
	defer file.Close()

	// access the dataset
	path := datasetName + "/" + attributeName
	dataset, err := file.OpenDataset(path)
	if err != nil {
		log.Fatal("Error opening dataset:", err)
	}
	defer dataset.Close()

	// Get the dimensions of the dataset
	space := dataset.Space()
	dims, _, err := space.SimpleExtentDims()
	if err != nil {
		log.Fatalf("Failed to get dataset dimensions: %v", err)
	}

	// Ensure it's a 1D dataset
	if len(dims) != 1 {
		log.Fatalf("Expected a 1-dimensional dataset, but got %d dimensions", len(dims))
	}

	// Create a slice to hold the data
	formationEnergy := make([]float64, dims[0])

	// Read the data from the dataset into the slice
	if err := dataset.Read(&formationEnergy); err != nil {
		log.Fatalf("Failed to read dataset: %v", err)
	}

	return formationEnergy
}

// GetDft_total_gradient() takes into a dataset name, and returns dft_total_gradient in the dataset.
// dft_total_gradient, Shape: (4, 38, 3), Dtype: float32
func GetDft_total_gradient(fileName, datasetName, attributeName string) [][][]float32 {
	file, err := hdf5.OpenFile(fileName, hdf5.F_ACC_RDONLY)
	if err != nil {
		log.Fatal("Error opening file:", err)
	}
	defer file.Close()

	// access the dataset
	path := datasetName + "/" + attributeName
	dataset, err := file.OpenDataset(path)
	if err != nil {
		log.Fatal("Error opening dataset:", err)
	}
	defer dataset.Close()

	// Get the dimensions of the 3-D matrix
	space := dataset.Space()
	dims, _, err := space.SimpleExtentDims()
	if err != nil {
		log.Fatalf("failed to get dataset dimensions: %v", err)
	}
	totalSize := int(dims[0] * dims[1] * dims[2])

	// Create a 1-D slice to hold the data
	flatMatrix := make([]float32, totalSize)

	// Read the data from the dataset into the flat slice
	if err := dataset.Read(&flatMatrix); err != nil {
		log.Fatalf("failed to read dataset: %v", err)
	}

	// Reshape the flat slice into a 3-D slice
	matrix := make([][][]float32, dims[0])
	for i := range matrix {
		matrix[i] = make([][]float32, dims[1])
		for j := range matrix[i] {
			matrix[i][j] = flatMatrix[(i*int(dims[1])*int(dims[2]))+(j*int(dims[2])) : (i*int(dims[1])*int(dims[2]))+(j+1)*int(dims[2])]
		}
	}

	// Print the 3-D matrix (or process it as needed)
	// fmt.Println("last element of 3-D matrix data:", matrix)

	return matrix
}

// GetAtomic_numbers() takes into a dataset name, and returns atomic_numbers in the dataset.
// atomic_numbers, Shape: (38,), Dtype: int16
func GetAtomic_numbers(fileName, datasetName, attributeName string) []int {

	file, err := hdf5.OpenFile(fileName, hdf5.F_ACC_RDONLY)
	if err != nil {
		log.Fatal("Error opening file:", err)
	}
	defer file.Close()

	// access the dataset
	path := datasetName + "/" + attributeName
	dataset, err := file.OpenDataset(path)
	if err != nil {
		log.Fatal("Error opening dataset:", err)
	}
	defer dataset.Close()

	// Get the dimensions of the dataset
	space := dataset.Space()
	dims, _, err := space.SimpleExtentDims()
	if err != nil {
		log.Fatalf("Failed to get dataset dimensions: %v", err)
	}

	// Ensure it's a 1D dataset
	if len(dims) != 1 {
		log.Fatalf("Expected a 1-dimensional dataset, but got %d dimensions", len(dims))
	}

	// Create a slice to hold the data
	atomicNumbers := make([]int, dims[0])

	// Read the data from the dataset into the slice
	if err := dataset.Read(&atomicNumbers); err != nil {
		log.Fatalf("Failed to read dataset: %v", err)
	}

	return atomicNumbers
}

// GetConformations() takes into a dataset name, and returns conformations in the dataset.
// conformations, Shape: (4, 38, 3), Dtype: float32
func GetConformations(fileName, datasetName, attributeName string) [][][]float32 {
	file, err := hdf5.OpenFile(fileName, hdf5.F_ACC_RDONLY)
	if err != nil {
		log.Fatal("Error opening file:", err)
	}
	defer file.Close()

	// access the dataset
	path := datasetName + "/" + attributeName
	dataset, err := file.OpenDataset(path)
	if err != nil {
		log.Fatal("Error opening dataset:", err)
	}
	defer dataset.Close()

	// Get the dimensions of the 3-D matrix
	space := dataset.Space()
	dims, _, err := space.SimpleExtentDims()
	if err != nil {
		log.Fatalf("failed to get dataset dimensions: %v", err)
	}
	totalSize := int(dims[0] * dims[1] * dims[2])

	// Create a 1-D slice to hold the data
	flatMatrix := make([]float32, totalSize)

	// Read the data from the dataset into the flat slice
	if err := dataset.Read(&flatMatrix); err != nil {
		log.Fatalf("failed to read dataset: %v", err)
	}

	// Reshape the flat slice into a 3-D slice
	matrix := make([][][]float32, dims[0])
	for i := range matrix {
		matrix[i] = make([][]float32, dims[1])
		for j := range matrix[i] {
			matrix[i][j] = flatMatrix[(i*int(dims[1])*int(dims[2]))+(j*int(dims[2])) : (i*int(dims[1])*int(dims[2]))+(j+1)*int(dims[2])]
		}
	}

	// Print the 3-D matrix (or process it as needed)
	// fmt.Println("last element of 3-D matrix data:", matrix)

	return matrix
}

// GetDft_total_energy() takes into a dataset name, and returns dft_total_energy in the dataset.
// dft_total_energy, Shape: (4,), Dtype: float64
func GetDft_total_energy(fileName, datasetName, attributeName string) []float64 {
	file, err := hdf5.OpenFile(fileName, hdf5.F_ACC_RDONLY)
	if err != nil {
		log.Fatal("Error opening file:", err)
	}
	defer file.Close()

	// access the dataset
	path := datasetName + "/" + attributeName
	dataset, err := file.OpenDataset(path)
	if err != nil {
		log.Fatal("Error opening dataset:", err)
	}
	defer dataset.Close()

	// Get the dimensions of the dataset
	space := dataset.Space()
	dims, _, err := space.SimpleExtentDims()
	if err != nil {
		log.Fatalf("Failed to get dataset dimensions: %v", err)
	}

	// Ensure it's a 1D dataset
	if len(dims) != 1 {
		log.Fatalf("Expected a 1-dimensional dataset, but got %d dimensions", len(dims))
	}

	// Create a slice to hold the data
	dftTotalEnergy := make([]float64, dims[0])

	// Read the data from the dataset into the slice
	if err := dataset.Read(&dftTotalEnergy); err != nil {
		log.Fatalf("Failed to read dataset: %v", err)
	}

	return dftTotalEnergy
}

// GetFeatures() takes into a dataset name, and returns all the 14 features in the dataset as output.
// Order of putput features: atomic_numbers, conformations, dft_total_energy, dft_total_gradient, formation_energy, mayer_indices, mbis_charges, mbis_dipoles, mbis_octupoles, mbis_quadrupoles, scf_dipole, scf_quadrupole,smiles, wiberg_lowdin_indices
func GetFeatures(file, datasetName, AttributeName string) ([]int, [][][]float32, []float64, [][][]float32, []float64, [][][]float32, [][][]float32, [][][]float32, [][][][][]float32, [][][][]float32, [][]float32, [][][]float32, string, [][][]float32) {
	atomic_numbers := GetAtomic_numbers(file, datasetName, AttributeName)
	conformations := GetConformations(file, datasetName, AttributeName)
	dft_total_energy := GetDft_total_energy(file, datasetName, AttributeName)
	dft_total_gradient := GetDft_total_gradient(file, datasetName, AttributeName)
	formation_energy := GetFormation_energy(file, datasetName, AttributeName)
	mayer_indices := GetMayer_indices(file, datasetName, AttributeName)
	mbis_charges := GetMbis_charges(file, datasetName, AttributeName)
	mbis_dipoles := GetMbis_dipoles(file, datasetName, AttributeName)
	mbis_octupoles := GetMbis_octupoles(file, datasetName, AttributeName)
	mbis_quadrupoles := GetMbis_quadrupoles(file, datasetName, AttributeName)
	scf_dipole := GetScf_dipole(file, datasetName, AttributeName)
	scf_quadrupole := GetScf_quadrupole(file, datasetName, AttributeName)
	smiles := GetSmiles(file, datasetName, AttributeName)
	wiberg_lowdin_indices := GetWiberg_lowdin_indices(file, datasetName, AttributeName)

	return atomic_numbers, conformations, dft_total_energy, dft_total_gradient, formation_energy, mayer_indices, mbis_charges, mbis_dipoles, mbis_octupoles, mbis_quadrupoles, scf_dipole, scf_quadrupole, smiles, wiberg_lowdin_indices

}

func main() {

	// fileName := "pubchem-1-2500.hdf5"
	// attributeName := "subset"

	// categories, err := categorizeDatasets(fileName, attributeName)
	// if err != nil {
	// 	log.Fatalf("Failed to categorize datasets: %v", err)
	// }

	// // Print the categorized datasets
	// for category, datasets := range categories {
	// 	fmt.Printf("Category: %s\n", category)
	// 	for _, dataset := range datasets {
	// 		fmt.Printf("  - Dataset: %s\n", dataset)
	// 	}
	// }

}
