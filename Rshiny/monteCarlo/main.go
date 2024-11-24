package main

// Here, I changed the num of iterations to make the process faster for APP!
import (
	"encoding/csv"
	"fmt"
	"log"
	"os"
	"runtime"
	"strconv"
)

func main() {
	if len(os.Args) < 3 {
		fmt.Println("Please provide the genome file path.")
		return
	}

	proteinFilePath := os.Args[1]
	ligandFilePaths := os.Args[2:] // All arguments after the first are ligand file paths

	protein, err2 := ParseMol2(proteinFilePath)
	Check(err2)
	fmt.Println("Protein atom count:", len(protein.atoms))

	results := make([]string, len(ligandFilePaths))

	for i, ligandFilePath := range ligandFilePaths {
		ligand, err := ParseMol2(ligandFilePath)
		Check(err)
		fmt.Println("Ligand atom count:", len(ligand.atoms))

		iterations := 3
		temperature := 310.15
		numProcs := runtime.NumCPU()

		// Perform energy minimization
		newLigand := SimulateEnergyMinimizationParallel(protein, ligand, iterations, temperature, numProcs)
		newEnergy := CalculateEnergy(protein, newLigand)

		// Convert energy to a formatted string
		floatStr := strconv.FormatFloat(newEnergy, 'f', 6, 64)

		results[i] = floatStr

	}

	// Create a CSV file
	outFile, err := os.Create("../output/simulations.csv")
	if err != nil {
		log.Fatalf("Failed to create output file: %v", err)
	}
	defer outFile.Close()

	// Create a CSV writer
	writer := csv.NewWriter(outFile)
	defer writer.Flush()

	// Write the header
	header := []string{"BindingEnergy"}
	if err := writer.Write(header); err != nil {
		log.Fatalf("Failed to write header: %v", err)
	}

	// Write the map data
	for _, energyStr := range results {
		row := []string{energyStr}
		if err := writer.Write(row); err != nil {
			log.Fatalf("Failed to write row: %v", err)
		}
	}

	fmt.Println("Binding Energy data written to simulations.csv")
}

// iterations := 3
// // I changed the num of iterations to make the process faster for APP!
// temperature := 310.15
// numProcs := runtime.NumCPU()
// // fmt.Println("Old Energy", CalculateEnergy(protein, ligand))
// //newLigand := SimulateEnergyMinimization(protein, ligand, iterations, temperature)
// newLigand2 := SimulateEnergyMinimizationParallel(protein, ligand, iterations, temperature, numProcs)
// //fmt.Println("New Energy", CalculateEnergy(protein, newLigand))
// fmt.Println("New Energy Parallel", CalculateEnergy(protein, newLigand2))

func Check(err error) {
	if err != nil {
		panic(err)
	}
}
