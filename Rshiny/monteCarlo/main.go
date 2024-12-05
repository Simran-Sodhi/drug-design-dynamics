package main

import (
	"encoding/csv"
	"fmt"
	"io"
	"log"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
)

func main() {
	if len(os.Args) < 3 {
		fmt.Println("Please provide right inputs.")
		return
	}

	proteinFilePath := os.Args[1]
	ligandFilePaths := os.Args[2:] // All arguments after the first are ligand file paths

	protein, err2 := ParseMol2(proteinFilePath)
	Check(err2)

	results := make([]string, len(ligandFilePaths))

	minEnergy := 1e20     // Set an initial large value for comparison
	minEnergyLigand := "" // To store the path of the ligand with minimum energy

	for i, ligandFilePath := range ligandFilePaths {
		ligand, err := ParseMol2(ligandFilePath)
		Check(err)

		iterations := 3000
		temperature := 310.15
		numProcs := runtime.NumCPU()

		// Perform energy minimization
		newLigand := SimulateEnergyMinimizationParallel(protein, ligand, iterations, temperature, numProcs)
		newEnergy := CalculateEnergy(protein, newLigand)

		// Update the minimum energy and ligand file path
		if newEnergy < minEnergy {
			minEnergy = newEnergy
			minEnergyLigand = ligandFilePath
		}

		// Convert energy to a formatted string
		floatStr := strconv.FormatFloat(newEnergy, 'f', 6, 64)

		results[i] = floatStr
		fmt.Println("Finish i ligand here!:", i)

	}

	// Copy the ligand with the minimum energy to the output directory
	if minEnergyLigand != "" {
		outputDir := "../output"
		baseName := filepath.Base(minEnergyLigand) // Get the original file name
		outputLigandPath := filepath.Join(outputDir, baseName)

		err := CopyFile(minEnergyLigand, outputLigandPath)
		if err != nil {
			log.Fatalf("Failed to copy the ligand file: %v", err)
		}

		fmt.Println("Ligand with minimum energy copied to:", outputLigandPath)
	}

	// Copy the protein to the output directory
	outputDir := "../output"
	// fmt.Println("Protein file path:", proteinFilePath)
	baseName := filepath.Base(proteinFilePath) // Get the original file name
	outputProteinPath := filepath.Join(outputDir, baseName)

	err := CopyFile(proteinFilePath, outputProteinPath)
	if err != nil {
		log.Fatalf("Failed to copy the protein file: %v", err)
	}

	fmt.Println("Protein with minimum energy copied to:", outputProteinPath)

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

func Check(err error) {
	if err != nil {
		panic(err)
	}
}

// CopyFile copies a file from src to dst
func CopyFile(src, dst string) error {
	sourceFile, err := os.Open(src)
	if err != nil {
		return err
	}
	defer sourceFile.Close()

	destFile, err := os.Create(dst)
	if err != nil {
		return err
	}
	defer destFile.Close()

	_, err = io.Copy(destFile, sourceFile)
	if err != nil {
		return err
	}

	return nil
}
