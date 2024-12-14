package main

import (
	"bufio"
	"encoding/csv"
	"fmt"
	"log"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
)

// plotRMSD plots RMSD values for various proteins, saves the plot as a PNG, and writes corresponding indices to a CSV file.
// Input: a slice of strings x, a slice of float64 y, a string fileName
// Output: none (saves a plot and CSV file)
func plotRMSD(x []string, y []float64, fileName string) {
	plotXY(x, y, fileName, "RMSD for various proteins", "Protein_index", "RMSD Value")
}

// plotEnergy plots ligand binding energy values, saves the plot as a PNG, and writes corresponding indices to a CSV file.
// Input: a slice of strings x, a slice of float64 y, a string fileName
// Output: none (saves a plot and CSV file)
func plotEnergy(x []string, y []float64, fileName string) {
	plotXY(x, y, fileName, "Energy of various ligands", "Ligand_index", "Protein Ligand Binding Energy")
}

// plotXY creates a plot using provided data and labels and saves it as a PNG
// Input: a slice of strings x, a slice of float64 y, a string fileName, a string title, a string xLabel, a string yLabel
// Output: none (saves a plot and CSV file)
func plotXY(x []string, y []float64, fileName string, title, xLabel, yLabel string) {
	points := make(plotter.XYs, len(x))
	for i := range x {
		points[i].X = float64(i) // Numeric representation of X
		points[i].Y = y[i]
	}
	x_ind := make([]string, len(x))
	for i := range x {
		x_ind[i] = strconv.Itoa(i)
	}
	// Create a new plot
	p := plot.New()

	// Set plot title and axis labels
	p.Title.Text = title
	p.X.Label.Text = xLabel
	p.Y.Label.Text = yLabel

	// Add a line to the plot
	line, err := plotter.NewLine(points)
	Check(err)
	p.Add(line)

	// Customize the X-axis with string labels
	p.NominalX(x_ind...)

	// Save the plot to a file
	err2 := p.Save(6*vg.Inch, 4*vg.Inch, fileName+".png")
	Check(err2)
	log.Printf("Plot saved as %s.png", fileName)
	saveMappingToCSV(fileName+".csv", x)
}

// saveMappingToCSV saves the mapping of plot indices to labels in a CSV file.
// Input: a string fileName, a slice of strings labels
// Output: none (saves a CSV file)
func saveMappingToCSV(fileName string, labels []string) {
	file, err := os.Create(fileName)
	Check(err)
	defer file.Close()

	writer := csv.NewWriter(file)
	defer writer.Flush()

	// Write header
	err2 := writer.Write([]string{"Index", "File Name"})
	Check(err2)

	// Write mapping rows
	for i, label := range labels {
		if err := writer.Write([]string{strconv.Itoa(i), label}); err != nil {
			log.Fatalf("Failed to write row: %v", err)
		}
	}
	log.Printf("Mapping of plot indices to labels saved to %s", fileName)
}

// ExtractFileLabel extracts the ligand label from a file path by removing the file extension and prefix.
// Input: a string filePath
// Output: a string file label
func ExtractFileLabel(filePath string) string {
	// Extract the base file name
	base := filepath.Base(filePath)
	// Split the base name by "_" and remove the extension
	parts := strings.Split(base, "_")
	if len(parts) > 0 {
		return strings.TrimSuffix(parts[0], ".mol2")
	}
	return base
}

// SaveMinimumEnergyLigand identifies the ligand with minimum energy and saves its structure to a specified MOL2 file.
// Input: a slice of float64 energyList, a slice of strings ligandFiles, a string fileName, a slice of Molecule minLigands
// Output: none (saves the best ligand structure)
func SaveMinimumEnergyLigand(energyList []float64, ligandFiles []string, fileName string, minLigands []Molecule) {
	minIndex := 0
	minEnergy := 0.0
	for i, energy := range energyList {
		if energy < minEnergy {
			minIndex = i
			minEnergy = energy
		}
	}
	err4 := UpdateMol2Coordinates(ligandFiles[minIndex], fileName, minLigands[minIndex])
	Check(err4)
	fmt.Printf("Wrote min ligand file to %s", fileName)
}

// UpdateMol2Coordinates updates a MOL2 file with new atomic coordinates from a specified Molecule and saves the updated file.
// Input: a string originalFile, a string updatedFile, a Molecule newMolecule
// Output: an error or nil
func UpdateMol2Coordinates(originalFile, updatedFile string, newMolecule Molecule) error {
	// Open the original file
	file, err := os.Open(originalFile)
	if err != nil {
		return fmt.Errorf("failed to open file: %v", err)
	}
	defer file.Close()

	// Create the updated file
	output, err := os.Create(updatedFile)
	if err != nil {
		return fmt.Errorf("failed to create output file: %v", err)
	}
	defer output.Close()

	scanner := bufio.NewScanner(file)
	writer := bufio.NewWriter(output)

	atomIndex := 0
	inAtomSection := false

	for scanner.Scan() {
		line := scanner.Text()

		// Check if we're in the ATOM section
		if strings.HasPrefix(line, "@<TRIPOS>ATOM") {
			inAtomSection = true
			writer.WriteString(line + "\n")
			continue
		} else if strings.HasPrefix(line, "@<TRIPOS>") && inAtomSection {
			// Exit ATOM section
			inAtomSection = false
		}

		if inAtomSection && atomIndex < len(newMolecule.atoms) {
			// Update the atom line with new coordinates
			parts := strings.Fields(line)
			if len(parts) >= 9 {
				atom := newMolecule.atoms[atomIndex]
				parts[2] = fmt.Sprintf("%.4f", atom.Position.X) // X coordinate
				parts[3] = fmt.Sprintf("%.4f", atom.Position.Y) // Y coordinate
				parts[4] = fmt.Sprintf("%.4f", atom.Position.Z) // Z coordinate
				writer.WriteString(strings.Join(parts, " ") + "\n")
				atomIndex++
			} else {
				writer.WriteString(line + "\n")
			}
		} else {
			// Copy other lines as is
			writer.WriteString(line + "\n")
		}
	}

	if err := scanner.Err(); err != nil {
		return fmt.Errorf("error reading file: %v", err)
	}

	if atomIndex != len(newMolecule.atoms) {
		return fmt.Errorf("mismatch between number of atoms in molecule and file")
	}

	writer.Flush()
	return nil
}
