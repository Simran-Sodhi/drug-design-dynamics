package main

import (
	"bufio"
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

func parsePDB(filename string) (Molecule, error) {
	file, err := os.Open(filename)
	molecule := Molecule{}
	if err != nil {
		return molecule, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, "ATOM") || strings.HasPrefix(line, "HETATM") {
			x, _ := strconv.ParseFloat(strings.TrimSpace(line[30:38]), 64)
			y, _ := strconv.ParseFloat(strings.TrimSpace(line[38:46]), 64)
			z, _ := strconv.ParseFloat(strings.TrimSpace(line[46:54]), 64)
			charge, _ := strconv.ParseFloat(strings.TrimSpace(line[78:]), 64) // Assuming charge is appended at the end

			atom := Atom{
				Position: Position3d{X: x, Y: y, Z: z},
				Charge:   charge,
			}
			molecule.atoms = append(molecule.atoms, atom)
		}
	}

	if err := scanner.Err(); err != nil {
		return Molecule{}, err
	}
	return molecule, nil
}

func ParseMol2(filename string) (Molecule, error) {
	file, err := os.Open(filename)
	molecule := Molecule{}
	if err != nil {
		return molecule, err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	inAtomSection := false

	for scanner.Scan() {
		line := scanner.Text()

		// Check for the ATOM section start
		if strings.HasPrefix(line, "@<TRIPOS>ATOM") {
			inAtomSection = true
			continue
		}

		// Check for the end of the ATOM section
		if inAtomSection && strings.HasPrefix(line, "@<TRIPOS>") {
			inAtomSection = false
			break
		}

		// Parse ATOM lines
		if inAtomSection {
			fields := strings.Fields(line)
			if len(fields) < 9 {
				continue
			}

			// Extract atom details
			//id, _ := strconv.Atoi(fields[0])                           // Atom ID
			x, _ := strconv.ParseFloat(fields[2], 64)                  // X coordinate
			y, _ := strconv.ParseFloat(fields[3], 64)                  // Y coordinate
			z, _ := strconv.ParseFloat(fields[4], 64)                  // Z coordinate
			charge, _ := strconv.ParseFloat(fields[len(fields)-1], 64) // Partial charge

			atom := Atom{
				//ID:       id,
				//Name:     fields[1], // Atom Name
				Position: Position3d{X: x, Y: y, Z: z},
				//Type:     fields[5], // Atom Type
				Charge: charge,
			}
			molecule.atoms = append(molecule.atoms, atom)
		}
	}

	if err := scanner.Err(); err != nil {
		return Molecule{}, err
	}

	return molecule, nil
}

//func getProteinAndLigands(id string)

func (m *Molecule) SaveToMol2(filename string) error {
	file, err := os.Create(filename)
	Check(err)
	defer file.Close()

	// Write MOL2 header
	fmt.Fprintf(file, "@<TRIPOS>MOLECULE\n")
	fmt.Fprintf(file, "Generated Molecule\n")
	fmt.Fprintf(file, "%d\n", len(m.atoms)) // Number of atoms
	fmt.Fprintf(file, "SMALL\n")
	fmt.Fprintf(file, "USER_CHARGES\n\n")

	// Write ATOM section
	fmt.Fprintf(file, "@<TRIPOS>ATOM\n")
	for i, atom := range m.atoms {
		fmt.Fprintf(file, "%6d %s %10.4f %10.4f %10.4f %s %4d %s %10.4f\n",
			i+1,                                               // Atom ID
			"C",                                               // Atom name (default to "C" for simplicity)
			atom.Position.X, atom.Position.Y, atom.Position.Z, // Coordinates
			"C.3",       // Atom type (default to "C.3" for simplicity)
			1,           // Substructure ID (default to 1)
			"MOLECULE",  // Substructure name
			atom.Charge) // Charge
	}

	return nil
}

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

func plotEnergy(x []string, y []float64) {

	points := make(plotter.XYs, len(x))
	for i := range x {
		points[i].X = float64(i) // Numeric representation of X
		points[i].Y = y[i]
	}

	// Create a new plot
	p := plot.New()

	// Set plot title and axis labels
	p.Title.Text = "Energy of various ligands"
	p.Y.Label.Text = "Protein Ligand Binding Energy"

	// Add a line to the plot
	line, err := plotter.NewLine(points)
	Check(err)
	p.Add(line)

	// Customize the X-axis with string labels
	p.NominalX(x...)

	// Save the plot to a file
	err2 := p.Save(6*vg.Inch, 4*vg.Inch, "line_plot_strings.png")
	Check(err2)

	log.Println("Line plot saved as line_plot_strings.png")
}

func findFilesWithSubstring(rootDir, searchString string) ([]string, error) {
	var matchingFiles []string

	// Walk through the directory
	err := filepath.Walk(rootDir, func(path string, info os.FileInfo, err error) error {
		if err != nil {
			return err
		}
		// Check if it's a file and contains the search string
		if !info.IsDir() && strings.Contains(info.Name(), searchString) {
			matchingFiles = append(matchingFiles, path)
		}
		return nil
	})
	if err != nil {
		return nil, err
	}
	return matchingFiles, nil
}
