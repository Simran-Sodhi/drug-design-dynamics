package main

import (
	"bufio"
	"fmt"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

// ParsePDB parses a PDB file to extract atomic coordinates and charges, returning a Molecule containing the parsed atoms.
// Input: a string filename
// Output: a Molecule and an error
func ParsePDB(filename string) (Molecule, error) {
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

// ParseMol2 parses a MOL2 file, extracting atomic information from the ATOM section, including coordinates and charges, and returns a Molecule.
// Input: a string filename
// Output: a Molecule and an error
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

// SaveToMol2 saves the current Molecule to a MOL2 file, writing its atoms with default properties such as type and substructure information.
// Input: a string filename
// Output: an error or nil
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

// findFilesWithSubstring searches for files in the specified root directory that contain the given substring in their filenames and returns a list of matching file paths.
// Input: a string rootDir, a string searchString
// Output: a slice of strings containing matching file paths, and an error or nil
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
