package main

import (
	"bufio"
	"os"
	"strconv"
	"strings"
)

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
