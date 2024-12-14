package main

import (
	"fmt"
	"math"
	"math/rand"
	"os"
	"strconv"
)

// MultipleProteinRMSD computes the RMSD for multiple proteins
// Input: a string dir, an int iterations, a bool rotate, an int numProteins, an int numProcs
// Output: none (prints the average RMSD and generates an RMSD curve plot)
func MultipleProteinRMSD(dir string, iterations int, rotate bool, numProteins int, numProcs int) {
	proteinFiles, err := findFilesWithSubstring(dir, "protein")
	proteinFiles = proteinFiles[50 : numProteins+50]
	proteinLabels := make([]string, len(proteinFiles))
	Check(err)
	rmsd := make([]float64, len(proteinFiles))
	for i := range proteinFiles {
		protein, err := ParseMol2(proteinFiles[i])
		Check(err)
		label := ExtractFileLabel(proteinFiles[i])
		proteinLabels[i] = label
		ligand, err2 := ParseMol2(dir + "/" + label + "_ligand.mol2")
		//ligand = RandomizeLigandPose(ligand)
		Check(err2)
		rmsd[i] = CompareRMSD(protein, ligand, iterations, rotate, TEMPERATURE, numProcs)
	}
	fmt.Println("The average RMSD value was:", average(rmsd))
	outputDir := "Output/rmsd_curve/"
	err2 := os.MkdirAll(outputDir, 0755)
	Check(err2)
	plotRMSD(proteinLabels, rmsd, outputDir+"rmsd_curve_"+strconv.Itoa(len(proteinFiles)))
}

// CompareRMSD simulates energy minimization of the ligand, then calculates the RMSD between the minimized and reference ligand positions.
// Input: a Molecule protein, a Molecule ligand, an int iterations, a bool rotate, a float64 temperature, an int numProcs
// Output: a float64 RMSD value
func CompareRMSD(protein Molecule, ligand Molecule, iterations int, rotate bool, temperature float64, numProcs int) float64 {
	reference := CopyLigand(ligand)
	ligand = RandomizeLigandPose(ligand)
	simulated := SimulateEnergyMinimizationParallel(protein, ligand, iterations, rotate, temperature, numProcs)
	return CalculateRMSD(simulated, reference)
}

// CalculateRMSD  calculates the root-mean-square deviation (RMSD) between corresponding atoms in the two given molecules
// Input: Molecules simulated and reference
// Output: a float64 RMSD value
func CalculateRMSD(simulated, reference Molecule) float64 {
	var sumSquaredDist float64
	numAtoms := float64(len(simulated.atoms))

	for i := 0; i < len(simulated.atoms); i++ {
		simPos := simulated.atoms[i].Position
		refPos := reference.atoms[i].Position

		dx := simPos.X - refPos.X
		dy := simPos.Y - refPos.Y
		dz := simPos.Z - refPos.Z

		sumSquaredDist += dx*dx + dy*dy + dz*dz
	}
	return math.Sqrt(sumSquaredDist / numAtoms)
}

// RandomizeLigandPose applies random translations and rotations to a ligand molecule to generate a randomized starting pose.
// Input: a Molecule ligand
// Output: a Molecule with randomized pose
func RandomizeLigandPose(ligand Molecule) Molecule {
	for i := range ligand.atoms {
		// Apply random translation (±5 Å)
		ligand.atoms[i].Position.X += (rand.Float64() - 0.5) * 10.0
		ligand.atoms[i].Position.Y += (rand.Float64() - 0.5) * 10.0
		ligand.atoms[i].Position.Z += (rand.Float64() - 0.5) * 10.0
	}
	// Apply random rotation
	return RotateLigand(ligand, math.Pi)
}

// average calculates and returns the average of a slice of float64 values
// Input: a slice of float64 values
// Output: a float64 average value
func average(arr []float64) float64 {
	if len(arr) == 0 {
		return 0
	}
	var sum float64
	for _, value := range arr {
		sum += value
	}
	return sum / float64(len(arr))
}
