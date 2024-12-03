package main

import (
	"math"
	"math/rand"
)

//Monte Carlo simulation
//try different ligand positions relative to protein
//calculate the energy, accept the new position or not
// can't have simple minima fn as can then get stuck in a local minima situation : "Metropolis criteria"
//also takes in a iterations parameter to measure
//more likely to find global minima if more iterations but computationally expensive
// try using Go's concurrency to make it faster

// Perturb the ligand by randomly changing its x, y and z position

//Metropolis criteria: If energy reduced, accept it; else if energy >, accept based on the probability given by e^(-deltaE / temperature)
// where deltaE = newEnergy - oldEnergy
//design choice: Take temp as a parameter or hard code (Preference take as a parameter and set it globally to body temp for now or with default value)

// Energy function
//Lennard-Jones potential or Coulomb interactions.
//Takes in a protein and a ligand
// Goes over all the protein's atoms and all the ligand's atoms and calculates the energy. (Coulomb: k q1, q2/r^2)
// returns an energy list correponding to energy of each ligand
//, evaluate each ligand

// Metropolis acceptance criterion
func AcceptMove(currentEnergy, newEnergy, temperature float64) bool {
	if newEnergy < currentEnergy {
		// fmt.Println("Identifying as less:", newEnergy)
		return true
	}
	// fmt.Println("Identifying as more")
	// fmt.Println(newEnergy)
	deltaE := newEnergy - currentEnergy
	probability := math.Exp(-deltaE / temperature)
	return rand.Float64() < probability
}

func PerturbLigand(ligand Molecule) Molecule {
	for i := range ligand.atoms {
		ligand.atoms[i].Position.X += (rand.Float64() - 0.5) * 0.1
		ligand.atoms[i].Position.Y += (rand.Float64() - 0.5) * 0.1
		ligand.atoms[i].Position.Z += (rand.Float64() - 0.5) * 0.1
	}
	return ligand
}

func CalculateEnergy(protein, ligand Molecule) float64 {
	energy := 0.0
	for _, atomP := range protein.atoms {
		for _, atomL := range ligand.atoms {
			atomPPos := atomP.Position
			atomLPos := atomL.Position
			distance := math.Sqrt(math.Pow(atomPPos.X-atomLPos.X, 2) + math.Pow(atomPPos.Y-atomLPos.Y, 2) + math.Pow(atomPPos.Z-atomLPos.Z, 2))
			energy += K * (atomP.Charge * atomL.Charge) / distance //kQ1,Q2/d^2
		}
	}
	return energy
}

func SimulateEnergyMinimization(protein, ligand Molecule, iterations int, temperature float64) Molecule {
	currentLigand := ligand
	currentEnergy := CalculateEnergy(protein, currentLigand)

	for i := 0; i < iterations; i++ {
		newLigand := PerturbLigand(currentLigand)
		newEnergy := CalculateEnergy(protein, newLigand)
		if AcceptMove(currentEnergy, newEnergy, temperature) {
			currentLigand = newLigand
			currentEnergy = newEnergy
		}
	}
	return currentLigand
}

func SimulateEnergyMinimizationParallel(protein, ligand Molecule, iterations int, temperature float64, numProcs int) Molecule {
	currentLigand := ligand
	currentEnergy := CalculateEnergy(protein, currentLigand)
	width := iterations / numProcs
	c := make(chan Molecule)
	for i := 0; i < numProcs; i++ {
		go SimulateEnergyMinimizationOneProc(protein, currentLigand, width, temperature, c)
	}
	for i := 0; i < numProcs; i++ {
		newLigand := <-c
		newEnergy := CalculateEnergy(protein, newLigand)
		//fmt.Println(newEnergy)
		if AcceptMove(currentEnergy, newEnergy, temperature) {
			currentLigand = newLigand
			currentEnergy = newEnergy
		}
	}
	return currentLigand
}

func SimulateEnergyMinimizationOneProc(protein, ligand Molecule, iterations int, temperature float64, c chan Molecule) {
	currentLigand := ligand
	currentEnergy := CalculateEnergy(protein, currentLigand)
	for i := 0; i < iterations; i++ {
		newLigand := PerturbLigand(currentLigand)
		newEnergy := CalculateEnergy(protein, newLigand)
		if AcceptMove(currentEnergy, newEnergy, temperature) {
			currentLigand = newLigand
			currentEnergy = newEnergy
		}
	}
	c <- currentLigand
}

func SimulateMultipleLigands(protein Molecule, ligands []Molecule, iterations int, temperature float64, numProcs int) ([]Molecule, []float64) {
	minEnergy := make([]float64, len(ligands))
	minLigands := make([]Molecule, len(ligands))
	for i, ligand := range ligands {
		minLigands[i] = SimulateEnergyMinimizationParallel(protein, ligand, iterations, temperature, numProcs)
		minEnergy[i] = CalculateEnergy(protein, minLigands[i])
	}
	return minLigands, minEnergy
}

func SimulateMultipleLigandsParallel(protein Molecule, ligands []Molecule, iterations int, temperature float64, numProcs int) ([]Molecule, []float64) {
	minEnergy := make([]float64, 0)
	minLigands := make([]Molecule, 0)
	ligandChannels := make([]chan []Molecule, numProcs)
	energyChannels := make([]chan []float64, numProcs)
	for i := range ligandChannels {
		ligandChannels[i] = make(chan []Molecule, len(ligands))
		energyChannels[i] = make(chan []float64, len(ligands))
	}
	width := len(ligands) / numProcs
	for i := 0; i < numProcs; i++ {
		startIndex := i * width
		var endIndex int
		if i != numProcs-1 {
			endIndex = startIndex + width
		} else {
			endIndex = len(ligands)
		}
		go SimulateLigandMinimizationOneProc(protein, ligands[startIndex:endIndex], iterations, temperature, ligandChannels[i], energyChannels[i], numProcs)
	}
	for i := 0; i < numProcs; i++ {
		minEnergy = append(minEnergy, <-energyChannels[i]...)
		minLigands = append(minLigands, <-ligandChannels[i]...)
	}
	return minLigands, minEnergy
}

func SimulateLigandMinimizationOneProc(protein Molecule, ligands []Molecule, iterations int, temperature float64, ligandChannel chan []Molecule, energyChannel chan []float64, numProcs int) ([]Molecule, []float64) {
	minEnergy := make([]float64, len(ligands))
	minLigands := make([]Molecule, len(ligands))
	for i, ligand := range ligands {
		minLigands[i] = SimulateEnergyMinimizationParallel(protein, ligand, iterations, temperature, numProcs)
		minEnergy[i] = CalculateEnergy(protein, minLigands[i])
	}
	ligandChannel <- minLigands
	energyChannel <- minEnergy
	return minLigands, minEnergy
}
