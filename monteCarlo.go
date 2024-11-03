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
		return true
	}
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
			distance := math.Sqrt(math.Pow(atomP.X-atomL.X, 2) + math.Pow(atomP.Y-atomL.Y, 2) + math.Pow(atomP.Z-atomL.Z, 2))
			energy += K * (atomP.Charge * atomL.Charge) / distance //kQ1,Q2/d^2
		}
	}
	return energy
}

func MonteCarloSimulation(protein, ligand Molecule, iterations int, temperature float64) Molecule {
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
