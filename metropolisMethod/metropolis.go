package main

import (
	"math"
	"math/rand"
)

// SimulateMultipleLigands simulates energy minimization for multiple ligands sequentially by calling the energy minimization function for each ligand.
// Input: a Molecule protein, a slice of Molecule ligands, an int iterations, a float64 temperature, an int numProcs
// Output: a slice of minimized Molecule ligands and corresponding float64 energies
func SimulateMultipleLigands(protein Molecule, ligands []Molecule, iterations int, rotate bool, temperature float64, numProcs int) ([]Molecule, []float64) {
	minEnergy := make([]float64, len(ligands))
	minLigands := make([]Molecule, len(ligands))
	for i, ligand := range ligands {
		minLigands[i] = SimulateEnergyMinimizationParallel(protein, ligand, iterations, rotate, temperature, numProcs)
		minEnergy[i] = CalculateEnergy(protein, minLigands[i])
	}
	return minLigands, minEnergy
}

// SimulateMultipleLigandsParallel simulates energy minimization for multiple ligands distributing ligands across processors.
// Input: a Molecule protein, a slice of Molecule ligands, an int iterations, a float64 temperature, an int numProcs
// Output: a slice of minimized Molecule ligands and corresponding float64 energies, calculated after having distributed them over numProcs
func SimulateMultipleLigandsParallel(protein Molecule, ligands []Molecule, iterations int, rotate bool, temperature float64, numProcs int) ([]Molecule, []float64) {
	minEnergy := make([]float64, 0)
	minLigands := make([]Molecule, 0)
	ligandChannels := make([]chan MultipleLigandSimulationOutput, numProcs)
	for i := range ligandChannels {
		ligandChannels[i] = make(chan MultipleLigandSimulationOutput, len(ligands))
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
		go SimulateLigandMinimizationOneProc(protein, ligands[startIndex:endIndex], iterations, rotate, temperature, numProcs, ligandChannels[i])
	}
	for i := 0; i < numProcs; i++ {
		minLigAndDelta := <-ligandChannels[i]
		minEnergy = append(minEnergy, minLigAndDelta.Energy...)
		minLigands = append(minLigands, minLigAndDelta.Ligand...)
	}
	return minLigands, minEnergy
}

// SimulateLigandMinimizationOneProc minimizes ligand energies in a single processor and sends results through a channel.
// Input: a Molecule protein, a slice of Molecule ligands, an int iterations, a float64 temperature, an int numProcs, a channel ligandChannel
// Output: none (but sends the minimized ligands and their energies sent through the channel ligandChannel)
func SimulateLigandMinimizationOneProc(protein Molecule, ligands []Molecule, iterations int, rotate bool, temperature float64, numProcs int, ligandChannel chan MultipleLigandSimulationOutput) {
	minEnergy := make([]float64, len(ligands))
	minLigands := make([]Molecule, len(ligands))
	for i, ligand := range ligands {
		minLigands[i] = SimulateEnergyMinimizationParallel(protein, ligand, iterations, rotate, temperature, numProcs)
		minEnergy[i] = CalculateEnergy(protein, minLigands[i])
	}
	ligandChannel <- MultipleLigandSimulationOutput{
		Ligand: minLigands,
		Energy: minEnergy,
	}
}

// SimulateEnergyMinimization performs energy minimization using the Metropolis criterion
// Input: a Molecule protein, a Molecule ligand, an int iterations, a float64 temperature
// Output: a minimized Molecule ligand
func SimulateEnergyMinimization(protein, ligand Molecule, iterations int, rotate bool, temperature float64) Molecule {
	currentLigand := ligand
	currentEnergy := CalculateEnergy(protein, currentLigand)
	minDistance := 1.5
	maxAngle := 0.79 //radians = ~45 degrees
	for i := 0; i < iterations; i++ {
		var newLigand Molecule
		if rotate {
			newLigand = JitterAndRotateLigand(currentLigand, minDistance, maxAngle)
		} else {
			newLigand = JitterLigand(currentLigand, minDistance)
		}
		newEnergy := CalculateEnergy(protein, newLigand)
		if AcceptMove(currentEnergy, newEnergy, temperature) {
			currentLigand = newLigand
			currentEnergy = newEnergy
		}
	}
	return currentLigand
}

// SimulateEnergyMinimizationParallel performs energy minimization using the Metropolis criterion distributed over processors
// Input: a Molecule protein, a Molecule ligand, an int iterations, a float64 temperature, an int numProcs
// Output: a minimized Molecule ligand
func SimulateEnergyMinimizationParallel(protein, ligand Molecule, iterations int, rotate bool, temperature float64, numProcs int) Molecule {
	currentLigand := ligand
	currentEnergy := CalculateEnergy(protein, currentLigand)
	width := iterations / numProcs
	c := make(chan Molecule)
	for i := 0; i < numProcs; i++ {
		go SimulateEnergyMinimizationOneProc(protein, currentLigand, width, rotate, temperature, c)
	}
	for i := 0; i < numProcs; i++ {
		newLigand := <-c
		newEnergy := CalculateEnergy(protein, newLigand)
		if AcceptMove(currentEnergy, newEnergy, temperature) {
			currentLigand = newLigand
			currentEnergy = newEnergy
		}
	}
	return currentLigand
}

// SimulateEnergyMinimizationOneProc minimizes energy of a protein ligand interaction and sends the minimized ligand through a channel
// Input: a Molecule protein, a Molecule ligand, an int iterations, a float64 temperature, a channel c
// Output: none (sends the minimized ligand results through channel c)
func SimulateEnergyMinimizationOneProc(protein, ligand Molecule, iterations int, rotate bool, temperature float64, c chan Molecule) {
	currentLigand := ligand
	currentEnergy := CalculateEnergy(protein, currentLigand)
	minDistance := 1.5
	maxAngle := 0.79 //radians = ~45 degrees
	for i := 0; i < iterations; i++ {
		var newLigand Molecule
		if rotate {
			newLigand = JitterAndRotateLigand(currentLigand, minDistance, maxAngle)
		} else {
			newLigand = JitterLigand(currentLigand, minDistance)
		}
		newEnergy := CalculateEnergy(protein, newLigand)
		if AcceptMove(currentEnergy, newEnergy, temperature) {
			currentLigand = newLigand
			currentEnergy = newEnergy
		}
	}
	c <- currentLigand
}

// AcceptMove determines whether to accept a new ligand state based on the Metropolis criterion.
// Input: two float64 values for current and new energy, and a float64 temperature
// Output: a bool indicating acceptance
func AcceptMove(currentEnergy, newEnergy, temperature float64) bool {
	if newEnergy < currentEnergy {
		return true
	}
	deltaE := newEnergy - currentEnergy
	// Metropolis acceptance criterion
	probability := math.Exp(-deltaE / temperature)
	return rand.Float64() < probability
}

// JitterLigand applies random changes to ligand atom positions while ensuring no atom collisions.
// Input: a Molecule ligand, a float64 minDistance
// Output: a perturbed Molecule ligand which follows no atom collision based on the provided minDistance
func JitterLigand(ligand Molecule, minDistance float64) Molecule {
	for {
		newLigand := ligand
		for i := range newLigand.atoms {
			newLigand.atoms[i].Position.X += (rand.Float64() - 0.5) * 0.1
			newLigand.atoms[i].Position.Y += (rand.Float64() - 0.5) * 0.1
			newLigand.atoms[i].Position.Z += (rand.Float64() - 0.5) * 0.1
		}
		if isCollisionFree(newLigand, minDistance) {
			return newLigand
		}

	}
}

// JitterAndRotateLigand applies random changes to ligand atom positions and rotates it while ensuring no atom collisions.
// Input: a Molecule ligand, a float64 minDistance
// Output: a perturbed Molecule ligand which follows no atom collision based on the provided minDistance
func JitterAndRotateLigand(ligand Molecule, minDistance float64, maxAngle float64) Molecule {
	for {
		newLigand := RotateLigand(ligand, maxAngle)
		for i := range newLigand.atoms {
			newLigand.atoms[i].Position.X += (rand.Float64() - 0.5) * 0.1
			newLigand.atoms[i].Position.Y += (rand.Float64() - 0.5) * 0.1
			newLigand.atoms[i].Position.Z += (rand.Float64() - 0.5) * 0.1
		}
		if isCollisionFree(newLigand, minDistance) {
			return newLigand
		}

	}
}

func RotateLigand(ligand Molecule, maxAngle float64) Molecule {
	newLigand := ligand

	// Randomly choose rotation axis
	axis := Position3d{
		X: rand.Float64()*2 - 1,
		Y: rand.Float64()*2 - 1,
		Z: rand.Float64()*2 - 1,
	}
	axis.Normalize()

	// Randomly choose rotation angle
	theta := (rand.Float64()*2 - 1) * maxAngle

	// Apply rotation to each atom
	for i := range newLigand.atoms {
		newLigand.atoms[i].Position = RotateAtom(newLigand.atoms[i].Position, axis, theta)
	}

	return newLigand
}

func RotateAtom(pos, axis Position3d, theta float64) Position3d {
	cosTheta := math.Cos(theta)
	sinTheta := math.Sin(theta)
	dot := pos.Dot(axis)

	// Rodrigues' rotation formula
	x := pos.X*cosTheta + sinTheta*(axis.Y*pos.Z-axis.Z*pos.Y) + axis.X*dot*(1-cosTheta)
	y := pos.Y*cosTheta + sinTheta*(axis.Z*pos.X-axis.X*pos.Z) + axis.Y*dot*(1-cosTheta)
	z := pos.Z*cosTheta + sinTheta*(axis.X*pos.Y-axis.Y*pos.X) + axis.Z*dot*(1-cosTheta)

	return Position3d{X: x, Y: y, Z: z}
}

// isCollisionFree checks for atomic collisions in the ligand based on a minimum allowed distance.
// Input: a Molecule ligand, a float64 minimum distance
// Output: a bool indicating if the ligand is collision-free
func isCollisionFree(ligand Molecule, minDistance float64) bool {
	for i := 0; i < len(ligand.atoms); i++ {
		for j := i + 1; j < len(ligand.atoms); j++ {
			dx := ligand.atoms[i].Position.X - ligand.atoms[j].Position.X
			dy := ligand.atoms[i].Position.Y - ligand.atoms[j].Position.Y
			dz := ligand.atoms[i].Position.Z - ligand.atoms[j].Position.Z
			distance := math.Sqrt(dx*dx + dy*dy + dz*dz)
			if distance < minDistance {
				return false
			}
		}
	}
	return true
}

// CalculateEnergy computes the electrostatic potential energy between the protein and ligand atoms using Coulomb's law.
// Input: a Molecule protein, a Molecule ligand
// Output: a float64 energy value
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

// Normalize scales the vector to have a magnitude of 1
func (v *Position3d) Normalize() {
	mag := v.Magnitude()
	if mag > 0 {
		v.X /= mag
		v.Y /= mag
		v.Z /= mag
	}
}

// Magnitude calculates the vector's magnitude
func (v Position3d) Magnitude() float64 {
	return math.Sqrt(v.X*v.X + v.Y*v.Y + v.Z*v.Z)
}

// Dot calculates the dot product between two vectors
func (v Position3d) Dot(other Position3d) float64 {
	return v.X*other.X + v.Y*other.Y + v.Z*other.Z
}
