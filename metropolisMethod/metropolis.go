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
		ligands[i] = ShiftLigandCloserByThreshold(ligand, protein, THRESHOLD)
	}
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
	for i := 0; i < iterations; i++ {
		var newLigand Molecule
		if rotate {
			newLigand = JitterAndRotateLigand(currentLigand, MINDISTANCE, MAXANGLE)
		} else {
			newLigand = JitterLigand(currentLigand, MINDISTANCE)
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
	for i := 0; i < iterations; i++ {
		var newLigand Molecule
		prevLigand := CopyLigand(currentLigand)
		if rotate {
			newLigand = JitterAndRotateLigand(currentLigand, MINDISTANCE, MAXANGLE)
		} else {
			newLigand = JitterLigand(currentLigand, MINDISTANCE)
		}
		newEnergy := CalculateEnergy(protein, newLigand)
		if AcceptMove(currentEnergy, newEnergy, temperature) {
			currentLigand = newLigand
			currentEnergy = newEnergy
		} else {
			currentLigand = prevLigand
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
		if IsCollisionFree(newLigand, minDistance) {
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
		if IsCollisionFree(newLigand, minDistance) {
			return newLigand
		}
	}
}

// RotateLigand rotates the ligand around a random axis by a random angle within the specified maxAngle.
// Input: a Molecule ligand, a float64 maxAngle
// Output: a rotated Molecule ligand
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

// RotateAtom rotates a 3D position around a specified axis by an angle using Rodrigues' rotation formula.
// Input: a Position3d pos, a Position3d axis, a float64 theta
// Output: a rotated Position3d
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

// IsCollisionFree checks for atomic collisions in the ligand based on a minimum allowed distance.
// Input: a Molecule ligand, a float64 minimum distance
// Output: a bool indicating if the ligand is collision-free
func IsCollisionFree(ligand Molecule, minDistance float64) bool {
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
			// to ensure non-zero
			if distance < 1e-6 {
				distance = 1e-6
			}
			energy += K * (atomP.Charge * atomL.Charge) / distance //kQ1,Q2/d^2
		}
	}
	return energy
}

// ShiftLigandCloserByThreshold shifts the ligand closer to the protein if its closest atom distance exceeds the threshold.
// Input: a Molecule ligand, a Molecule protein, a float64 threshold
// Output: a shifted Molecule ligand
func ShiftLigandCloserByThreshold(ligand, protein Molecule, threshold float64) Molecule {
	closestDistance, ligandAtom, proteinAtom := FindClosestAtomDistance(ligand, protein)
	if closestDistance > threshold {
		// Calculate shift vector
		shiftVector := Position3d{
			X: proteinAtom.X - ligandAtom.X,
			Y: proteinAtom.Y - ligandAtom.Y,
			Z: proteinAtom.Z - ligandAtom.Z,
		}
		// Normalize shift vector and scale to close the gap
		shiftVector.Normalize()
		shiftVector = shiftVector.Scale(closestDistance - threshold)
		// Apply the shift to the entire ligand
		for i := range ligand.atoms {
			ligand.atoms[i].Position = ligand.atoms[i].Position.Add(shiftVector)
		}
	}
	return ligand
}

// FindClosestAtomDistance finds the minimum distance and the corresponding atom positions between the ligand and the protein.
// Input: a Molecule ligand, a Molecule protein
// Output: a float64 minimum distance and corresponding closest atom positions
func FindClosestAtomDistance(ligand, protein Molecule) (float64, Position3d, Position3d) {
	minDistance := math.MaxFloat64
	var ligandAtomPos, proteinAtomPos Position3d
	for _, ligAtom := range ligand.atoms {
		for _, protAtom := range protein.atoms {
			dist := Distance(ligAtom.Position, protAtom.Position)
			if dist < minDistance {
				minDistance = dist
				ligandAtomPos = ligAtom.Position
				proteinAtomPos = protAtom.Position
			}
		}
	}
	return minDistance, ligandAtomPos, proteinAtomPos
}

// Distance calculates the Euclidean distance between two 3D vectors.
// Input: two Position3d vectors a and b
// Output: a float64 distance
func Distance(a, b Position3d) float64 {
	dx := a.X - b.X
	dy := a.Y - b.Y
	dz := a.Z - b.Z
	return math.Sqrt(dx*dx + dy*dy + dz*dz)
}

// CopyLigand creates a deep copy of the given ligand to avoid reference sharing.
// Input: a Molecule ligand
// Output: a deep copy of the Molecule ligand
func CopyLigand(ligand Molecule) Molecule {
	newAtoms := make([]Atom, len(ligand.atoms))
	for i, atom := range ligand.atoms {
		newAtoms[i] = Atom{
			Position: Position3d{
				X: atom.Position.X,
				Y: atom.Position.Y,
				Z: atom.Position.Z,
			},
			Charge: atom.Charge,
		}
	}
	return Molecule{
		atoms: newAtoms,
	}
}
