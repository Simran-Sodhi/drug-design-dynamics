package main

import (
	"math"
	"testing"
)

func TestSimulateMultipleLigands(t *testing.T) {
	protein := createMockProtein(1.0, -1.0)
	ligands := createMockLigands(3)
	iterations := 1000
	temperature := 300.0
	numProcs := 4

	minLigands, minEnergies := SimulateMultipleLigands(protein, ligands, iterations, true, temperature, numProcs)

	if len(minLigands) != len(ligands) || len(minEnergies) != len(ligands) {
		t.Fatalf("Expected %d ligands and energies, got %d ligands and %d energies", len(ligands), len(minLigands), len(minEnergies))
	}
}

func TestSimulateMultipleLigandsParallel(t *testing.T) {
	protein := createMockProtein(1.0, -1.0)
	ligands := createMockLigands(5)
	iterations := 2000
	temperature := 300.0
	numProcs := 2

	minLigands, minEnergies := SimulateMultipleLigandsParallel(protein, ligands, iterations, true, temperature, numProcs)

	if len(minLigands) != len(ligands) || len(minEnergies) != len(ligands) {
		t.Fatalf("Expected %d ligands and energies, got %d ligands and %d energies", len(ligands), len(minLigands), len(minEnergies))
	}
}

func TestCalculateEnergyPositive(t *testing.T) {
	protein := createMockProtein(1.0, 1.0)
	ligand := createMockLigandWithCustomCharge(1.0, 1.0)

	energy := CalculateEnergy(protein, ligand)

	if energy <= 0 {
		t.Errorf("Expected positive energy, got %f", energy)
	}
}

func TestCalculateEnergyNegative(t *testing.T) {
	protein := createMockProtein(1.0, 1.0)
	ligand := createMockLigandWithCustomCharge(-1.0, -1.0)

	energy := CalculateEnergy(protein, ligand)

	if energy >= 0 {
		t.Errorf("Expected negative energy, got %f", energy)
	}
}

func TestAcceptMove(t *testing.T) {
	currentEnergy := 10.0
	newEnergy := 5.0
	temperature := 300.0

	accepted := AcceptMove(currentEnergy, newEnergy, temperature)
	if !accepted {
		t.Errorf("Expected move to be accepted")
	}
}

func TestJitterLigand(t *testing.T) {
	ligand := createMockLigand()
	minDistance := 0.5

	jittered := JitterLigand(ligand, minDistance)
	if !IsCollisionFree(jittered, minDistance) {
		t.Errorf("Expected collision-free ligand after jittering")
	}
}

func TestRotateLigand(t *testing.T) {
	ligand := createMockLigand()
	oldLigand := createMockLigand()
	maxAngle := math.Pi / 4

	rotated := RotateLigand(ligand, maxAngle)
	if rotated.atoms[0].Position == oldLigand.atoms[0].Position {
		t.Errorf("Expected ligand to be rotated")
	}
}

func TestShiftLigandCloserByThreshold(t *testing.T) {
	protein := createMockProtein(1.0, -1.0)
	ligand := createMockLigand()
	threshold := 5.0

	shifted := ShiftLigandCloserByThreshold(ligand, protein, threshold)
	distance1 := Distance(shifted.atoms[0].Position, protein.atoms[0].Position)
	distance2 := Distance(shifted.atoms[1].Position, protein.atoms[1].Position)
	if distance1 > threshold && distance2 > threshold {
		t.Errorf("Expected shifted ligand to be within threshold distance")
	}
}

func TestCopyLigand(t *testing.T) {
	ligand := createMockLigandWithCustomCharge(2.3, -8.7)
	copied := CopyLigand(ligand)

	if &ligand.atoms[0] == &copied.atoms[0] {
		t.Errorf("Expected ligand to be deep copied")
	}
	if ligand.atoms[0] != ligand.atoms[0] || ligand.atoms[1] != ligand.atoms[1] {
		t.Errorf("Ligand atom fields don't match")
	}
}

func TestDistance(t *testing.T) {
	pos1 := Position3d{X: 1, Y: 2, Z: 3}
	pos2 := Position3d{X: 4, Y: 6, Z: 8}
	expected := math.Sqrt(50)
	result := Distance(pos1, pos2)

	if result != expected {
		t.Errorf("Expected distance %f, got %f", expected, result)
	}
}

func TestRotateAtom(t *testing.T) {
	pos := Position3d{X: 1, Y: 0, Z: 0}
	axis := Position3d{X: 0, Y: 0, Z: 1}
	theta := math.Pi / 2
	expected := Position3d{X: 0, Y: 1, Z: 0}
	result := RotateAtom(pos, axis, theta)

	if !almostEqual(result.X, expected.X, 1e-6) ||
		!almostEqual(result.Y, expected.Y, 1e-6) ||
		!almostEqual(result.Z, expected.Z, 1e-6) {
		t.Errorf("Expected rotated position %v, got %v", expected, result)
	}
}

func TestFindClosestAtomDistance(t *testing.T) {
	protein := createMockProtein(2.0, 4.0)
	ligand := createMockLigand()
	expected := math.Sqrt(12.0)
	distance, _, _ := FindClosestAtomDistance(ligand, protein)
	if distance != expected {
		t.Errorf("Expected distance %f, got %f", expected, distance)
	}
}

// createMockProtein creates a mock protein molecule with two atoms at specified charges and fixed positions.
// Input: two float64 charges charge1 and charge2
// Output: a Molecule
func createMockProtein(charge1, charge2 float64) Molecule {
	return Molecule{
		atoms: []Atom{
			{Position: Position3d{X: 0, Y: 0, Z: 0}, Charge: charge1},
			{Position: Position3d{X: 1, Y: 1, Z: 1}, Charge: charge2},
		},
	}
}

// createMockLigands generates a specified number of mock ligands with default atomic properties.
// Input: an int count
// Output: a slice of Molecule
func createMockLigands(count int) []Molecule {
	ligands := make([]Molecule, count)
	for i := 0; i < count; i++ {
		ligands[i] = createMockLigand()
	}
	return ligands
}

// createMockLigand creates a mock ligand molecule with two atoms of fixed charges and positions
// Input: none
// Output: a Molecule
func createMockLigand() Molecule {
	return Molecule{
		atoms: []Atom{
			{Position: Position3d{X: 4, Y: 4, Z: 4}, Charge: 1},
			{Position: Position3d{X: 3, Y: 3, Z: 3}, Charge: -1},
		},
	}
}

// createMockLigandWithCustomCharge creates a mock ligand molecule with two atoms having specified charges and fixed positions.
// Input: two float64 charges charge1 and charge2
// Output: a Molecule
func createMockLigandWithCustomCharge(charge1, charge2 float64) Molecule {
	return Molecule{
		atoms: []Atom{
			{Position: Position3d{X: 4, Y: 4, Z: 4}, Charge: charge1},
			{Position: Position3d{X: 3, Y: 3, Z: 3}, Charge: charge2},
		},
	}
}

// almostEqual compares two float64 values and checks if they are within a specified tolerance epsilon.
// Input: two float64 values a and b, a float64 epsilon
// Output: a bool indicating if the values are approximately equal
func almostEqual(a, b float64, epsilon float64) bool {
	return math.Abs(a-b) <= epsilon
}
