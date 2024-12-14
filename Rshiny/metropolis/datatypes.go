package main

import "math"

// universal constant
const K = 9e9

// these constants are not universal but adjustable parameters for the simulation
const THRESHOLD = 5        // max distance between protein and ligand in Å to start the simulation with
const MINDISTANCE = 0.5    // min Distance to maintain between atoms of the ligand for it to be considered a valid perturbation
const MAXANGLE = 0.79      //max angle in radians (= ~45 degrees) for which atom rotation is allowed when perturbing ligands
const TEMPERATURE = 310.15 // body temperature

type Molecule struct {
	atoms []Atom
}

type Atom struct {
	Position Position3d //Coordinates
	Charge   float64    // Charge
}

type Position3d struct {
	X, Y, Z float64
}

type MultipleLigandSimulationOutput struct {
	Ligand []Molecule
	Energy []float64
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

// Add performs vector addition
func (v Position3d) Add(other Position3d) Position3d {
	return Position3d{X: v.X + other.X, Y: v.Y + other.Y, Z: v.Z + other.Z}
}

// Scale multiplies the vector by a scalar value
func (v Position3d) Scale(scalar float64) Position3d {
	return Position3d{X: v.X * scalar, Y: v.Y * scalar, Z: v.Z * scalar}
}
