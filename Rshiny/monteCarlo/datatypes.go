package main

const K = 9e9

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
