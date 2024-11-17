package main

import (
	"fmt"
	"runtime"
)

func main() {
	ligand, err := ParseMol2("Data/Ligands/ligand_charged.mol2")
	Check(err)
	fmt.Println(len(ligand.atoms))
	protein, err2 := ParseMol2("Data/Proteins/protein-charge.mol2")
	Check(err2)
	fmt.Println(len(protein.atoms))
	iterations := 3000
	temperature := 310.15
	numProcs := runtime.NumCPU()
	fmt.Println("Old Energy", CalculateEnergy(protein, ligand))
	//newLigand := SimulateEnergyMinimization(protein, ligand, iterations, temperature)
	newLigand2 := SimulateEnergyMinimizationParallel(protein, ligand, iterations, temperature, numProcs)
	//fmt.Println("New Energy", CalculateEnergy(protein, newLigand))
	fmt.Println("New Energy Parallel", CalculateEnergy(protein, newLigand2))
}

func Check(err error) {
	if err != nil {
		panic(err)
	}
}
