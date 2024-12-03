package main

import (
	"fmt"
	"runtime"
	"time"
)

func main() {
	dir := "Data/mol2_files"
	ligandFiles, err := findFilesWithSubstring(dir, "ligand")
	Check(err)
	//ligandFiles := []string{"1a0i_ligand.mol2", "1a4h_ligand.mol2", "1ajp_ligand.mol2", "1gi8_ligand.mol2"}
	// for i := range ligandFiles {

	// }
	//ligandList := [4]string{dir + "/1a0i_ligand.mol2", dir + "/1a4h_ligand.mol2", dir + "/1ajp_ligand.mol2", dir + "/1gi8_ligand.mol2"}
	ligandFiles = ligandFiles[:200]
	ligands := make([]Molecule, len(ligandFiles))
	//ligand, err := ParseMol2("Data/Ligands/ligand_charged.mol2")
	for i := range ligandFiles {
		//dir + "/" +
		ligand, err := ParseMol2(ligandFiles[i])
		Check(err)
		ligands[i] = ligand
	}
	// ligand, err := ParseMol2("Data/Ligands/ligand_charged.mol2")
	// Check(err)
	//fmt.Println(len(ligand.atoms))
	//protein, err2 := ParseMol2("Data/Proteins/protein-charge.mol2")
	proteinFile := "1a0i_protein.mol2"
	protein, err2 := ParseMol2(dir + "/" + proteinFile)
	Check(err2)
	//fmt.Println(len(protein.atoms))
	iterations := 3000
	temperature := 310.15
	numProcs := runtime.NumCPU()
	//fmt.Println("Old Energy", CalculateEnergy(protein, ligand))
	//newLigand := SimulateEnergyMinimization(protein, ligand, iterations, temperature)
	start := time.Now()
	_, energyList := SimulateMultipleLigands(protein, ligands, iterations, temperature, numProcs)
	end := time.Since(start)
	fmt.Println("Time taken: ", end)
	// minIndex := 0
	// minEnergy := 0.0
	// for i, energy := range energyList {
	// 	if energy < minEnergy {
	// 		minIndex = i
	// 		minEnergy = energy
	// 	}
	// }
	// minLigands[minIndex].SaveToMol2("minLigand_" + proteinFile)
	//dir+"/"+
	// err3 := UpdateMol2Coordinates(ligandFiles[minIndex], "minLigand_"+proteinFile, minLigands[minIndex])
	// Check(err3)
	// fmt.Println(energyList)
	// start2 := time.Now()
	// //newLigand2 := SimulateEnergyMinimizationParallel(protein, ligands[0], iterations, temperature, numProcs)
	// newLigand2 := SimulateEnergyMinimization(protein, ligands[0], iterations, temperature)
	// end2 := time.Since(start2)
	// fmt.Println("Time taken: ", end2)
	//fmt.Println("New Energy", CalculateEnergy(protein, newLigand))
	// fmt.Println("Old Energy Parallel", CalculateEnergy(protein, ligands[0]))
	// fmt.Println("New Energy Parallel", CalculateEnergy(protein, newLigand2))
	plotEnergy(ligandFiles, energyList)
}

func Check(err error) {
	if err != nil {
		panic(err)
	}
}
