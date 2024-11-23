package main

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io/fs"
	"os"
	"path/filepath"
	"strings"
)

func main() {
	outputDir := "./PDB_splitted"

	pdbs := ReadDirectory("./PDB_origional")

	for _, pdb := range pdbs {
		inputFile := "./PDB_origional/" + pdb.Name()

		WriteProtein(inputFile, outputDir)
		WriteLigand(inputFile, outputDir)
	}

	// for test
	// inputFile := "./PDBdata/1a4h.pdb.gz"

	// WriteProtein(inputFile, outputDir)
	// WriteLigand(inputFile, outputDir)
}

func ReadDirectory(dir string) []fs.DirEntry {
	//read in all files in the given directory
	files, err := os.ReadDir(dir)
	if err != nil {
		panic(err)
	}

	return files
}

func WriteProtein(inputFile, outputDir string) {
	// Open the gzipped PDB file
	gzipFile, err := os.Open(inputFile)
	if err != nil {
		fmt.Println("Error opening input file:", err)
		return
	}
	defer gzipFile.Close()

	// Create a gzip reader
	reader, err := gzip.NewReader(gzipFile)
	if err != nil {
		fmt.Println("Error creating gzip reader:", err)
		return
	}
	defer reader.Close()

	// Get the base name of the input file (e.g., "input.pdb.gz" -> "input")
	inputBase := strings.TrimSuffix(strings.TrimSuffix(filepath.Base(inputFile), ".gz"), ".pdb")

	// Construct the output file path in the given directory
	outputFilePath := filepath.Join(outputDir, inputBase+"_protein.pdb")

	// Open the output file for writing
	outFile, err := os.Create(outputFilePath)
	if err != nil {
		fmt.Println("Error creating output file:", err)
		return
	}
	defer outFile.Close()

	// Create a buffered writer for the output file
	writer := bufio.NewWriter(outFile)

	// Flags to track whether we're in Model 1 or in a single-model PDB file
	inModel1 := false
	multiModel := false // Detect if the file contains multiple models

	// Process the input file line by line
	scanner := bufio.NewScanner(reader)
	for scanner.Scan() {
		line := scanner.Text()

		// Detect MODEL start (multi-model case)
		if strings.HasPrefix(line, "MODEL") {
			if strings.TrimSpace(line) == "MODEL        1" {
				inModel1 = true
			} else if multiModel && inModel1 {
				// Stop processing after ENDMDL of Model 1 in multi-model PDB
				break
			}
			multiModel = true
			continue
		}

		// Stop processing after ENDMDL (multi-model case)
		if strings.HasPrefix(line, "ENDMDL") && inModel1 {
			break
		}

		// If no MODEL keyword is present, assume a single-model PDB
		if !multiModel {
			inModel1 = true
		}

		// Write ATOM records only if in Model 1 or single-model case
		if inModel1 && strings.HasPrefix(line, "ATOM") {
			_, err := writer.WriteString(line + "\n")
			if err != nil {
				fmt.Println("Failed to write to output file:", err)
				return
			}
		}
	}

	// Flush the writer to ensure all data is written
	err = writer.Flush()
	if err != nil {
		fmt.Println("Error flushing output file:", err)
		return
	}

	// Check for any errors during scanning
	if err := scanner.Err(); err != nil {
		fmt.Println("Error reading input file:", err)
		return
	}

	fmt.Println("Filtered ATOM records written to:", outputFilePath)
}

func WriteLigand(inputFile, outputDir string) {
	// Open the gzipped PDB file
	gzipFile, err := os.Open(inputFile)
	if err != nil {
		fmt.Println("Error opening input file:", err)
		return
	}
	defer gzipFile.Close()

	// Create a gzip reader
	reader, err := gzip.NewReader(gzipFile)
	if err != nil {
		fmt.Println("Error creating gzip reader:", err)
		return
	}
	defer reader.Close()

	// Get the base name of the input file (e.g., "input.pdb.gz" -> "input")
	inputBase := strings.TrimSuffix(strings.TrimSuffix(filepath.Base(inputFile), ".gz"), ".pdb")

	// Construct the output file path in the given directory
	outputFilePath := filepath.Join(outputDir, inputBase+"_ligand.pdb")

	// Open the output file for writing
	outFile, err := os.Create(outputFilePath)
	if err != nil {
		fmt.Println("Error creating output file:", err)
		return
	}
	defer outFile.Close()

	// Create a buffered writer for the output file
	writer := bufio.NewWriter(outFile)

	// Flags to track whether we're in Model 1 or in a single-model PDB file
	inModel1 := false
	multiModel := false // Detect if the file contains multiple models

	// Read the gzipped PDB file line by line
	scanner := bufio.NewScanner(reader)
	for scanner.Scan() {
		line := scanner.Text()

		// Detect MODEL start (multi-model case)
		if strings.HasPrefix(line, "MODEL") {
			if strings.TrimSpace(line) == "MODEL        1" {
				inModel1 = true
			} else if multiModel && inModel1 {
				// Stop processing after ENDMDL of Model 1 in multi-model PDB
				break
			}
			multiModel = true
			continue
		}

		// Stop processing after ENDMDL (multi-model case)
		if strings.HasPrefix(line, "ENDMDL") && inModel1 {
			break
		}

		// If no MODEL keyword is present, assume a single-model PDB
		if !multiModel {
			inModel1 = true
		}

		// Check if the line starts with "HETATM" and exclude water molecules (HOH)
		if inModel1 && strings.HasPrefix(line, "HETATM") {
			resName := strings.TrimSpace(line[17:20]) // Residue name is in columns 18-20
			if resName != "HOH" {                     // Exclude water molecules
				_, err := writer.WriteString(line + "\n")
				if err != nil {
					fmt.Println("Failed to write to output file:", err)
					return
				}
			}
		}
	}

	// Flush the writer to ensure all data is written
	err = writer.Flush()
	if err != nil {
		fmt.Println("Error flushing output file:", err)
		return
	}

	// Check for any errors during scanning
	if err := scanner.Err(); err != nil {
		fmt.Println("Error reading input file:", err)
		return
	}

	fmt.Println("Filtered HETATM records from Model 1 written to:", outputFilePath)
}
