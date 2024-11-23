package main

import (
	"bufio"
	"fmt"
	"os"
	"strings"
)

func main() {
	// Open the input CSV file
	inputFile, err := os.Open("extended_PLAS20K.csv")
	if err != nil {
		fmt.Println("Error opening input file:", err)
		return
	}
	defer inputFile.Close()

	// Create the output text file
	outputFile, err := os.Create("PLAS20K_pdb_ids.txt")
	if err != nil {
		fmt.Println("Error creating output file:", err)
		return
	}
	defer outputFile.Close()

	// Initialize scanners and writers
	scanner := bufio.NewScanner(inputFile)
	writer := bufio.NewWriter(outputFile)

	// Skip the header line
	scanner.Scan() // Reads the first line but does nothing with it

	// Process each subsequent line in the input file
	// Collect PDB IDs
	var pdbIDs []string

	// Process each subsequent line in the input file
	for scanner.Scan() {
		line := scanner.Text()

		// Split the line by commas and get the first column (PDB_ID)
		columns := strings.Split(line, ",")
		if len(columns) > 0 {
			pdbID := columns[0] // First column
			pdbIDs = append(pdbIDs, pdbID)
		}
	}

	// Join all PDB IDs with commas and write to the file
	_, err = writer.WriteString(strings.Join(pdbIDs, ",") + "\n")
	if err != nil {
		fmt.Println("Error writing to output file:", err)
		return
	}

	// Flush the writer to ensure all data is written to the file
	err = writer.Flush()
	if err != nil {
		fmt.Println("Error flushing output file:", err)
		return
	}

	// Check for any errors while reading the input file
	if err := scanner.Err(); err != nil {
		fmt.Println("Error reading input file:", err)
		return
	}

	fmt.Println("PDB IDs have been successfully written to PLAS20K_pdb_ids.txt")
}
