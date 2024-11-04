package main

import (
	"fmt"
	"os"
)

func main() {
	// Check if a sequence was provided as a command-line argument
	if len(os.Args) < 2 {
		fmt.Println("Please provide a sequence as a command-line argument.")
		return
	}

	// Get the sequence from the command-line arguments
	sequence := os.Args[1]

	// Print the prediction for the provided sequence
	fmt.Println("Chou-Fasman Prediction:")
	fmt.Printf("Sequence: %s\n", sequence)
	fmt.Printf("Structure: %s\n", PredictStructure(sequence))

	//fmt.Println("\nGOR Prediction:")
	//fmt.Printf("Sequence: %s\n", sequence)
	//fmt.Printf("Structure: %s\n", GORPredictStructure(sequence))
}

//MVLSPADKTNVKAAW
