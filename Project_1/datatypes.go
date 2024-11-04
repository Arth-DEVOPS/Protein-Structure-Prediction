package main

type AminoAcidPropensities struct {
	alphaHelix float64
	betaSheet  float64
	turn       float64
}

type BendProbabilities struct {
	Position0 float64 // Probability at position 0
	Position1 float64 // Probability at position 1
	Position2 float64 // Probability at position 2
	Position3 float64 // Probability at position 3
}

// GORPropensities represents context-based propensities for secondary structures
type GORPropensities struct {
	helixMatrix map[rune]float64
	sheetMatrix map[rune]float64
	turnMatrix  map[rune]float64
}

// GORModelPredictor stores the probability matrices used in the GOR method
var GORModel = GORPropensities{
	helixMatrix: map[rune]float64{
		// Simplified example values (actual values would be based on a statistical model)
		'A': 1.2, 'C': 0.5, 'D': 0.8, 'E': 1.3, 'F': 0.9,
		// Add more amino acids here...
	},
	sheetMatrix: map[rune]float64{
		'A': 0.8, 'C': 1.3, 'D': 0.7, 'E': 0.6, 'F': 1.4,
		// Add more amino acids here...
	},
	turnMatrix: map[rune]float64{
		'A': 0.9, 'C': 0.7, 'D': 1.4, 'E': 0.8, 'F': 0.5,
		// Add more amino acids here...
	},
}
