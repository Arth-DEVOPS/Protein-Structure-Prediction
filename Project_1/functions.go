package main

import "strings"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// CHAU FASMAN MODEL FOR SECONDARY STRUCTURE PREDICTION /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Initialize the propensities for each amino acid
var propensities = map[rune]AminoAcidPropensities{
	'A': {1.42, 0.83, 0.66}, // Alanine- strong α former; β indifferent
	'R': {0.98, 0.93, 0.95}, // Arginine- α indifferent; β indifferent
	'N': {0.67, 0.89, 1.56}, // Asparagine- α breaker; β indifferent
	'D': {1.01, 0.54, 1.46}, // Aspartic Acid- weak α former; strong β breaker
	'C': {0.70, 1.19, 1.19}, // Cysteine- α indifferent; β former
	'Q': {1.11, 1.10, 0.98}, // Glutamine- α former; β former
	'E': {1.51, 0.37, 0.74}, // Glutamic Acid- strong α former; strong β breaker
	'G': {0.57, 0.75, 1.64}, // Glycine- strong α breaker; β breaker
	'H': {1.00, 0.87, 0.95}, // Histidine- weak α former; β indifferent
	'I': {1.08, 1.60, 0.47}, // Isoleucine- α former; strong β former
	'L': {1.21, 1.30, 0.59}, // Leucine- strong α former; β former
	'K': {1.16, 0.74, 1.01}, // Lysine- α former; β breaker
	'M': {1.45, 1.05, 0.60}, // Methionine- strong α former; β former
	'F': {1.13, 1.38, 0.60}, // Phenylalanine- α former; β former
	'P': {0.57, 0.55, 1.52}, // Proline- strong α breaker; strong β breaker
	'S': {0.77, 0.75, 1.43}, // Serine- α indifferent; β breaker
	'T': {0.83, 1.19, 0.96}, // Threonine- α indifferent; β former
	'W': {1.08, 1.37, 0.96}, // Tryptophan- α former; β former
	'Y': {0.69, 1.47, 1.14}, // Tyrosine- α breaker; strong β former
	'V': {1.06, 1.70, 0.50}, // Valine- α former; strong β former
}

// Initialise bend propensities for each amino acid
var bendProbabilitiesTable = map[rune]BendProbabilities{
	'A': {0.060, 0.076, 0.035, 0.058},
	'R': {0.070, 0.106, 0.099, 0.085},
	'N': {0.161, 0.083, 0.191, 0.091},
	'D': {0.147, 0.110, 0.179, 0.081},
	'C': {0.149, 0.050, 0.117, 0.128},
	'Q': {0.074, 0.098, 0.037, 0.098},
	'E': {0.056, 0.060, 0.077, 0.064},
	'G': {0.102, 0.085, 0.190, 0.152},
	'H': {0.140, 0.047, 0.093, 0.054},
	'I': {0.043, 0.034, 0.013, 0.056},
	'L': {0.061, 0.025, 0.036, 0.070},
	'K': {0.055, 0.115, 0.072, 0.095},
	'M': {0.068, 0.082, 0.014, 0.055},
	'F': {0.059, 0.041, 0.065, 0.065},
	'P': {0.102, 0.301, 0.034, 0.068},
	'S': {0.120, 0.139, 0.125, 0.106},
	'T': {0.086, 0.108, 0.065, 0.079},
	'W': {0.089, 0.073, 0.064, 0.167},
	'Y': {0.082, 0.065, 0.114, 0.125},
	'V': {0.062, 0.048, 0.028, 0.053},
}

// PredictHelix()
// Input: a String sequence
// Output: a boolean after checking if a sequence segment has a high propensity for alpha helix formation
func PredictHelix(sequence string) bool {
	helixNucleationCount := 0.0
	propensityScore := 0.0
	AACount := 0.0
	breakers := 0.0
	formers := 0.0

	for _, aminoAcid := range sequence {
		if propensities[aminoAcid].alphaHelix >= 1.05 {
			propensityScore += propensities[aminoAcid].alphaHelix
			helixNucleationCount++ //increase counter when encountered with amino acid with tendency to form alpha helix
			formers++
		} else if propensities[aminoAcid].alphaHelix >= 1.0 {
			propensityScore += propensities[aminoAcid].alphaHelix
			helixNucleationCount += 0.5 // weak helix former counts as half
			formers++
		} else if propensities[aminoAcid].alphaHelix <= 0.69 {
			propensityScore += propensities[aminoAcid].alphaHelix
			breakers++
		}
		AACount++
	}
	// Calculate average alpha helix nucleation count
	average := propensityScore / AACount

	// Check if average exceeds threshold of 1.03
	if average >= 1.03 {
		// Check if the number of breakers exceeds one-third of the amino acid count
		if breakers > 2 {
			return false
		}
		// Check if number of formers is less than 3 or if the nucleation count is less than 4
		if helixNucleationCount < 4 || formers < 3 {
			return false
		}
		return true
	}

	return false
}

// PredictBetaSheet()
// Input: a String sequence
// Output: a boolean after checking if a sequence segment has a high propensity for beta sheet formation
func PredictBetaSheet(sequence string) bool {
	propensityScore := 0.0
	formers := 0.0
	breakers := 0.0
	AACount := 0.0

	for _, aminoAcid := range sequence {
		if propensities[aminoAcid].betaSheet >= 1 {
			propensityScore += propensities[aminoAcid].betaSheet
			formers++ // Increase counter when encountered with amino acid with tendency to form beta sheet
		} else if propensities[aminoAcid].betaSheet <= 0.75 {
			propensityScore += propensities[aminoAcid].betaSheet
			breakers++ // Add to the beta sheet breaker count
		}
		AACount++
	}

	// Calculate average beta sheet nucleation count
	average := propensityScore / AACount

	// Check if average exceeds threshold of 1.05
	if average >= 1.05 {
		// Check if the number of breakers is more than 1
		if breakers > 1 {
			return false
		}
		return true
	}

	return false
}

// PredictTurn()
// Input: a String sequence
// Output: a boolean after checking if a sequence segment has a high propensity for forming turns/bends
func PredictTurn(sequence string) bool {
	// Calculate the average bend probability for the four-residue segment
	posProbability := 1.0 //initilaise positional probability
	turnPropensity := 0.0 //initialise AA propensity to turn

	for i, aminoAcid := range sequence[:4] {
		bendProb := bendProbabilitiesTable[aminoAcid]
		turnPropensity += propensities[aminoAcid].turn

		switch i {
		case 0:
			posProbability *= bendProb.Position0
		case 1:
			posProbability *= bendProb.Position1
		case 2:
			posProbability *= bendProb.Position2
		case 3:
			posProbability *= bendProb.Position3
		}
	}
	averageTurnPropensity := turnPropensity / 4.0
	// Check if product of the positional probabilities and averate turn propensity exceeds the threshold
	if posProbability < 0.000075 && averageTurnPropensity <= 1.0 {
		return false
	}

	// Ensure the two middle residues meet a minimum bend probability (e.g., 0.5)
	if propensities[rune(sequence[1])].turn < 0.5 || propensities[rune(sequence[2])].turn < 0.5 {
		return false
	}

	return true
}

// ClassifyOverlap()
// Input: a String sequence
// Output: a string corresponding to whether an overlapping region favors helix or sheet structure
func ClassifyOverlap(sequence string) string {
	helixSum := 0.0
	sheetSum := 0.0
	turnSum := 0.0

	// Loop through each amino acid in the sequence to sum probabilities
	for _, aminoAcid := range sequence {
		helixSum += propensities[aminoAcid].alphaHelix
		sheetSum += propensities[aminoAcid].betaSheet
		turnSum += propensities[aminoAcid].turn // Assuming you have a turn property in the propensities
	}

	// Determine classification based on summed probabilities
	if turnSum > helixSum && turnSum > sheetSum {
		return "Turn"
	} else if helixSum > sheetSum {
		return "Helix"
	} else if sheetSum > helixSum {
		return "Beta Sheet"
	} else {
		return "Undefined"
	}
}

// PredictStructure()
// Input: a String sequence
// Output: a string that predicts the secondary structure of a protein by employing the Chou-Fasman model
func PredictStructure(sequence string) string {
	structure := make([]rune, len(sequence))

	// Loop needs to account for different sliding windows for different structure predictions
	for i := 0; i < len(sequence); i++ {
		if i <= len(sequence)-6 && PredictHelix(sequence[i:i+6]) {
			for j := i; j < i+6 && j < len(sequence); j++ {
				structure[j] = 'H' // Helix prediction
			}
		}
		if i <= len(sequence)-5 && PredictBetaSheet(sequence[i:i+5]) {
			for j := i; j < i+5 && j < len(sequence); j++ {
				structure[j] = 'E' // Beta-sheet prediction
			}
		}
		if i <= len(sequence)-4 && PredictTurn(sequence[i:i+4]) {
			for j := i; j < i+4 && j < len(sequence); j++ {
				structure[j] = 'T' // Turn prediction
			}
		}
		if i <= len(sequence)-4 {
			overlapClass := ClassifyOverlap(sequence[i : i+4])
			for j := i; j < i+4 && j < len(sequence); j++ {
				if overlapClass == "Helix" {
					structure[j] = 'H'
				} else if overlapClass == "Beta Sheet" {
					structure[j] = 'E'
				} else if overlapClass == "Turn" {
					structure[j] = 'T'
				}
			}
		}
	}
	// Replace any unclassified positions with 'C' for coil
	for i := range structure {
		if structure[i] == 0 {
			structure[i] = 'C' // 'C' stands for Coil
		}
	}

	return string(structure)
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// GOR MODEL FOR SECONDARY STRUCTURE PREDICTION /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// GORPredictStructure implements a simplified version of the GOR model
func GORPredictStructure(sequence string) string {
	sequence = strings.ToUpper(sequence)
	structure := make([]rune, len(sequence))

	// Predict structure for each amino acid based on its neighbors
	for i := range sequence {
		helixScore := GORScore(sequence, i, GORModel.helixMatrix)
		sheetScore := GORScore(sequence, i, GORModel.sheetMatrix)
		turnScore := GORScore(sequence, i, GORModel.turnMatrix)

		// Classify based on highest score
		if helixScore > sheetScore && helixScore > turnScore {
			structure[i] = 'H' // Helix
		} else if sheetScore > helixScore && sheetScore > turnScore {
			structure[i] = 'E' // Sheet
		} else {
			structure[i] = 'T' // Turn
		}
	}
	return string(structure)
}

// GORScore calculates the score of a specific position considering its neighbors
func GORScore(sequence string, index int, matrix map[rune]float64) float64 {
	windowSize := 2 // Consider 2 residues on each side
	score := 0.0

	// Loop through the window around the current position
	for offset := -windowSize; offset <= windowSize; offset++ {
		pos := index + offset
		if pos >= 0 && pos < len(sequence) {
			score += matrix[rune(sequence[pos])]
		}
	}

	return score
}
