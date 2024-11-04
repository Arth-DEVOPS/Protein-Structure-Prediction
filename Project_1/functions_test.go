package main

import (
	"testing"
)

// TestPredictHelix tests the predictHelix function with various amino acid sequences
func TestPredictHelix(t *testing.T) {
	tests := []struct {
		sequence string
		expected bool
	}{
		{"AAAAAA", true},  // Example of a sequence likely to form a helix
		{"PPPPPP", false}, // Example of a sequence unlikely to form a helix
		{"AKLVRA", true},  // Mixed sequence with sufficient helix-forming residues
	}

	for _, tt := range tests {
		result := PredictHelix(tt.sequence)
		if result != tt.expected {
			t.Errorf("predictHelix(%q) = %v; want %v", tt.sequence, result, tt.expected)
		}
	}
}

// TestPredictSheet tests the predictSheet function with various amino acid sequences
func TestPredictSheet(t *testing.T) {
	tests := []struct {
		sequence string
		expected bool
	}{
		{"VVVVV", true},  // Example of a sequence likely to form a sheet
		{"DDDDD", false}, // Example of a sequence unlikely to form a sheet
		{"ILVYV", true},  // Mixed sequence with sufficient sheet-forming residues
	}

	for _, tt := range tests {
		result := PredictBetaSheet(tt.sequence)
		if result != tt.expected {
			t.Errorf("predictSheet(%q) = %v; want %v", tt.sequence, result, tt.expected)
		}
	}
}

// TestPredictTurn tests the predictTurn function with various amino acid sequences
func TestPredictTurn(t *testing.T) {
	tests := []struct {
		sequence string
		expected bool
	}{
		{"GGGG", true},  // Example of a sequence likely to form a turn
		{"LLLL", false}, // Example of a sequence unlikely to form a turn
		{"NGST", true},  // Mixed sequence with turn-forming residues
	}

	for _, tt := range tests {
		result := PredictTurn(tt.sequence)
		if result != tt.expected {
			t.Errorf("predictTurn(%q) = %v; want %v", tt.sequence, result, tt.expected)
		}
	}
}
