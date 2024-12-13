package main

import (
	"encoding/csv"
	"fmt"
	"math"
	"math/rand"
	"os"
	"sort"
	"strconv"
	"strings"
	"time"
)

type ProteinData struct {
	PDB_ID       string
	Experimental float64
	DELTA_TOTAL  float64
	VDWAALS      float64
	EEL          float64
	EPB          float64
	ENPOLAR      float64
}

type DecisionTree struct {
	value     float64
	feature   int
	threshold float64
	left      *DecisionTree
	right     *DecisionTree
}

type RandomForest struct {
	trees  []*DecisionTree
	nTrees int
}

func engineerFeatures(data []ProteinData) [][]float64 {
	features := make([][]float64, len(data))
	for i, d := range data {
		// Create interaction terms and polynomial features
		vdwaalsEel := d.VDWAALS * d.EEL
		epbEnpolar := d.EPB * d.ENPOLAR
		deltaSquared := d.DELTA_TOTAL * d.DELTA_TOTAL
		totalEnergy := d.DELTA_TOTAL + d.VDWAALS + d.EEL
		solvationEnergy := d.EPB + d.ENPOLAR

		features[i] = []float64{
			d.DELTA_TOTAL,
			d.VDWAALS,
			d.EEL,
			d.EPB,
			d.ENPOLAR,
			vdwaalsEel,
			epbEnpolar,
			deltaSquared,
			totalEnergy,
			solvationEnergy,
			math.Abs(d.DELTA_TOTAL),
			math.Abs(d.EEL),
		}
	}
	return features
}

func findBestSplit(X [][]float64, y []float64, features []int) (int, float64, float64) {
	bestScore := math.Inf(1)
	bestFeature := 0
	bestThreshold := 0.0
	bestValue := 0.0

	for _, feature := range features {
		values := make([]float64, len(X))
		for i := range X {
			values[i] = X[i][feature]
		}
		sort.Float64s(values)

		for i := 0; i < len(values)-1; i++ {
			threshold := (values[i] + values[i+1]) / 2
			leftSum, rightSum := 0.0, 0.0
			leftCount, rightCount := 0, 0
			leftSquaredSum, rightSquaredSum := 0.0, 0.0

			for j := range X {
				if X[j][feature] <= threshold {
					leftSum += y[j]
					leftSquaredSum += y[j] * y[j]
					leftCount++
				} else {
					rightSum += y[j]
					rightSquaredSum += y[j] * y[j]
					rightCount++
				}
			}

			if leftCount > 0 && rightCount > 0 {
				leftMean := leftSum / float64(leftCount)
				rightMean := rightSum / float64(rightCount)
				leftVar := leftSquaredSum/float64(leftCount) - leftMean*leftMean
				rightVar := rightSquaredSum/float64(rightCount) - rightMean*rightMean
				score := float64(leftCount)*leftVar + float64(rightCount)*rightVar

				if score < bestScore {
					bestScore = score
					bestFeature = feature
					bestThreshold = threshold
					bestValue = (leftSum + rightSum) / float64(leftCount+rightCount)
				}
			}
		}
	}
	return bestFeature, bestThreshold, bestValue
}

func (dt *DecisionTree) train(X [][]float64, y []float64, depth int) {
	if len(X) == 0 || depth >= 25 {
		dt.value = 0
		return
	}

	sum := 0.0
	count := 0
	for _, val := range y {
		if !math.IsNaN(val) && !math.IsInf(val, 0) {
			sum += val
			count++
		}
	}

	if count > 0 {
		dt.value = sum / float64(count)
	}

	if len(X) < 5 || depth >= 25 {
		return
	}

	nFeatures := len(X[0])
	featureCount := int(math.Sqrt(float64(nFeatures))) + 1
	features := rand.Perm(nFeatures)[:featureCount]

	feature, threshold, value := findBestSplit(X, y, features)
	dt.feature = feature
	dt.threshold = threshold
	dt.value = value

	leftX := make([][]float64, 0)
	leftY := make([]float64, 0)
	rightX := make([][]float64, 0)
	rightY := make([]float64, 0)

	for i := range X {
		if X[i][feature] <= threshold {
			leftX = append(leftX, X[i])
			leftY = append(leftY, y[i])
		} else {
			rightX = append(rightX, X[i])
			rightY = append(rightY, y[i])
		}
	}

	if len(leftX) > 0 && len(rightX) > 0 {
		dt.left = &DecisionTree{}
		dt.right = &DecisionTree{}
		dt.left.train(leftX, leftY, depth+1)
		dt.right.train(rightX, rightY, depth+1)
	}
}

func (dt *DecisionTree) predict(x []float64) float64 {
	if dt.left == nil || dt.right == nil {
		return dt.value
	}
	if x[dt.feature] <= dt.threshold {
		return dt.left.predict(x)
	}
	return dt.right.predict(x)
}

func NewRandomForest(nTrees int) *RandomForest {
	return &RandomForest{
		trees:  make([]*DecisionTree, nTrees),
		nTrees: nTrees,
	}
}

func (rf *RandomForest) train(X [][]float64, y []float64) {
	fmt.Println("Training Random Forest...")
	nFolds := 5
	foldSize := len(X) / nFolds

	for i := 0; i < rf.nTrees; i++ {
		foldIndex := i % nFolds
		startIdx := foldIndex * foldSize
		endIdx := startIdx + foldSize

		validX := X[startIdx:endIdx]
		validY := y[startIdx:endIdx]

		trainX := make([][]float64, 0, len(X)-foldSize)
		trainY := make([]float64, 0, len(y)-foldSize)

		for j := 0; j < len(X); j++ {
			if j < startIdx || j >= endIdx {
				trainX = append(trainX, X[j])
				trainY = append(trainY, y[j])
			}
		}

		n := len(trainX)
		sampleSize := int(float64(n) * 0.8)
		sampleX := make([][]float64, sampleSize)
		sampleY := make([]float64, sampleSize)

		for j := 0; j < sampleSize; j++ {
			idx := rand.Intn(n)
			sampleX[j] = trainX[idx]
			sampleY[j] = trainY[idx]
		}

		tree := &DecisionTree{}
		tree.train(sampleX, sampleY, 0)
		rf.trees[i] = tree

		if (i+1)%50 == 0 {
			validError := 0.0
			validCount := 0
			for j := range validX {
				pred := tree.predict(validX[j])
				if !math.IsNaN(pred) && !math.IsInf(pred, 0) {
					diff := math.Abs(pred - validY[j])
					validError += diff
					validCount++
				}
			}
			if validCount > 0 {
				validError /= float64(validCount)
				fmt.Printf("Trained %d trees, Validation MAE: %.4f\n", i+1, validError)
			}
		}
	}
}

func (rf *RandomForest) predict(x []float64) float64 {
	if len(rf.trees) == 0 {
		return 0
	}

	predictions := make([]float64, len(rf.trees))
	for i, tree := range rf.trees {
		predictions[i] = tree.predict(x)
	}

	sort.Float64s(predictions)
	q1 := predictions[len(predictions)/4]
	q3 := predictions[3*len(predictions)/4]
	iqr := q3 - q1
	lowerBound := q1 - 1.5*iqr
	upperBound := q3 + 1.5*iqr

	sum := 0.0
	count := 0
	for _, pred := range predictions {
		if pred >= lowerBound && pred <= upperBound {
			sum += pred
			count++
		}
	}

	if count > 0 {
		return sum / float64(count)
	}
	return predictions[len(predictions)/2]
}

// Main function and other utility functions remain the same as in the original code

func normalizeFeatures(X [][]float64) [][]float64 {
	if len(X) == 0 {
		return X
	}

	nFeatures := len(X[0])
	means := make([]float64, nFeatures)
	stdDevs := make([]float64, nFeatures)

	for j := 0; j < nFeatures; j++ {
		sum := 0.0
		for i := range X {
			sum += X[i][j]
		}
		means[j] = sum / float64(len(X))
	}

	for j := 0; j < nFeatures; j++ {
		sumSquared := 0.0
		for i := range X {
			diff := X[i][j] - means[j]
			sumSquared += diff * diff
		}
		stdDevs[j] = math.Sqrt(sumSquared / float64(len(X)))
		if stdDevs[j] == 0 {
			stdDevs[j] = 1
		}
	}

	normalized := make([][]float64, len(X))
	for i := range X {
		normalized[i] = make([]float64, nFeatures)
		for j := range X[i] {
			normalized[i][j] = (X[i][j] - means[j]) / stdDevs[j]
		}
	}

	return normalized
}

func readCSV(filename string) ([]ProteinData, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, fmt.Errorf("error opening file: %v", err)
	}
	defer file.Close()

	reader := csv.NewReader(file)
	reader.Comma = ','

	records, err := reader.ReadAll()
	if err != nil {
		return nil, fmt.Errorf("error reading CSV: %v", err)
	}

	var data []ProteinData
	for i, record := range records {
		if i == 0 || len(record) != 7 || strings.Contains(record[0], "Note") {
			continue
		}

		pdbID := strings.TrimSpace(record[0])
		if pdbID == "" {
			continue
		}

		// Parse experimental value first
		experimental := 0.0
		if strings.TrimSpace(record[1]) != "" {
			if exp, err := strconv.ParseFloat(strings.TrimSpace(record[1]), 64); err == nil {
				experimental = exp
			}
		}

		if experimental == 0 {
			continue
		}

		deltaTotal, err1 := strconv.ParseFloat(strings.TrimSpace(record[2]), 64)
		vdwaals, err2 := strconv.ParseFloat(strings.TrimSpace(record[3]), 64)
		eel, err3 := strconv.ParseFloat(strings.TrimSpace(record[4]), 64)
		epb, err4 := strconv.ParseFloat(strings.TrimSpace(record[5]), 64)
		enpolar, err5 := strconv.ParseFloat(strings.TrimSpace(record[6]), 64)

		if err1 != nil || err2 != nil || err3 != nil || err4 != nil || err5 != nil {
			continue
		}

		proteinData := ProteinData{
			PDB_ID:       pdbID,
			Experimental: experimental,
			DELTA_TOTAL:  deltaTotal,
			VDWAALS:      vdwaals,
			EEL:          eel,
			EPB:          epb,
			ENPOLAR:      enpolar,
		}
		data = append(data, proteinData)
	}

	if len(data) == 0 {
		return nil, fmt.Errorf("no valid data records found in file")
	}

	fmt.Printf("Successfully parsed %d valid records\n", len(data))
	return data, nil
}

func main() {
	rand.Seed(time.Now().UnixNano())

	fmt.Println("Starting protein binding prediction model...")

	data, err := readCSV("data.csv")
	if err != nil {
		fmt.Printf("Error reading CSV: %v\n", err)
		return
	}

	X := make([][]float64, len(data))
	y := make([]float64, len(data))

	for i, d := range data {
		X[i] = []float64{
			d.DELTA_TOTAL,
			d.VDWAALS,
			d.EEL,
			d.EPB,
			d.ENPOLAR,
		}
		y[i] = d.Experimental
	}

	fmt.Println("Normalizing features...")
	X = normalizeFeatures(X)

	trainSize := int(float64(len(data)) * 0.8)
	indices := rand.Perm(len(data))

	trainX := make([][]float64, trainSize)
	trainY := make([]float64, trainSize)
	testX := make([][]float64, len(data)-trainSize)
	testY := make([]float64, len(data)-trainSize)

	for i := 0; i < trainSize; i++ {
		trainX[i] = X[indices[i]]
		trainY[i] = y[indices[i]]
	}
	for i := trainSize; i < len(data); i++ {
		testX[i-trainSize] = X[indices[i]]
		testY[i-trainSize] = y[indices[i]]
	}

	rf := NewRandomForest(500) // Increased number of trees
	rf.train(trainX, trainY)

	fmt.Println("\nMaking predictions...")
	predictions := make([]float64, len(testX))
	for i, x := range testX {
		predictions[i] = rf.predict(x)
	}

	var mse, totalError float64
	errorBuckets := make(map[string]int)

	fmt.Println("\nTest Set Results:")
	fmt.Println("==================")
	for i := range predictions {
		diff := math.Abs(predictions[i] - testY[i])
		mse += diff * diff
		totalError += diff

		fmt.Printf("PDB_ID: %s, Predicted: %.3f, Actual: %.3f, Diff: %.3f\n",
			data[indices[i+trainSize]].PDB_ID,
			predictions[i],
			testY[i],
			diff)

		switch {
		case diff < 0.3:
			errorBuckets["Good (< 0.3)"]++
		case diff < 0.6:
			errorBuckets["Moderate (0.3-0.6)"]++
		case diff < 0.9:
			errorBuckets["Fair (0.6-0.9)"]++
		default:
			errorBuckets["Poor (> 0.9)"]++
		}
	}

	mse /= float64(len(predictions))
	rmse := math.Sqrt(mse)
	mae := totalError / float64(len(predictions))

	fmt.Printf("\nModel Performance:\n")
	fmt.Printf("Mean Squared Error: %.4f\n", mse)
	fmt.Printf("Root Mean Squared Error: %.4f\n", rmse)
	fmt.Printf("Mean Absolute Error: %.4f\n", mae)

	fmt.Printf("\nError Distribution:\n")
	for category, count := range errorBuckets {
		percentage := float64(count) * 100 / float64(len(predictions))
		fmt.Printf("%s: %d predictions (%.1f%%)\n", category, count, percentage)
	}

	fmt.Printf("\nFeature Importance Analysis:\n")
	features := []string{"DELTA_TOTAL", "VDWAALS", "EEL", "EPB", "ENPOLAR"}
	for i, feature := range features {
		var importance float64
		for _, x := range trainX {
			importance += math.Abs(x[i])
		}
		importance /= float64(len(trainX))
		fmt.Printf("%s: %.4f\n", feature, importance)
	}

	// In your main function, after calculating predictions:
	accuracy, precision := calculateAccuracyMetrics(predictions, testY, 0.5)
	fmt.Printf("Accuracy: %.4f\n", accuracy)
	fmt.Printf("Precision: %.4f\n", precision)

}

// Add these functions to your existing code

func calculateAccuracyMetrics(predictions []float64, actual []float64, threshold float64) (float64, float64) {
	if len(predictions) != len(actual) {
		return 0, 0
	}

	var truePositives, falsePositives, trueNegatives, falseNegatives float64

	for i := range predictions {
		diff := math.Abs(predictions[i] - actual[i])

		if diff <= threshold {
			if predictions[i] >= actual[i] {
				truePositives++
			} else {
				trueNegatives++
			}
		} else {
			if predictions[i] >= actual[i] {
				falsePositives++
			} else {
				falseNegatives++
			}
		}
	}

	accuracy := (truePositives + trueNegatives) / float64(len(predictions))
	precision := truePositives / (truePositives + falsePositives)

	return accuracy, precision
}
