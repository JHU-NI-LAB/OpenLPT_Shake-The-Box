/*
 * MergeSort.h
 *
 *  Created on: May 13, 2018
 *      Author: sut210
 */

#ifndef MERGESORT_H_
#define MERGESORT_H_


void CopyArray(int* A, int iBegin, int iEnd, int* B) {
	for (int k = iBegin; k < iEnd; ++k) {
		B[k] = A[k];
	}
}

void Merge(int* A, int iBegin, int iMiddle, int iEnd, int* B) {
	int i = iBegin, j = iMiddle;
	for(int k = iBegin; k < iEnd; ++k) {
		if (i < iMiddle && (j >= iEnd || A[i] >= A[j])) {
			B[k] = A[i];
			i = i + 1;
		} else {
			B[k] = A[j];
			j = j + 1;
		}
	}
}

void SplitMerge(int* B, int iBegin, int iEnd, int* A) {
	if (iEnd - iBegin < 2) return;
	int iMiddle = (iEnd + iBegin) / 2;
	SplitMerge(A, iBegin, iMiddle, B);
	SplitMerge(A, iMiddle, iEnd, B);
	Merge(B, iBegin, iMiddle, iEnd, A);
}

void MergeSort(int* A, int* B, int n) {
	CopyArray(A, 0, n, B);
	SplitMerge(B, 0, n, A);
}


#endif /* MERGESORT_H_ */
