#include "my_solver.h"

#define MATDIM 1024

double *pMatA, *pMatB, *pMatC, *pMatTrs;

void loopUnrolling() {
  initMat();

  std::chrono::high_resolution_clock::time_point start, stop;
  std::chrono::milliseconds duration;

  // Time measurement test
  CHECK_TIME_START;
  std::this_thread::sleep_for(std::chrono::seconds(1));
  CHECK_TIME_END;
  OUTPUT_DURATION("Time measurement test");

  // Original Matrix multiplication
  CHECK_TIME_START;
  multiplySquareMatrices1(pMatC, pMatA, pMatB, MATDIM);
  CHECK_TIME_END;
  OUTPUT_DURATION("multiplySquareMatrices1");

  // Transpose
  for (int i = 0; i < MATDIM * MATDIM; ++i)
    pMatC[i] = 0;
  CHECK_TIME_START;
  transposeMatrix(pMatB, MATDIM);
  multiplySquareMatrices2(pMatC, pMatA, pMatTrs, MATDIM);
  CHECK_TIME_END;
  OUTPUT_DURATION("multiplySquareMatrices2");

  // Loop unrolling increases the program's speed by eliminating loop control instruction and loop test instructions
  // Loop unrolling - 8 times
  for (int i = 0; i < MATDIM * MATDIM; ++i)
    pMatC[i] = 0;
  CHECK_TIME_START;
  transposeMatrix(pMatB, MATDIM);
  multiplySquareMatrices3(pMatC, pMatA, pMatTrs, MATDIM);
  CHECK_TIME_END;
  OUTPUT_DURATION("multiplySquareMatrices3");

  // Reduce performance:
  // Cause an increase in instruction cache misses
  // Increased usage of register in a single iteration to store temporary variables
  // Loop unrolling - 32 times
  for (int i = 0; i < MATDIM * MATDIM; ++i)
    pMatC[i] = 0;
  CHECK_TIME_START;
  transposeMatrix(pMatB, MATDIM);
  multiplySquareMatrices4(pMatC, pMatA, pMatTrs, MATDIM);
  CHECK_TIME_END;
  OUTPUT_DURATION("multiplySquareMatrices4");

  delete[] pMatA;
  delete[] pMatB;
  delete[] pMatC;
  delete[] pMatTrs;
}

void initMat() {
  pMatA = new double[MATDIM * MATDIM];
  pMatB = new double[MATDIM * MATDIM];
  pMatC = new double[MATDIM * MATDIM];
  pMatTrs = new double[MATDIM * MATDIM];

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0);

  for (int i = 0; i < MATDIM * MATDIM; ++i) {
    pMatA[i] = (double) dis(gen);
    pMatB[i] = (double) dis(gen);
  }
}

void transposeMatrix(double *pMat, int MatSize) {
  for (int i = 0; i < MatSize; i++)
    for (int j = 0; j < MatSize; j++)
      pMatTrs[i * MatSize + j] = pMat[j * MatSize + i];
}

void multiplySquareMatrices1(double *pDestMatrix, double *pLeftMatrix, double *pRightMatrix, int MatSize) {
  for (int i = 0; i < MatSize; ++i) {
    for (int j = 0; j < MatSize; ++j) {
      pDestMatrix[i * MatSize + j] = 0;
      for (int k = 0; k < MatSize; ++k)
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k] * pRightMatrix[k * MatSize + j];
    }
  }
}

void multiplySquareMatrices2(double *pDestMatrix, double *pLeftMatrix, double *pRightMatrix, int MatSize) {
  for (int i = 0; i < MatSize; ++i) {
    for (int j = 0; j < MatSize; ++j) {
      pDestMatrix[i * MatSize + j] = 0;
      for (int k = 0; k < MatSize; ++k)
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k] * pRightMatrix[j * MatSize + k];
    }
  }
}

void multiplySquareMatrices3(double *pDestMatrix, double *pLeftMatrix, double *pRightMatrix, int MatSize) {
  for (int i = 0; i < MatSize; ++i) {
    for (int j = 0; j < MatSize; ++j) {
      pDestMatrix[i * MatSize + j] = 0;
      for (int k = 0; k < MatSize / 8; ++k) {
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 8 + 0] * pRightMatrix[j * MatSize + k * 8 + 0];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 8 + 1] * pRightMatrix[j * MatSize + k * 8 + 1];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 8 + 2] * pRightMatrix[j * MatSize + k * 8 + 2];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 8 + 3] * pRightMatrix[j * MatSize + k * 8 + 3];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 8 + 4] * pRightMatrix[j * MatSize + k * 8 + 4];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 8 + 5] * pRightMatrix[j * MatSize + k * 8 + 5];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 8 + 6] * pRightMatrix[j * MatSize + k * 8 + 6];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 8 + 7] * pRightMatrix[j * MatSize + k * 8 + 7];
      }
    }
  }
}

void multiplySquareMatrices4(double *pDestMatrix, double *pLeftMatrix, double *pRightMatrix, int MatSize) {
  for (int i = 0; i < MatSize; i++) {
    for (int j = 0; j < MatSize; j++) {
      pDestMatrix[i * MatSize + j] = 0;
      for (int k = 0; k < MatSize / 32; k++) {
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 +  0] * pRightMatrix[j * MatSize + k * 32 +  0];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 +  1] * pRightMatrix[j * MatSize + k * 32 +  1];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 +  2] * pRightMatrix[j * MatSize + k * 32 +  2];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 +  3] * pRightMatrix[j * MatSize + k * 32 +  3];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 +  4] * pRightMatrix[j * MatSize + k * 32 +  4];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 +  5] * pRightMatrix[j * MatSize + k * 32 +  5];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 +  6] * pRightMatrix[j * MatSize + k * 32 +  6];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 +  7] * pRightMatrix[j * MatSize + k * 32 +  7];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 +  8] * pRightMatrix[j * MatSize + k * 32 +  8];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 +  9] * pRightMatrix[j * MatSize + k * 32 +  9];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 + 10] * pRightMatrix[j * MatSize + k * 32 + 10];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 + 11] * pRightMatrix[j * MatSize + k * 32 + 11];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 + 12] * pRightMatrix[j * MatSize + k * 32 + 12];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 + 13] * pRightMatrix[j * MatSize + k * 32 + 13];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 + 14] * pRightMatrix[j * MatSize + k * 32 + 14];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 + 15] * pRightMatrix[j * MatSize + k * 32 + 15];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 + 16] * pRightMatrix[j * MatSize + k * 32 + 16];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 + 17] * pRightMatrix[j * MatSize + k * 32 + 17];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 + 18] * pRightMatrix[j * MatSize + k * 32 + 18];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 + 19] * pRightMatrix[j * MatSize + k * 32 + 19];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 + 20] * pRightMatrix[j * MatSize + k * 32 + 20];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 + 21] * pRightMatrix[j * MatSize + k * 32 + 21];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 + 22] * pRightMatrix[j * MatSize + k * 32 + 22];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 + 23] * pRightMatrix[j * MatSize + k * 32 + 23];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 + 24] * pRightMatrix[j * MatSize + k * 32 + 24];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 + 25] * pRightMatrix[j * MatSize + k * 32 + 25];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 + 26] * pRightMatrix[j * MatSize + k * 32 + 26];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 + 27] * pRightMatrix[j * MatSize + k * 32 + 27];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 + 28] * pRightMatrix[j * MatSize + k * 32 + 28];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 + 29] * pRightMatrix[j * MatSize + k * 32 + 29];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 + 30] * pRightMatrix[j * MatSize + k * 32 + 30];
        pDestMatrix[i * MatSize + j] += pLeftMatrix[i * MatSize + k * 32 + 31] * pRightMatrix[j * MatSize + k * 32 + 31];
      }
    }
  }
}
