//
//  main.cpp
//  Block_Code_15_11
//
//  Created by Stefan Węgrzyn on 20/10/2019.
//  Copyright © 2019 Stefan Węgrzyn. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <fstream>

#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#define PI 3.141592654

using namespace std;

#define n 15    // codeword length
#define k 11    // information sequence length

#define MIN_INPUT_VECTOR_LEN    100000   // 10^(5)
#define MAX_INPUT_VECTOR_LEN    10000000 // 10^(7)

const int G[k][n] = { {1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                      {1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                      {0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                      {1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
                      {1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0},
                      {0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
                      {1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0},
                      {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0},
                      {1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0},
                      {0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0},
                      {1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1} };

const int H_T[n][n-k] = { {1, 0, 0, 0},
                          {0, 1, 0, 0},
                          {1, 1, 0, 0},
                          {0, 0, 1, 0},
                          {1, 0, 1, 0},
                          {0, 1, 1, 0},
                          {1, 1, 1, 0},
                          {0, 0, 0, 1},
                          {1, 0, 0, 1},
                          {0, 1, 0, 1},
                          {1, 1, 0, 1},
                          {0, 0, 1, 1},
                          {1, 0, 1, 1},
                          {0, 1, 1, 1},
                          {1, 1, 1, 1} };

const string inputFilePath = "input.txt";
const string outputFilePath = "output.txt";

//************************************************************************

float gauss(float mean, float sigma);
void kanal(float es_n0, long dl_kan, int *wej, float *wyj);

void getCodeword(int* infoSeq, int* codeWord);
void getSyndrome(int* codeWord, int* syndrome);

//************************************************************************

int main(int argc, const char * argv[])
{
    srand((unsigned int)time(NULL));
    
    ifstream file_reader(inputFilePath, ios::in);
    if ( !file_reader.is_open() ) {
        cout << "Could not open file: " << inputFilePath << endl;
        return -1;
    }
    ofstream file_writer(outputFilePath, ios::trunc);
    if ( !file_writer.is_open() ) {
        cout << "Could not open file: " << outputFilePath << endl;
        return -1;
    }

    // Acquire input parameters from the user
    
    long inputLen = 0;
    long extInputLen = 0;
    
    float min_Eb_N0 = -5;  // -5
    float max_Eb_N0 = 10;  // 6-11
    float step_Eb_N0 = 0.25; // 0.5
    
    cout << "Determine number of sent bits (min " << MIN_INPUT_VECTOR_LEN;
    cout << " max " << MAX_INPUT_VECTOR_LEN << "): ";
    cin >> inputLen;
    inputLen = inputLen - (inputLen % k);   // make inputLen divisible by k
    extInputLen = (inputLen / k) * n;
    
    cout << "Determine MINIMAL value of 'Eb/N0'[dB]: ";
    cin >> min_Eb_N0;
    
    cout << "Determine MAXIMAL value of 'Eb/N0'[dB]: ";
    cin >> max_Eb_N0;
    
    cout << "Determine change step of 'Eb/N0'[dB]: ";
    cin >> step_Eb_N0;
    
    // Read input vector
    char bit[4];
    int* inputVector = new int[inputLen];
    for (long idx = 0; idx < inputLen; idx++)
    {
        file_reader.getline(bit,4);     // "bit + \t\n\r"
        inputVector[idx] = atoi(bit);
    }
    
    for (float step = min_Eb_N0; step <= max_Eb_N0; step += step_Eb_N0)
    {
        double BER_no_correction = 0, BER_with_correction = 0;
        
        // NO CORRECTION
        
        // no additional control bits added to the input vector
        
        float* noisyInput = new float[inputLen];
        kanal(step, inputLen, inputVector, noisyInput); // add noise to input bit stream

        for (long idx = 0; idx < inputLen; idx++)
        {
            if ( (noisyInput[idx] >= 0) && (inputVector[idx] == 0) ) BER_no_correction++;
            else if  ( (noisyInput[idx] < 0) && (inputVector[idx] == 1) ) BER_no_correction++;
        }
        
        // WITH CORRECTION
        
        // adding control bits to every information sequence in the input vector
        
        int* extInputVector = new int[extInputLen];
        for (long seqIdx = 0; seqIdx < (inputLen/k); seqIdx++)
        {
            // determine codeword for given information sequence
            getCodeword(&inputVector[seqIdx*k], &extInputVector[seqIdx*n]);
        }
        
        delete[] noisyInput;
        noisyInput = new float[extInputLen];
        kanal(step, extInputLen, extInputVector, noisyInput); // add noise to input bit stream
        
        // decode noisy input
        int* outputVector = new int[extInputLen];
        for (long idx = 0; idx < extInputLen; idx++)
        {
            if (noisyInput[idx] >= 0) outputVector[idx] = 1;
            else outputVector[idx] = 0;
        }
        
        // apply correction algorithm
        for (long wordIdx = 0; wordIdx < (extInputLen/n); wordIdx++)
        {
            // determine syndrome for given codeword
            int* syndrome = new int[n-k];
            getSyndrome(&outputVector[wordIdx*n], syndrome);
            
            int errPos = 1*syndrome[0] + 2*syndrome[1] + 4*syndrome[2] + 8*syndrome[3];
            if ( (errPos != 0) && (errPos <= n) )
            {
                outputVector[wordIdx*n + errPos - 1] = abs(outputVector[wordIdx*n + errPos - 1] - 1);
            }
            
            delete[] syndrome;
        }
        
        // extract information sequence from codeword
        int* outputInfo = new int[inputLen];
        for (long seqIdx = 0; seqIdx < (inputLen/k); seqIdx++)
        {
            // factual information sequence is stored at
            // positions 2,4,5,6,8,9,10,11,12,13,14 in codeword
            outputInfo[seqIdx*k + 0]  = outputVector[seqIdx*n + 2];
            outputInfo[seqIdx*k + 1]  = outputVector[seqIdx*n + 4];
            outputInfo[seqIdx*k + 2]  = outputVector[seqIdx*n + 5];
            outputInfo[seqIdx*k + 3]  = outputVector[seqIdx*n + 6];
            outputInfo[seqIdx*k + 4]  = outputVector[seqIdx*n + 8];
            outputInfo[seqIdx*k + 5]  = outputVector[seqIdx*n + 9];
            outputInfo[seqIdx*k + 6]  = outputVector[seqIdx*n + 10];
            outputInfo[seqIdx*k + 7]  = outputVector[seqIdx*n + 11];
            outputInfo[seqIdx*k + 8]  = outputVector[seqIdx*n + 12];
            outputInfo[seqIdx*k + 9]  = outputVector[seqIdx*n + 13];
            outputInfo[seqIdx*k + 10] = outputVector[seqIdx*n + 14];
        }
        
        // determine BER
        for (long idx = 0; idx < inputLen; idx++)
        {
            // comparing only information sequence
            if ( outputInfo[idx] != inputVector[idx] ) BER_with_correction++;
        }
        
        file_writer << "Eb/N0= ";
        file_writer << fixed << setprecision(2) << setw(5) << setfill(' ') << step;
        file_writer << " BER_no_correction= ";
        file_writer << fixed << setprecision(10) << setw(12) << setfill(' ') << BER_no_correction/inputLen;
        file_writer << " BER_with_correction= ";
        file_writer << fixed << setprecision(10) << setw(12) << setfill(' ') << BER_with_correction/inputLen;
        file_writer << endl;
        
        delete[] outputInfo;
        delete[] outputVector;
        delete[] extInputVector;
        delete[] noisyInput;
    }
    
    
    delete[] inputVector;
    file_reader.close();
    file_writer.close();
    
    return 0;
}




//******************************************************************

// Function kanal changes binary values into bipolar ones (-1/+1) and adds noise
// *wej - Input vector of binary values (0/1)
// *wyj - Output vector of real numbers
// es_n0 - Es/N0
// dl_kan - the number of input bits
void kanal(float es_n0, long dl_kan, int *wej, float *wyj)
{
    float mean=0;
    float es=1;
    float sygnal;
    float sigma;
    float s_n;
    long y;
    
    s_n=(float) pow(10, (es_n0/10));
    sigma=(float) sqrt (es/(2*s_n));
    
    for (y=0; y<dl_kan; y++)
    {
        sygnal = 2 * *(wej+y)-1; // change the binary value (0/1) into symbol (-1/+1)
        *(wyj+y)=sygnal+gauss(mean,sigma);  // noise addition
    }
}

//*******************************************************************

float gauss(float mean, float sigma)
{
    double x;
    double z;
    
    z=(double)rand()/RAND_MAX;
    if (z==1.0) z=0.9999999;
    x=sigma*sqrt(2.0*log(1.0/(1.0-z)));
    
    z=(double)rand()/RAND_MAX;
    if (z==1.0) z=0.9999999;
    return((float)(mean+x*cos(2*PI*z)));
}

//*******************************************************************

void getCodeword(int* infoSeq, int* codeWord)
{
    // U = m * G
    for (int i = 0; i < n; i++)
    {
        codeWord[i] = 0;
        for (int j = 0; j < k; j++)
        {
            codeWord[i] += infoSeq[j] * G[j][i];
        }
        codeWord[i] = codeWord[i] % 2;  // modulo 2
    }
}

//*******************************************************************

void getSyndrome(int* codeWord, int* syndrome)
{
    // s = r * H_T
    for (int i = 0; i < n-k; i++)
    {
        syndrome[i] = 0;
        for (int j = 0; j < n; j++)
        {
            syndrome[i] += codeWord[j] * H_T[j][i];
        }
        syndrome[i] = syndrome[i] % 2;  // modulo 2
    }
}

//*******************************************************************
