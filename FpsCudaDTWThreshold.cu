//#include "cuda_def.cuh"
//#include "cuda_proc.h"
#include<time.h>
#include "cuda_def.cuh"
#include <iostream>
#include <vector>
//#include "cuda_kernels.h"
#include <omp.h>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <cmath>
#include <numeric>
#include <iomanip>
#include <algorithm>
// #include <queue>
// #include "clipp.h"
#include <iterator>
#include "time.h"

using namespace std;
#define DBL_MAX 1.7976931348623158e+308 // max value
#define PI 3.1415926535898
#define MAX_NUM 2147483647
#define MAXLENGTH 1023
#define CPUTHREADS 10

using namespace std;

//////////////////////////////////
//function
//////////////////////////////////

__global__ void cuDTW_2048(float *g_allColData, unsigned int *g_allColLength, float *g_allRowData,
                           unsigned int *g_allRowLength, float *g_odata);

float gpuCalculatemnDynamicTimeWarping_2048(const vector<vector<float>> &siga,
                                            const vector<vector<float>> &sigb,
                                            vector<vector<float>> &gpudistResult);
__global__ void cuDTW_ultimate(float *g_allColData,
                               unsigned int *g_allColLength,
                               float *g_allRowData,
                               unsigned int *g_allRowLength, float *g_odata);

float gpuCalculate1nDynamicTimeWarping(const vector<vector<float>> &siga,
                                       const vector<vector<float>> &sigb,
                                       vector<vector<float>> &gpudistResult);
float cpuCalculate1nDynamicTimeWarping(const vector<vector<float>> &siga,
                                       const vector<vector<float>> &sigb,
                                       vector<vector<float>> &cpudistResult);
float cpuDynamicTimeWarping(const std::vector<float> &seq1,
                            const std::vector<float> &seq2);

// using namespace std;
vector<float> GetSingleSignalData(int &SignalId, const string &argv_1);

void FPS(vector<int> &indexList, vector<int > &final_indexList, const string &dataFolder, vector<vector<float>> &TimeSeriesData, vector<vector<float>> &sample_res, vector<vector<float>> &dist_matrix, float threshold);

std::string getCmdRes(const std::string sc);

// const string* parseArgv(int argc, char **argv);

void getTimeSeriesData(const string &dataFolder, vector<vector<float>> &TimeSeriesData);

void delete_ele_in_vector(vector<vector<float> > &init_vector, vector<vector<float> > &delete_Vector);

template <class T>
void convertFromString(T &value, string &s);

void check(vector<vector<float>> sample_res, vector<vector<float>> check_res, float threshold);

void WriteLog();

void delete_ele_inVectorIntvesion(vector<int > &init_vector, vector<int > &delete_Vector);

int main(int argc, char **argv)
{
  // cout << argc << endl;
  // cout << argv[0] << argv[1] << argv[2] << endl;
  
  if(argc != 9){cout << "usage: ./FpsCudaDTWThreshold -i FolderPath -l SampleLength -t Threshold_factor -o Output_file" << endl;
                cout << "-i Folder_path: 'FolderPath' is the folder path of time series, the file in 'FolderPath' must be timeSeries_0.txt, timeSeries_1.txt,...,timeSeries_n.txt. Here, n-1 is the number of time series." << endl;
                cout << "-l Sample_length: 'Sample_length' is the length of sample(for automated threshold determination)." << endl;
                cout << "-o Output_file: 'Output_file' is used to store the selected time series[format: txt]." << endl;
                cout << "-t Threshold_factor: 'Threshold_factor' is used to control the threshold of the FPS algorithm, the recommended value is 0~30." << endl;
                cout << "Program description: Farthest point sampling algorithm based on CUDA DTW distance and our designed threshold for selecting specific time series." << endl;
                exit(0);}

  else{

    // string true_argv[3] = string* parseArgv(argv);
    // int sample_num = atoi(true_argv[1]);
    clock_t start, finish;
    double  duration;

    printf("program start time: ");
    WriteLog();
    cout << endl;

    cout << "###Parameter Information###" << endl;
    cout << "Folder Path: " << argv[2] << endl;
    cout << "Sample Length: " << argv[4] << endl;
    cout << "Output File: " << argv[6] << endl;
    cout << "Threshold factor: " << argv[8] << endl;
    cout << "###Parameter Information###" << endl;
    cout << endl;

    start = clock();
    int sample_length = atoi(argv[4]);
    float Threshold_factor = atof(argv[8]);
    float threshold = 5.4159*sample_length - 27.886 - Threshold_factor;
    cout << "###Threshold Information###" << endl;
    cout << "Threshold: " << threshold << endl;
    cout << "###Threshold Information###" << endl;
    cout << endl;

    vector<vector<float>> sample_res{};
    vector<vector<float>> TimeSeriesData{};
    // const string dataFolder = true_argv[0];
    const string dataFolder = argv[2];
    vector<vector<float> > dist_matrix{};
    vector<int > indexList;
    vector<int > final_indexList;

    FPS(indexList, final_indexList, dataFolder, TimeSeriesData, sample_res, dist_matrix, threshold);
    finish = clock();

    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    printf( "###The run time of sampling is: %f seconds###\n", duration);
    cout << endl;

    // ofstream disfile(true_argv[2], ios::out);
    ofstream disfile(argv[6], ios::out);
    disfile << "###Selected Time Series###" << endl;
    string temp;
    for (int i = 0; i < sample_res.size(); i++)
    {
      for (int j = 0; j < sample_res[i].size(); j++)
      {
        temp = std::to_string(sample_res[i][j]);
        disfile << temp;
        disfile << " ";
      }
      disfile << endl;
    }

    disfile << "###Selected Index###" << endl;

    for(int i = 0; i < final_indexList.size(); i++)
    {
      // cout << final_indexList[i];
      temp = std::to_string(final_indexList[i]);
      disfile << temp;
      disfile << " ";
    }

    printf("program end time: ");
    WriteLog();
    cout << endl;
    cout << "#######################" << endl;
    cout << "Final Pick Quantity: " << final_indexList.size() << endl;
    cout << "#######################" << endl;

    disfile.close();

    // check program
    vector<vector<float>> check_res{};
    check(sample_res,check_res, threshold);

    return 0;
  }
  
}

//////////////////////////////////
//function
//////////////////////////////////

void WriteLog()
{
	time_t timeReal;
	time(&timeReal);
	timeReal = timeReal + 8*3600;
	tm* t = gmtime(&timeReal); 
  // cout << timeReal << endl;
	printf("%d-%02d-%02d %02d:%02d:%02d\n", t->tm_year + 1900, t->tm_mon + 1, t->tm_mday, t->tm_hour, t->tm_min, t->tm_sec); 
}


void getTimeSeriesData(const string &dataFolder, vector<vector<float>> &TimeSeriesData)
{
    // input a folder that load all time series, output a vector that load them.
    // Note: the format of time series file in folder must be '%s_%d.txt', '%s' indicates the 
    // prefix name of file and '%d' indicates the index(consecutive integers) of file.
    // the file format example: timeSeries_0.txt, timeSeries_1.txt, timeSeries_2.txt,...,timeSeries_999.txt
    // The above example indicates folder has 1000 time series.

    // string shell_command = "ls " + dataFolder + "/timeSeries_*" + " | wc -l"; // get shell command.
    string shell_command = "find " + dataFolder + " -name " + "timeSeries_* " + "| wc -l";
    const int file_num = atoi((getCmdRes(shell_command)).c_str()); //get the file number.
  
    // the following loop for getting all time series data.
    for(int i = 0; i < file_num; i++)
    {
      vector<float> single_data = GetSingleSignalData(i,dataFolder);
      TimeSeriesData.push_back(single_data);
    }
    
}

void delete_ele_in_vector(vector<vector<float> > &init_vector, vector<vector<float> > &delete_Vector)
{
  
  for(int i = 0; i < delete_Vector.size(); i++)
  {
    init_vector.erase(remove(init_vector.begin(),init_vector.end(),delete_Vector[i]),init_vector.end());
  }
}

void delete_ele_inVectorIntvesion(vector<int > &init_vector, vector<int > &delete_Vector)
{
 
  for(int i = 0; i < delete_Vector.size(); i++)
  {
    init_vector.erase(remove(init_vector.begin(),init_vector.end(),delete_Vector[i]),init_vector.end());
  }
}

void check(vector<vector<float>> sample_res, vector<vector<float>> check_res, float threshold)
{
  gpuCalculatemnDynamicTimeWarping_2048(sample_res,sample_res,check_res);
  ofstream disfile("test.txt", ios::out);
  string temp;
  for (int i = 0; i < check_res.size(); i++)
  {
    for (int j = 0; j < check_res[i].size(); j++)
    {
      temp = std::to_string(check_res[i][j]);
      disfile << temp;
      disfile << " ";
      if(check_res[i][j] < threshold && i != j)
      { 
        cout << check_res[i][j] << endl;
        cout << "Program Error!" << endl;
        exit(0);
      }
    }
    disfile << endl;
  }
}


void FPS(vector<int> &indexList, vector<int > &final_indexList, const string &dataFolder, vector<vector<float>> &TimeSeriesData, vector<vector<float>> &sample_res, vector<vector<float>> &dist_matrix, float threshold)
{
  // input: sample number and time series data. output: final sample result.

  // Load data from dataFolder into TimeSeriesData
  getTimeSeriesData(dataFolder, TimeSeriesData);
  cout << "###Total Data Volume: " << TimeSeriesData.size() << "###" << endl;
  cout << endl;
  indexList.clear();
  
  // initialization
  for(int i = 0; i < TimeSeriesData.size(); i++)
  {
    indexList.push_back(i);
  }
 
  // select a random number
  int random_num = rand()%TimeSeriesData.size(); 
  
  sample_res.push_back(TimeSeriesData[random_num]);
  final_indexList.push_back(indexList[random_num]);
  delete_ele_in_vector(TimeSeriesData, sample_res);
  cout << "###Current Data Volume: " << TimeSeriesData.size() << "###" << endl;
  cout << "###Updating Sample List times: " << 1 << "###" << endl;
  cout << endl; 

  vector<int > tempIndexList{random_num};
  delete_ele_inVectorIntvesion(indexList, tempIndexList);

  const int loopNum = TimeSeriesData.size() - 1;

  cout << "######Start Iteration######" << endl;
  for(int k = 0; k < loopNum; k++)
  {
    cout << "###Current Data Volume: " << TimeSeriesData.size() << "###" << endl;
    // cout << endl;

    if(TimeSeriesData.size() == 0)
    {
      break;
    }

    vector<vector<float>> tempTimeseriesDataList{sample_res[sample_res.size()-1]};
    gpuCalculatemnDynamicTimeWarping_2048(tempTimeseriesDataList, TimeSeriesData, dist_matrix);

    float maxValue = *max_element(dist_matrix[0].begin(), dist_matrix[0].end());

    float minValue = *min_element(dist_matrix[0].begin(), dist_matrix[0].end());
    cout << "###minDistValue: "  << minValue << ", " << "maxDistValue: " << maxValue << "###" << endl;
    // cout << endl;

    if(maxValue < threshold)
    {
      break;
    }

    else
    {

      int maxPosition = max_element(dist_matrix[0].begin(), dist_matrix[0].end()) - dist_matrix[0].begin();
      cout << "###Updating Sample List times: " << k+2 << "###" << endl;
      cout << endl;
      sample_res.push_back(TimeSeriesData[maxPosition]);
      final_indexList.push_back(indexList[maxPosition]);
      
      vector<vector<float>> tempDeleteEle{};
      vector<int> tempIndexList2{};

      vector<int > delete_list1{};
      // vector<int > delete_list2{};

      for(int i = 0; i < dist_matrix[0].size(); i++)
      {
        if(dist_matrix[0][i] < threshold)
        {
          // tempDeleteEle.push_back(TimeSeriesData[i]);
          delete_list1.push_back(i);

          // tempIndexList2.push_back(indexList[i]);
          // delete_list2.push_back(i);
        }
      }

      // vector<vector<float>> temp_dele_data1; 
      // temp_dele_data1.push_back(TimeSeriesData[maxPosition]);
      vector<vector<float>> newTimeseriesData{};
      vector<int > newIndexList{};
      delete_list1.push_back(maxPosition);

      // cout << "Updating time series data list: " << k << endl;
      // delete_ele_in_vector(TimeSeriesData, tempDeleteEle);
      // delete_ele_in_vector(TimeSeriesData, temp_dele_data1);
      for(int i = 0; i < TimeSeriesData.size(); i++)
      {
        // cout << "delete iteration: " << i << endl;
        if (std::find(delete_list1.begin(), delete_list1.end(), i) == delete_list1.end())
        {
          newTimeseriesData.push_back(TimeSeriesData[i]);
          newIndexList.push_back(indexList[i]);
        }
      }

      TimeSeriesData.swap(newTimeseriesData);
      indexList.swap(newIndexList);

      // vector<int > temp_dele_data2; 
      // temp_dele_data2.push_back(indexList[maxPosition]);

      // delete_ele_inVectorIntvesion(indexList, tempIndexList2);
      // delete_ele_inVectorIntvesion(indexList, temp_dele_data2);

      dist_matrix.clear();

    }
  }
}

std::string getCmdRes(const std::string sc)

{ 
  /*
  a function that construct a shell command in the program and return its result.
  */
    FILE* crs = popen(sc.c_str(), "r"); // execute the shell command. 

    char result[1024] = "0";

    fread(result, sizeof(char), sizeof(result), crs);
    
    if (NULL != crs)
    {

    fclose(crs);
    crs = NULL;
    }

    std::string res = result;

    return res;

}

template <class T>
void convertFromString(T &value, string &s)
{
  std::stringstream ss(s);
  ss >> value;
}

void ZScoreNormalize(vector<float> &signals)
{
  float sum = accumulate(signals.begin(), signals.end(), 0.0);
  float mean = sum / signals.size();

  float acc = 0.0;
  for (size_t i = signals.size(); i--;)
  {
    signals[i] = signals[i] - mean;
    acc += signals[i] * signals[i];
  }

  float deviation = sqrt(acc / signals.size());

  for (size_t i = signals.size(); i--;)
  {
    signals[i] /= deviation;
  }
}

vector<float> GetSingleSignalData(int &SignalId, const string &argv_1)
{
  
  ifstream SignalFile;

  string TempId = std::to_string(SignalId);
 
  string TempString = argv_1 + "/" + "timeSeries_" + TempId + ".txt";
  const char *FileName = TempString.data();
  SignalFile.open(FileName, ios::in);

  if (!SignalFile.is_open())
  {
    cout << "Signal file open error!" << endl;
    cout << argv_1 << endl;
    cout << SignalId << endl;
    exit(0);
  }

  string FileLine;
  vector<float> SignalData;

  while (getline(SignalFile, FileLine))
  {
    float SignalValue;
    convertFromString(SignalValue, FileLine);
    SignalData.push_back(SignalValue);
  }

  ZScoreNormalize(SignalData);
  SignalFile.close();
  
  return SignalData;
}

float gpuCalculate1nDynamicTimeWarping(const vector<vector<float>> &siga,
                                       const vector<vector<float>> &sigb,
                                       vector<vector<float>> &gpudistResult)
{
  // return 0;
  int siga_n = siga.size();
  int sigb_n = sigb.size();
  gpudistResult.resize(siga_n);
  for (int i = 0; i < siga_n; i++)
  {
    gpudistResult[i].resize(sigb_n);
  }
  // int siga_length = 0;
  float *d_distResult = NULL;
  float *d_allColData = NULL;
  float *d_allRowData = NULL;
  unsigned int *d_allRowLength;
  unsigned int *d_allColLength;
  // vector<float *> rowDataList(sigb_n);
  cudaMalloc((void **)&d_distResult, siga_n * sigb_n * sizeof(float));
  CUERR
  cudaMalloc((void **)&d_allColData, siga_n * 1024 * sizeof(float));
  CUERR
  cudaMalloc((void **)&d_allRowData, sigb_n * 1024 * sizeof(float));
  CUERR

  vector<unsigned int> h_allRowLength(sigb_n);
  for (int i = 0; i < sigb_n; i++)
  {
    h_allRowLength[i] = sigb[i].size();
    cudaMemcpy(&d_allRowData[1024 * i], &sigb[i][0],
               h_allRowLength[i] * sizeof(float), cudaMemcpyHostToDevice);
    CUERR
  }
  cudaMalloc((void **)&d_allRowLength, sigb_n * sizeof(unsigned int));
  CUERR
  cudaMemcpy(d_allRowLength, &h_allRowLength[0], sigb_n * sizeof(unsigned int),
             cudaMemcpyHostToDevice);
  CUERR

  vector<unsigned int> h_allColLength(siga_n);
  for (int i = 0; i < siga_n; i++)
  {
    h_allColLength[i] = siga[i].size();
    cudaMemcpy(&d_allColData[1024 * i], &siga[i][0],
               h_allColLength[i] * sizeof(float), cudaMemcpyHostToDevice);
    CUERR
  }
  cudaMalloc((void **)&d_allColLength, siga_n * sizeof(unsigned int));
  CUERR
  cudaMemcpy(d_allColLength, &h_allColLength[0], siga_n * sizeof(unsigned int),
             cudaMemcpyHostToDevice);
  CUERR

  float timesum = 0;
  dim3 threadsPerBlock(1024);
  dim3 blocksPerGrid(sigb_n, siga_n);
  cuDTW_ultimate<<<blocksPerGrid, threadsPerBlock>>>(
      d_allColData, d_allColLength, d_allRowData, d_allRowLength, d_distResult);
  CUERR

  for (int i = 0; i < siga_n; i++)
  {
    cudaMemcpy(&gpudistResult[i][0], &d_distResult[sigb_n * i],
               sigb_n * sizeof(float), cudaMemcpyDeviceToHost);
    CUERR
  }

  cudaFree(d_allColData);
  CUERR
  cudaFree(d_distResult);
  CUERR
  cudaFree(d_allRowData);
  CUERR
  cudaFree(d_allRowLength);
  CUERR

  return timesum;
}

float cpuCalculate1nDynamicTimeWarping(const vector<vector<float>> &siga,
                                       const vector<vector<float>> &sigb,
                                       vector<vector<float>> &cpudistResult)
{

  int siga_n = siga.size();
  int sigb_n = sigb.size();
  cpudistResult.resize(siga_n);
  for (int i = 0; i < siga_n; i++)
  {
    cpudistResult[i].resize(sigb_n);
  }
  omp_set_num_threads(CPUTHREADS);
#pragma omp parallel for
  for (int i = 0; i < siga_n; i++)
  {
    for (int j = 0; j < sigb_n; j++)
    {
      cpudistResult[i][j] = cpuDynamicTimeWarping(siga[i], sigb[j]);
    }
  }

  float timesum = 0;
  // printf("cpu Average time use of DTW= %f sec\n", timesum / sigb_n);
  return timesum;
}

float cpuDynamicTimeWarping(const std::vector<float> &seq1,
                            const std::vector<float> &seq2)
{
  vector<vector<float>> score(seq1.size());

  for (int i = 0; i < seq1.size(); i++)
  {
    score[i].resize(seq2.size());
  }

  for (int i = 0; i < seq1.size(); i++)
  {
    for (int j = 0; j < seq2.size(); j++)
    {
      score[i][j] = std::fabs(seq1[i] - seq2[j]);
    }
  }

  for (int i = 1; i < seq1.size(); i++)
  {
    score[i][0] += score[i - 1][0];
  }

  for (int j = 1; j < seq2.size(); j++)
  {
    score[0][j] += score[0][j - 1];
  }

  for (int i = 1; i < seq1.size(); i++)
  {
    for (int j = 1; j < seq2.size(); j++)
    {
      score[i][j] += std::min(std::min(score[i - 1][j], score[i][j - 1]),
                              score[i - 1][j - 1]);
    }
  }

  float diff = score[seq1.size() - 1][seq2.size() - 1];

  return diff;
}

__global__ void cuDTW_ultimate(float *g_allColData,
                               unsigned int *g_allColLength,
                               float *g_allRowData,
                               unsigned int *g_allRowLength, float *g_odata)
{

  unsigned int inblockThreadIdx = threadIdx.x;
  unsigned int rowLength =
      g_allRowLength[blockIdx.x]; 
  unsigned int colLength = g_allColLength[blockIdx.y];
  float myNum = 0, myColNum;
  __shared__ unsigned int s_turn;
  __shared__ float preNum[1024], prepreNum[1024], rowData[1024];
  rowData[inblockThreadIdx] =
      g_allRowData[blockIdx.x * blockDim.x + threadIdx.x];
  if (inblockThreadIdx == 0)
  {
    s_turn = 0;
  }
  __syncthreads();
  if (inblockThreadIdx < colLength)
  {
    myColNum = g_allColData[blockIdx.y * 1024 + inblockThreadIdx];
    prepreNum[inblockThreadIdx] = preNum[inblockThreadIdx] = 0;
    int col;
    while (s_turn < colLength + rowLength)
    {
      col = s_turn - inblockThreadIdx;
      if (col >= 0 && col < rowLength)
      {
        if (inblockThreadIdx == 0)
        {
          myNum = preNum[inblockThreadIdx] + fabs(myColNum - rowData[col]);
        }
        else
        {
          if (col == 0)
          {
            myNum =
                preNum[inblockThreadIdx - 1] + fabs(myColNum - rowData[col]);
          }
          else
          {
            myNum = min(min(prepreNum[inblockThreadIdx - 1],
                            preNum[inblockThreadIdx - 1]),
                        preNum[inblockThreadIdx]) +
                    fabs(myColNum - rowData[col]);
          }
        }
      }
      __syncthreads();
      prepreNum[inblockThreadIdx] = preNum[inblockThreadIdx];
      preNum[inblockThreadIdx] = myNum;
      if (inblockThreadIdx == 0)
      {
        s_turn++;
      }
      __syncthreads();
    }
  }
  if (inblockThreadIdx == colLength - 1)
  {
    g_odata[blockIdx.y * gridDim.x + blockIdx.x] = myNum;
  }
}

float gpuCalculatemnDynamicTimeWarping_2048(const vector<vector<float>> &siga,
                                            const vector<vector<float>> &sigb,
                                            vector<vector<float>> &gpudistResult)
{
  // return 0;
  int siga_n = siga.size();
  int sigb_n = sigb.size();
  gpudistResult.resize(siga_n);
  for (int i = 0; i < siga_n; i++)
  {
    gpudistResult[i].resize(sigb_n);
  }
  // int siga_length = 0;
  float *d_distResult = NULL;
  float *d_allColData = NULL;
  float *d_allRowData = NULL;
  unsigned int *d_allRowLength;
  unsigned int *d_allColLength;
  // vector<float *> rowDataList(sigb_n);
  cudaMalloc((void **)&d_distResult, siga_n * sigb_n * sizeof(float));
  CUERR
  cudaMalloc((void **)&d_allColData, siga_n * 2048 * sizeof(float));
  CUERR
  cudaMalloc((void **)&d_allRowData, sigb_n * 2048 * sizeof(float));
  CUERR

  
  vector<unsigned int> h_allRowLength(sigb_n);
  for (int i = 0; i < sigb_n; i++)
  {
    h_allRowLength[i] = min(int(sigb[i].size()), 2048);
    cudaMemcpy(&d_allRowData[2048 * i], &sigb[i][0], h_allRowLength[i] * sizeof(float),
               cudaMemcpyHostToDevice);
    CUERR
  }
  cudaMalloc((void **)&d_allRowLength, sigb_n * sizeof(unsigned int));
  CUERR
  cudaMemcpy(d_allRowLength, &h_allRowLength[0], sigb_n * sizeof(unsigned int),
             cudaMemcpyHostToDevice);
  CUERR


  vector<unsigned int> h_allColLength(siga_n);
  for (int i = 0; i < siga_n; i++)
  {
    h_allColLength[i] = min(int(siga[i].size()), 2048);
    cudaMemcpy(&d_allColData[2048 * i], &siga[i][0], h_allColLength[i] * sizeof(float),
               cudaMemcpyHostToDevice);
    CUERR
  }
  cudaMalloc((void **)&d_allColLength, siga_n * sizeof(unsigned int));
  CUERR
  cudaMemcpy(d_allColLength, &h_allColLength[0], siga_n * sizeof(unsigned int),
             cudaMemcpyHostToDevice);
  CUERR

  float timesum = 0;
  dim3 threadsPerBlock(1024);
  dim3 blocksPerGrid(sigb_n, siga_n); 
  cuDTW_2048<<<blocksPerGrid, threadsPerBlock>>>(d_allColData, d_allColLength, d_allRowData,
                                                 d_allRowLength, d_distResult);
  CUERR

  for (int i = 0; i < siga_n; i++)
  {
    cudaMemcpy(&gpudistResult[i][0], &d_distResult[sigb_n * i], sigb_n * sizeof(float),
               cudaMemcpyDeviceToHost);
    CUERR
  }

  cudaFree(d_allColData);
  CUERR
  cudaFree(d_distResult);
  CUERR
  cudaFree(d_allRowData);
  CUERR
  cudaFree(d_allRowLength);
  CUERR

  return timesum;
}

__global__ void cuDTW_2048(float *g_allColData, unsigned int *g_allColLength, float *g_allRowData,
                           unsigned int *g_allRowLength, float *g_odata)
{

  unsigned int rowLength = g_allRowLength[blockIdx.x]; 
  unsigned int colLength = g_allColLength[blockIdx.y];
  float myNum1 = 0, myNum2 = 0, myColNum1, myColNum2;
  __shared__ unsigned int s_turn;
  __shared__ float preNum1[1024], preNum2[1024], prepreNum2[1024], rowData[2048];
  
  rowData[threadIdx.x] = g_allRowData[blockIdx.x * 2048 + threadIdx.x];
  __syncthreads();
  rowData[threadIdx.x + 1024] = g_allRowData[blockIdx.x * 2048 + threadIdx.x + 1024];
  if (threadIdx.x == 0)
  {
    s_turn = 0;
  }
  __syncthreads();
  if (threadIdx.x < (colLength - 1) / 2 + 1)
  {
    
    myColNum1 = g_allColData[blockIdx.y * 2048 + (threadIdx.x) * 2];
    myColNum2 = g_allColData[blockIdx.y * 2048 + (threadIdx.x) * 2 + 1];
    prepreNum2[threadIdx.x] = preNum2[threadIdx.x] = preNum1[threadIdx.x] = 0;
    int col;
    while (s_turn < (colLength - 1) / 2 + 1 + rowLength)
    {                             
      col = s_turn - threadIdx.x; 
      if (col >= 0 && col < rowLength)
      {
        if (threadIdx.x == 0)
        {
          myNum1 = preNum1[0] + fabs(myColNum1 - rowData[col]);

          if (col == 0)
          { 
            myNum2 = myNum1 + fabs(myColNum2 - rowData[col]);
          }
          else
          {
            myNum2 = min(min(myNum1, preNum1[0]), preNum2[0]) +
                     fabs(myColNum2 - rowData[col]);
          }
        }
        else
        {
          if (col == 0)
          {
            myNum1 = preNum2[threadIdx.x - 1] + fabs(myColNum1 - rowData[col]);
            myNum2 = myNum1 + fabs(myColNum2 - rowData[col]);
          }
          else
          {
            myNum1 = min(min(prepreNum2[threadIdx.x - 1], preNum2[threadIdx.x - 1]),
                         preNum1[threadIdx.x]) +
                     fabs(myColNum1 - rowData[col]);
            myNum2 = min(min(myNum1, preNum1[threadIdx.x]), preNum2[threadIdx.x]) +
                     fabs(myColNum2 - rowData[col]);
          }
        }
      }
      __syncthreads();
      prepreNum2[threadIdx.x] = preNum2[threadIdx.x];
      preNum2[threadIdx.x] = myNum2;
      preNum1[threadIdx.x] = myNum1;
      if (threadIdx.x == 0)
      {
        s_turn++;
      }
      __syncthreads();
    }
  }
  if (threadIdx.x == (colLength - 1) / 2)
  {
    if (colLength % 2 == 0)
      g_odata[blockIdx.y * gridDim.x + blockIdx.x] = myNum2;
    else
      g_odata[blockIdx.y * gridDim.x + blockIdx.x] = myNum1;
  }
}
