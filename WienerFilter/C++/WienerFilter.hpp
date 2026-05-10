#pragma once

#include "main.hpp"
#include "WienerFilter.hpp"

class WienerFilter
{
private:
    vector<double> desired_signal;
    vector<double> input_signal;
    vector<double> output_signal;
    vector<double> optimize_coefficient;
    double mmse;

    int M;
    const int MAX_SIZE = 10;
    vector<double> solveLinearSystem(vector<vector<double>> R_M, vector<double> gamma_dx);
    string getRoundDouble(const double &num);
    bool saveData(const string &outputPath, const string &content);

public:
    WienerFilter(int M);

    bool loadData(const string &inputFile, const string &desiredPath);
    void process();
    void saveAndDisplay(const string &outputPath);
};