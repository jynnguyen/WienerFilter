#include "WienerFilter.hpp"

WienerFilter::WienerFilter(int filterLen) : mmse(0.0), M(filterLen)
{
}

bool WienerFilter::loadData(const string &inputFile, const string &desiredFile)
{
    auto loadFunc = [&](const string &path, vector<double> &vec) -> bool
    {
        ifstream file(path);
        if (!file.is_open())
            return false;

        double value;
        while (file >> value)
        {
            vec.push_back(value);
        }
        file.close();
        return true;
    };

    if (!loadFunc(desiredFile, desired_signal) || !loadFunc(inputFile, input_signal))
    {
        cout << "[Load Data]: Could not open files." << endl;
        return false;
    }

    if ((int)desired_signal.size() != MAX_SIZE || (int)input_signal.size() != MAX_SIZE)
    {
        cout << "[Load Data]: size not match." << endl;
        return false;
    }

    return true;
}

bool WienerFilter::saveData(const string &outputPath, const string &content)
{
    ofstream file(outputPath);
    if (file.is_open())
    {
        file << content;
        file.close();
        return true;
    }
    return false;
}

void WienerFilter::process()
{
    int inputSize = input_signal.size();
    vector<double> gamma_xx(M, 0.0);
    vector<double> gamma_dx(M, 0.0);

    // Colleration matrix
    for (int k = 0; k < M; k++)
    {
        for (int n = k; n < inputSize; n++)
        {
            gamma_xx[k] += input_signal[n] * input_signal[n - k];
            gamma_dx[k] += desired_signal[n] * input_signal[n - k];
        }
        gamma_xx[k] /= (inputSize);
        gamma_dx[k] /= (inputSize);
    }

    // Hermitian Toeplitz matrix
    vector<vector<double>> R(M, vector<double>(M));
    for (int l = 0; l < M; l++)
    {
        for (int k = 0; k < M; k++)
            R[l][k] = gamma_xx[abs(l - k)];
    }
    optimize_coefficient = solveLinearSystem(R, gamma_dx);

    output_signal.assign(inputSize, 0.0);
    for (int n = 0; n < inputSize; n++)
    {
        for (int k = 0; k < M; k++)
        {
            if (n - k >= 0)
                output_signal[n] += optimize_coefficient[k] * input_signal[n - k];
        }
    }

    double sumError = 0.0;
    for (int n = 0; n < inputSize; n++)
    {
        double error = desired_signal[n] - output_signal[n];
        sumError += error * error;
    }
    mmse = sumError / inputSize;
}

void WienerFilter::saveAndDisplay(const string &outputPath)
{
    stringstream ss;
    ss << "Filtered output: ";
    for (int i = 0; i < (int)output_signal.size(); ++i)
    {
        ss << getRoundDouble(output_signal[i]);
        if (i != (int)output_signal.size() - 1)
            ss << " ";
    }
    ss << "\nMMSE: " << getRoundDouble(mmse) << endl;

    if (!saveData(outputPath, ss.str()))
        cout << "[Save Data]: could not save\n";

    cout << ss.str();
}

vector<double> WienerFilter::solveLinearSystem(vector<vector<double>> R_M, vector<double> gamma_dx)
{
    int n = gamma_dx.size();

    for (int i = 0; i < n; i++)
    {
        // 3. Thực hiện quá trình khử xuôi (Forward Elimination)
        for (int k = i + 1; k < n; k++)
        {
            double factor = R_M[k][i] / R_M[i][i];
            for (int j = i; j < n; j++)
            {
                R_M[k][j] -= factor * R_M[i][j];
            }
            gamma_dx[k] -= factor * gamma_dx[i];
        }
    }

    // 4. Giải ngược (Back Substitution) để tìm vector hệ số h_opt
    vector<double> h_opt(n, 0);
    for (int i = n - 1; i >= 0; i--)
    {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++)
        {
            sum += R_M[i][j] * h_opt[j];
        }
        h_opt[i] = (gamma_dx[i] - sum) / R_M[i][i];
    }

    return h_opt;
}

string WienerFilter::getRoundDouble(const double &num)
{
    double newNum = num;
    stringstream ss;
    ss << fixed << setprecision(1) << newNum;
    return ss.str();
}