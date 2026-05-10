#include "main.hpp"
#include "WienerFilter.hpp"

// void test()
// {
//     const int M = 10;
//     string inputPath = "data/in.txt";
//     string outputPath = "data/results/output.txt";

//     cout << endl
//          << string(60, '=') << endl;
//     cout << " *======= TESTCASE WITH INPUT: <" << inputPath << "> =======* \n";

//     WienerFilter wf(M);
//     wf.loadData(inputPath, "data/desired_1.txt");
//     wf.process();
//     wf.saveAndDisplay(outputPath);
// }

void mainTest()
{
    const int M = 10;

    for (int testNumber = 1; testNumber <= 3; testNumber++)
    {
        string inputPath = "data/input_" + to_string(testNumber) + ".txt";
        string outputPath = "data/results/output_" + to_string(testNumber) + ".txt";

        cout << endl
             << string(60, '=') << endl;
        cout << " *======= TESTCASE WITH INPUT: <" << inputPath << "> =======* \n";

        WienerFilter wf(M);
        wf.loadData(inputPath, "data/desired.txt");
        wf.process();
        wf.saveAndDisplay(outputPath);
    }

    cout << "\n>> NOTE: You can recheck 'Your Result' and compare with 'Expected Results' (provided by lecturer) in folder <data/result> \n\n";
}

int main()
{
    mainTest();
}