#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <limits>

using namespace std;

const int INF = 0x3f3f3f3f;
const double eps = 1e-6;

vector<vector<double>> tableau;
vector<vector<double>> tableauAux;
vector<double> viableSolution;
vector<double> certificate;
vector<vector<int>> A;
vector<int> b, c, B;
int n = 0, m = 0;

void print(vector<vector<double>> &t)
{
    cout << endl
         << "tableau " << endl;
    for (int i = 0; i < t.size(); i++)
    {
        for (int j = 0; j < t[0].size(); j++)
            cout << t[i][j] << ' ';
        cout << endl;
    }
    cout << endl
         << endl;
}

void readData(char *argv[])
{
    string filePath = argv[1];
    ifstream inputFile(filePath);

    inputFile >> n >> m;
    tableau.resize(n + 1, vector<double>(m + n + 1, 0));
    tableauAux.resize(n + 1, vector<double>(m + (2 * n) + 1, 0));
    A.resize(n, vector<int>(m, 0));
    b.resize(n, 0);
    c.resize(m, 0);
    B.resize(n, 0);
    viableSolution.resize(m, 0);
    certificate.resize(n, 0);

    for (int i = 0; i < m; ++i)
        inputFile >>
            c[i];

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
            inputFile >> A[i][j];
        inputFile >> b[i];
    }
}

void initTableau()
{
    for (int i = 1; i <= n; i++)
    {
        for (int j = 0; j < n + m; j++)
            tableau[i][j] = tableauAux[i][j];
        tableau[i][m + n] = tableauAux[i][m + (2 * n)];
    }
    for (int i = 0; i < m; i++)
        tableau[0][i + n] = -c[i];
}

void initTableauAux()
{
    for (int i = 0; i < n; i++)
    {
        tableauAux[i + 1][i] = 1;
        tableauAux[0][i + m + n] = 1;
        tableauAux[i + 1][i + m + n] = 1;
        B[i] = i + m + n;
        for (int j = 0; j < m; j++)
            tableauAux[i + 1][j + n] = A[i][j];
        tableauAux[i + 1][m + 2 * n] = b[i];
        if (b[i] < 0)
        {
            for (int j = 0; j <= m + (2 * n); j++)
            {
                if (tableauAux[i + 1][j] != 0)
                    tableauAux[i + 1][j] *= -1;
            }
        }
    }

    for (int i = 1; i <= n; i++)
    {
        bool op = tableauAux[i][B[i - 1]] > 0;
        for (int j = 0; j <= m + (2 * n); j++)
        {
            if (op)
                tableauAux[0][j] -= tableauAux[i][j];
            else
                tableauAux[0][j] += tableauAux[i][j];
        }
    }
}

int findPivot(vector<vector<double>> &t, int col)
{
    int last = t[0].size() - 1;

    if (n == 1)
        return ((t[1][col] != 0) && ((t[1][last] / t[1][col]) >= 0)) ? 1 : -1;
    else
    {
        int l = -1;
        double bigger = numeric_limits<double>::max();
        for (int i = 1; i <= n; i++)
        {
            if (t[i][col] > 0)
            {
                double factor = t[i][last] / t[i][col];
                if (factor < bigger)
                {
                    bigger = factor;
                    l = i;
                }
            }
        }
        return l;
    }
}

void pivot(vector<vector<double>> &t, int i, int j)
{
    print(t);
    cout << i << ' ' << j << endl;
    double pivotValue = t[i][j];
    cout << "pivotValue " << pivotValue << endl
         << endl;
    if (pivotValue > 0 || pivotValue < 0)
    {
        int cols = t[0].size();

        for (int k = 0; k < cols; k++)
            t[i][k] /= pivotValue;

        for (int k = 0; k <= n; ++k)
        {
            if (k != i)
            {
                double factor = t[k][j];
                for (int c = 0; c < cols; c++)
                    t[k][c] -= factor * t[i][c];
            }
        }
    }
}

void canonical(vector<vector<double>> &t)
{
    for (int i = 0; i < n; i++)
        pivot(t, i + 1, B[i]);
}

int getNegativeColumn(vector<vector<double>> &t)
{
    int index = -1;
    int value = 0;

    for (int i = n; i < n + m; i++)
    {
        if (t[0][i] < value)
        {
            value = t[0][i];
            index = i;
        }
    }

    return index;
}

void generateViableSolution()
{
    for (int i = 0; i < n; i++)
        if (B[i] < m + n)
            viableSolution[B[i] - n] = tableau[i + 1][m + n];
}

void generateCertificate(bool isInf, int col)
{
    if (isInf)
    {
        certificate.resize(m, 0);
        certificate[col - n] = 1;
        for (int i = 0; i < n; i++)
            if (B[i] < m + n)
                certificate[B[i] - n] = -tableau[i + 1][col];
    }
    else
        for (int i = 0; i < n; i++)
            certificate[i] = tableau[0][i];
}

void optimalSolution()
{
    cout << "otima" << endl;
    cout << tableau[0][m + n] << endl;

    generateViableSolution();
    for (int i = 0; i < m; i++)
        cout << viableSolution[i] << ' ';
    cout << endl;

    generateCertificate(false, 0);
    for (int i = 0; i < n; i++)
        cout << certificate[i] << ' ';
    cout << endl;
}

void ilimitedSolution(int col)
{
    cout << "ilimitada" << endl;

    generateViableSolution();
    for (int i = 0; i < m; i++)
        cout << viableSolution[i] << ' ';
    cout << endl;

    generateCertificate(true, col);
    for (int i = 0; i < m; i++)
        cout << certificate[i] << ' ';
    cout << endl;
}

void unviableSolution()
{
    cout << "inviavel" << endl;

    // generateCertificate(false, 0);
    // for (int i = 0; i < n; i++)
    //     cout << certificate[i] << ' ';
    // cout << endl;
}

int objective(vector<vector<double>> &t)
{
    canonical(t);

    int col = getNegativeColumn(t);
    while (col > 0)
    {
        int l = findPivot(t, col);

        if (l > 0)
            B[l - 1] = col;
        else
        {
            ilimitedSolution(col);
            return INF;
        }

        canonical(t);
        col = getNegativeColumn(t);
    }

    return t[0][t.size() - 1];
}

int main(int argc, char *argv[])
{
    cout << fixed << setprecision(0);

    readData(argv);
    if (n <= 0 || m <= 0)
    {
        cout << "ERRO: matriz com dimensões menores ou iguais a zero" << endl;
        return 0;
    }

    initTableauAux();

    double sol = objective(tableauAux);
    if (sol < -eps)
    {
        unviableSolution();
        return 0;
    }
    else if (sol == INF)
        return 0;

    initTableau();

    sol = objective(tableau);
    if (sol == INF)
        return 0;
    else
        optimalSolution();

    return 0;
}
