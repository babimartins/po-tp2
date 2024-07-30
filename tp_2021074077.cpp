#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>

using namespace std;

const int INF = 0x3f3f3f3f;
const double eps = 1e-6;

vector<vector<double>> tableau;
vector<vector<double>> tableauAux;
vector<double> viableSolution;
vector<double> certificate;
vector<vector<int>> A;
vector<int> b, c, B;
int n, m, N, M;

void initTableauAux()
{
    vector<int> negIndexes;

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
            negIndexes.push_back(i);
    }

    for (int i = 0; i < negIndexes.size(); i++)
    {
        for (int j = 0; j < M; j++)
        {
            if (tableauAux[negIndexes[i] + 1][j] != 0)
                tableauAux[negIndexes[i] + 1][j] *= -1;
        }
    }

    for (int i = 1; i < N; i++)
    {
        bool op = tableauAux[i][B[i - 1]] > 0;
        for (int j = 0; j < M; j++)
        {
            if (op)
                tableauAux[0][j] -= tableauAux[i][j];
            else
                tableauAux[0][j] += tableauAux[i][j];
        }
    }

    // for (int i = 0; i < s.size(); i++)
    //     s[i] = tableauAux[0][i + n];
}

int findPivot(int col)
{
    int l = 1;
    if (N == 1)
        return ((tableauAux[l][col] != 0) && ((tableauAux[l][M - 1] / tableauAux[l][col]) >= 0)) ? l : -1;
    else
    {
        for (int i = 2; i < N; i++)
        {
            // cout << "tableauAux[" << i << "][" << col << "] " << tableauAux[i][col] << endl;
            if (tableauAux[i][col] > 0)
            {
                if (tableauAux[l][col] <= 0 || ((tableauAux[i][M - 1] / tableauAux[i][col]) < (tableauAux[l][M - 1] / tableauAux[l][col])))
                    l = i;
            }
        }
        return ((tableauAux[l][M - 1] != 0) && ((tableauAux[l][M - 1] / tableauAux[l][col]) >= 0)) ? l : -1;
    }
}

void gauss(int i, int j)
{
    double pivotValue = tableauAux[i][j];
    if (pivotValue == 0)
        return;

    for (int col = 0; col < M; ++col)
        tableauAux[i][col] /= pivotValue;

    for (int k = 0; k < N; ++k)
    {
        if (k != i)
        {
            double factor = tableauAux[k][j];
            for (int col = 0; col < M; ++col)
                tableauAux[k][col] -= factor * tableauAux[i][col];
        }
    }
}

void pivot(int cols)
{
    for (int col = n; col < (cols + n); col++)
    {
        // cout << "col " << col << endl;
        if (tableauAux[0][col] < 0)
        {
            int l = findPivot(col);
            // cout << "pivot " << l << ' ' << col << endl;
            if (l > 0)
            {
                B[l - 1] = col;
                gauss(l, col);
            }
        }
    }
}

void initTableau()
{
    for (int i = 0; i < n; i++)
    {
        tableau[i + 1][i] = 1;
        for (int j = 0; j < m; j++)
            tableau[i + 1][j + n] = tableauAux[i + 1][j + n];
        tableau[i + 1][m + n] = tableauAux[i + 1][M - 1]; // TODO: inclui o zero final da primeira linha ou não? agora está não
    }
    for (int i = 0; i < m; i++)
        tableau[0][i + n] = -c[i];
}

int main(int argc, char *argv[])
{
    cout << setprecision(3);

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
    N = tableauAux.size();
    M = tableauAux[0].size();

    for (int i = 0; i < m; ++i)
        inputFile >>
            c[i];

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
            inputFile >> A[i][j];
        inputFile >> b[i];
    }

    initTableauAux();

    pivot(m + n);
    // TODO: se B não mudou, então é inviavel
    // mas ainda é preciso continuar e conferir novamente, pode ser inviavel em outros passos

    /*
        TODO: identificar se ela é inviável quando c é diferente de zero acima da identidade
        for (int i = 0; i < n; i++)
        {
            if (tableauAux[0][i + m + n] != 0)
            {
                cout << i + m + n << endl;
                cout << tableauAux[0][i + m + n] << endl;
                cout << "inviável" << endl;
                // return 0;
            }
        }
    */

    int aux = 0;
    for (int i = 0; i < m; i++)
    {
        if ((B[aux] - n) != i)
            viableSolution[i] = 0;
        else
        {
            viableSolution[i] = b[aux];
            aux++;
        }
    }

    cout << "tableauAux depois de pivotear" << endl;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
            cout << tableauAux[i][j] << ' ';
        cout << endl;
    }
    cout << endl;

    initTableau();

    cout << "tableau new" << endl;
    for (int i = 0; i < tableau.size(); i++)
    {
        for (int j = 0; j < tableau[0].size(); j++)
            cout << tableau[i][j] << ' ';
        cout << endl;
    }
    cout << endl;

    // cout << "B" << endl;
    // for (int i = 0; i < n; i++)
    // {
    //     cout << B[i] << ' ';
    // }
    // cout << endl;

    // cout << "viable solution" << endl;
    // for (int i = 0; i < m; i++)
    // {
    //     cout << viableSolution[i] << ' ';
    // }
    // cout << endl;

    return 0;
}
